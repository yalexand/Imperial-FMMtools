 
classdef FMMtools_data_controller < handle 
    
    % Copyright (C) 2013 Imperial College London.
    % All rights reserved.
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along
    % with this program; if not, write to the Free Software Foundation, Inc.,
    % 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    %
    % This software tool was developed with support from the UK 
    % Engineering and Physical Sciences Council 
   
    properties(Constant)
        data_settings_filename = 'FMMtools_data_settings.xml';
    end
    
    properties(SetObservable = true)
            
        ADC_trails_features_data = []; 
        ADC_feature_names = [];
        % feature vector
        ADC_fv_all = {'trail_length', 'entropy','energy','Ea','R1','R2','R3','Ed1','Ed2','Ed3','Ed4','Ed5','Ed6'};
        ADC_fv_selected = {'entropy','energy'};
        % available "conditions", i.e. types of motion
        groups_all = {'breathe','general','head','limb','startle','other','T1','T2','T3','T4'};
        groups_available = {'breathe','general','head','limb','startle','other','T1','T2','T3','T4'};
        groups_selected = {'breathe','general','startle'};
        %
        supervised_learning_method = {'Linear','Quadratic','kNN','PCA->Linear','PCA->Quadratic','PCA->kNN'};
                                
        corr_map_W = 80;
        
    end                    
    
    properties(Transient)
        
        DefaultDirectory = ['C:' filesep];
        RootDirectory = ['C:' filesep]; 
        ExcelDirectory = [];
                
        current_filename = [];
        
        current_data = [];
        
        current_IMU_OK = [];
        current_ADC_OK = [];
        
        current_ADC_pre_processed = [];        
        current_ADC_segmented = [];                
        current_ADC_PSD = [];                
        current_ADC_preprocessed_PSD = [];         
        
        current_IMU_pre_processed = [];        
        current_IMU_segmented = [];          
        current_IMU_PSD = []; 
        current_IMU_preprocessed_PSD = []; 
        
        current_omega_ADC; % frequencies
        current_omega_IMU;

        current_annotation = [];
        current_annotation_time = [];
        
        %
        subj_data = [];
        subj_filenames = [];
        %
                
        %
        Fs_ADC = 1024; %Hz
        Fs_IMU = 64; %Hz

        preprocessing_types = {'freq.notch-comb + median','none'};
        ADC_prprss_Median_Size = 3;
        ADC_prprss_notch_comb_fb = 39.1328; % base frequency, Hz
        ADC_prprss_notch_comb_bw = 0.0022; % bandwidth        
        %
        segmentation_types = {'moving average subtraction + thresholding','ICA assisted thresholding','annotator"s b/g/s'};
        ADC_segm_Moving_Average_Window = 5; % seconds!
        ADC_segm_Minimal_Trail_Duration = 1.2;  %seconds       
        %
        supervised_learning_types = {'annotator"s + segmentation', ...
                                            'annotator"s b/g/s', ...
                                            'auto annotated'};
        %
        exclude_IMU = true;
        
        IMU_sgm_low_signal_quantile = 0.1;
        IMU_sgm_signal_cutoff = 0.05; % weaker signals rejected
        IMU_sgm_SN = 200; % signal-to-noise
        %
        ADC_sgm_low_signal_quantile = 0.1;
        ADC_sgm_signal_cutoff = 0.004; % weaker signals rejected
        ADC_sgm_SN = 150; % signal-to-noise
        %
        annotators_delay = 2; % seconds
        annotators_reaction = 0.4; % seconds
        %
        SGM_MODE = 'OR'; % 'AND2' is for AND with threshold on 2 detections
                        
    end    
        
    properties(Transient,Hidden)
        % Properties that won't be saved to a data_settings_file etc.
                
    end
    
    events
        
    end
            
    methods
        
        function obj = FMMtools_data_controller(varargin)            
            %   
            handles = args2struct(varargin);
            assign_handles(obj,handles);            
                                
            try 
                obj.load_settings([pwd filesep obj.data_settings_filename]);
            catch            
            end
                    
            obj.RootDirectory = pwd;   
            
           if isempty(obj.ExcelDirectory)
                hw = waitbar(0,'looking for Excel directory..');
                waitbar(0.1,hw);                
                if ispc
                       prevdir = pwd;
                       cd('c:\');
                       [~,b] = dos('dir /s /b excel.exe');
                       if ~strcmp(b,'File Not Found')
                            filenames = textscan(b,'%s','delimiter',char(10));
                            s = char(filenames{1});
                            s = s(1,:);
                            s = strsplit(s,'excel.exe');
                            obj.ExcelDirectory = s{1};
                       end
                       cd(prevdir);
                elseif ismac
                    % to do
                else
                    % to do
                end                
                delete(hw); drawnow;                               
            end            
                        
        end
                
%-------------------------------------------------------------------------%
        function delete(obj)
            obj.save_settings([obj.RootDirectory filesep obj.data_settings_filename]);
        end       
%-------------------------------------------------------------------------%                
         function clear_current_data(obj,~)
                         
            obj.current_filename = [];

            obj.current_data = [];

            obj.current_IMU_OK = [];
            obj.current_ADC_OK = [];

            obj.current_ADC_pre_processed = [];        
            obj.current_ADC_segmented = [];                
            obj.current_ADC_PSD = [];                
            obj.current_ADC_preprocessed_PSD = [];         

            obj.current_IMU_pre_processed = [];        
            obj.current_IMU_segmented = [];          
            obj.current_IMU_PSD = []; 
            obj.current_IMU_preprocessed_PSD = []; 
            
            obj.current_annotation = [];
            obj.current_annotation_time = [];
            
            obj.current_omega_ADC = [];
            obj.current_omega_IMU = [];
                        
         end
%-------------------------------------------------------------------------%        
        function res = load_single_subject(obj,prp_type,sgm_type,~)
            
            res = false;
                        
            [filename,pathname] = uigetfile({'*.mat','Subject Records Files'},'Select data file',obj.DefaultDirectory);
            if filename == 0, return, end;       

            obj.clear_current_data;
            obj.subj_data = [];
            obj.subj_filenames = [];            
            
            obj.ADC_trails_features_data = [];
            obj.open_datafile([pathname filesep filename]);
                 
            obj.current_filename = filename;                    
            obj.DefaultDirectory = pathname;
            %
            obj.pre_process_ADC(prp_type);
            obj.segment_ADC(sgm_type);
            obj.pre_process_IMU([]);
            obj.segment_IMU([]);
            obj.calculate_PSDs;
            %
            obj.possibly_exclude_ADC_findings_with_simultaneous_IMU_response;
            %                            
            %
            % fill 1-st element in data struct
                    subj_elem.filename = obj.current_filename;                                
                    subj_elem.data = obj.current_data;
                    %
                    subj_elem.IMU_OK = obj.current_IMU_OK;
                    subj_elem.ADC_OK = obj.current_ADC_OK;
                    %
                    subj_elem.ADC_pre_processed = obj.current_ADC_pre_processed;        
                    subj_elem.ADC_segmented = obj.current_ADC_segmented;                
                    subj_elem.ADC_PSD = obj.current_ADC_PSD;                
                    subj_elem.ADC_preprocessed_PSD = obj.current_ADC_preprocessed_PSD;
                    %
                    subj_elem.IMU_pre_processed = obj.current_IMU_pre_processed;        
                    subj_elem.IMU_segmented = obj.current_IMU_segmented;
                    subj_elem.IMU_PSD = obj.current_IMU_PSD;
                    subj_elem.IMU_preprocessed_PSD = obj.current_IMU_preprocessed_PSD;
                    %
                    subj_elem.annotation  = obj.current_annotation;
                    subj_elem.annotation_time = obj.current_annotation_time;
                    %
                    subj_elem.omega_ADC = obj.current_omega_ADC;
                    subj_elem.omega_IMU = obj.current_omega_IMU;                    
                    %
                    obj.subj_data = [obj.subj_data; subj_elem];
                    obj.subj_filenames = cellstr(filename);
                    %
                    % set up available groups
                    available_types_ind = unique(obj.current_annotation);
                    obj.groups_available = cell(1,length(available_types_ind));
                    for k=1:numel(available_types_ind)
                        obj.groups_available(k) = obj.groups_all(available_types_ind(k));
                    end 
                    
                    res = true;
        end
%-------------------------------------------------------------------------%        
        function res = load_multiple_subjects(obj,prp_type,sgm_type,verbose,~)
                     
        res = false;
            
        % get list of files...
        [filenames, pathname] = uigetfile('*.mat','Select data files',obj.DefaultDirectory,'MultiSelect','on');
        if pathname == 0, return, end;
                     
        filenames = sort_nat(cellstr(filenames));

        hw = [];
        waitmsg = 'Loading ...';
        if verbose
            hw = waitbar(0,waitmsg);
        end 
                                            
        obj.subj_data = [];
        
        all_annotations = [];
        
        for k = 1:numel(filenames)                                                            
            if ~isempty(hw), waitbar(k/numel(filenames),hw); drawnow, end;
                obj.clear_current_data;                   
                obj.current_filename = char(filenames{k});                
                obj.open_datafile([pathname filesep obj.current_filename]);
                obj.pre_process_ADC(prp_type);
                obj.segment_ADC(sgm_type);
                obj.pre_process_IMU([]); 
                obj.segment_IMU([]);
                obj.calculate_PSDs;                                            
                %
                obj.possibly_exclude_ADC_findings_with_simultaneous_IMU_response;
                %                
                    %
                    subj_elem.filename = obj.current_filename;                                
                    subj_elem.data = obj.current_data;
                    %
                    subj_elem.IMU_OK = obj.current_IMU_OK;
                    subj_elem.ADC_OK = obj.current_ADC_OK;
                    %
                    subj_elem.ADC_pre_processed = obj.current_ADC_pre_processed;        
                    subj_elem.ADC_segmented = obj.current_ADC_segmented;                
                    subj_elem.ADC_PSD = obj.current_ADC_PSD;                
                    subj_elem.ADC_preprocessed_PSD = obj.current_ADC_preprocessed_PSD;
                    %
                    subj_elem.IMU_pre_processed = obj.current_IMU_pre_processed;        
                    subj_elem.IMU_segmented = obj.current_IMU_segmented;
                    subj_elem.IMU_PSD = obj.current_IMU_PSD;
                    subj_elem.IMU_preprocessed_PSD = obj.current_IMU_preprocessed_PSD;
                    %
                    subj_elem.annotation  = obj.current_annotation;
                    subj_elem.annotation_time = obj.current_annotation_time;
                    %
                    subj_elem.omega_ADC = obj.current_omega_ADC;
                    subj_elem.omega_IMU = obj.current_omega_IMU;                    
                    %
                obj.subj_data = [obj.subj_data; subj_elem];
                
                all_annotations = [all_annotations; obj.current_annotation];
        end
        if ~isempty(hw), delete(hw), drawnow; end;
                
                    % set up available groups
                    available_types_ind = unique(all_annotations);
                    obj.groups_available = cell(1,length(available_types_ind));
                    for k=1:numel(available_types_ind)
                        obj.groups_available(k) = obj.groups_all(available_types_ind(k));
                    end
                
        obj.subj_filenames = filenames;
        
        obj.DefaultDirectory = pathname;
        obj.ADC_trails_features_data = [];
        
        res = true;

        end     
%-------------------------------------------------------------------------%        
function switch_current_to_subject(obj,data_filename,~)
    
    [~,ind]=ind2sub(size(obj.subj_filenames),strmatch(data_filename,obj.subj_filenames,'exact'));
    
    subj_elem = obj.subj_data(ind);
    
    obj.current_filename = subj_elem.filename;                                
    obj.current_data = subj_elem.data;
    %
    obj.current_IMU_OK = subj_elem.IMU_OK;
    obj.current_ADC_OK = subj_elem.ADC_OK;
    %
    obj.current_ADC_pre_processed = subj_elem.ADC_pre_processed;        
    obj.current_ADC_segmented = subj_elem.ADC_segmented;                
    obj.current_ADC_PSD = subj_elem.ADC_PSD;                
    obj.current_ADC_preprocessed_PSD = subj_elem.ADC_preprocessed_PSD;
    %
    obj.current_IMU_pre_processed = subj_elem.IMU_pre_processed;        
    obj.current_IMU_segmented = subj_elem.IMU_segmented;
    obj.current_IMU_PSD = subj_elem.IMU_PSD;
    obj.current_IMU_preprocessed_PSD = subj_elem.IMU_preprocessed_PSD;
    %
    obj.current_annotation  = subj_elem.annotation;
    obj.current_annotation_time = subj_elem.annotation_time;
    %
    obj.current_omega_ADC = subj_elem.omega_ADC;
    obj.current_omega_IMU = subj_elem.omega_IMU;
        
end
        
%-------------------------------------------------------------------------%
        function open_datafile(obj,full_filename,~)
            %
            % to do
            %     
            load(full_filename);
            [~, fname, ~] = fileparts(full_filename);
            subjname = char(strrep(fname,'.mat',''));
            obj.current_data = eval(subjname);
        end       
%-------------------------------------------------------------------------%                    
        function save_settings(obj,fname,~)        
                settings = [];
                settings.DefaultDirectory = obj.DefaultDirectory;
                settings.Fs_ADC = obj.Fs_ADC;
                settings.Fs_IMU = obj.Fs_IMU;
                %
                settings.ADC_prprss_Median_Size = obj.ADC_prprss_Median_Size;
                settings.ADC_prprss_notch_comb_fb = obj.ADC_prprss_notch_comb_fb;
                settings.ADC_prprss_notch_comb_bw = obj.ADC_prprss_notch_comb_bw;
                %
                settings.ADC_segm_Moving_Average_Window = obj.ADC_segm_Moving_Average_Window;
                settings.ADC_segm_Minimal_Trail_Duration = obj.ADC_segm_Minimal_Trail_Duration; 
                %
                settings.ADC_fv_selected = obj.ADC_fv_selected;
                settings.exclude_IMU = obj.exclude_IMU;
                %
                settings.IMU_sgm_low_signal_quantile = obj.IMU_sgm_low_signal_quantile;
                settings.IMU_sgm_signal_cutoff = obj.IMU_sgm_signal_cutoff;
                settings.IMU_sgm_SN = obj.IMU_sgm_SN;
                %
                settings.ADC_sgm_low_signal_quantile = obj.ADC_sgm_low_signal_quantile;
                settings.ADC_sgm_signal_cutoff = obj.ADC_sgm_signal_cutoff;
                settings.ADC_sgm_SN = obj.ADC_sgm_SN;                
                %
                settings.annotators_delay = obj.annotators_delay;
                settings.annotators_reaction = obj.annotators_reaction;                
                %
                settings.corr_map_W = obj.corr_map_W; 
                %
                settings.ExcelDirectory = obj.ExcelDirectory;
                settings.SGM_MODE = obj.SGM_MODE;
            try
                xml_write(fname,settings);
            catch
                disp('xml_write: error, settings were not saved');
            end
        end
%-------------------------------------------------------------------------%
        function load_settings(obj,fname,~)
             if exist(fname,'file') 
                [ settings, ~ ] = xml_read (fname);                                 
                obj.DefaultDirectory = settings.DefaultDirectory;                                                                  
                obj.Fs_ADC = settings.Fs_ADC;
                obj.Fs_IMU = settings.Fs_IMU;                
                %
                obj.ADC_prprss_Median_Size = settings.ADC_prprss_Median_Size;
                obj.ADC_prprss_notch_comb_fb = settings.ADC_prprss_notch_comb_fb;
                obj.ADC_prprss_notch_comb_bw = settings.ADC_prprss_notch_comb_bw;
                %
                obj.ADC_segm_Moving_Average_Window = settings.ADC_segm_Moving_Average_Window;
                obj.ADC_segm_Minimal_Trail_Duration = settings.ADC_segm_Minimal_Trail_Duration;
                %
                obj.ADC_fv_selected = settings.ADC_fv_selected;
                obj.exclude_IMU = settings.exclude_IMU;                
                %
                obj.IMU_sgm_low_signal_quantile = settings.IMU_sgm_low_signal_quantile;
                obj.IMU_sgm_signal_cutoff = settings.IMU_sgm_signal_cutoff;
                obj.IMU_sgm_SN = settings.IMU_sgm_SN;
                %
                obj.ADC_sgm_low_signal_quantile = settings.ADC_sgm_low_signal_quantile;
                obj.ADC_sgm_signal_cutoff = settings.ADC_sgm_signal_cutoff;
                obj.ADC_sgm_SN = settings.ADC_sgm_SN;                                
                %
                obj.annotators_delay = settings.annotators_delay;
                obj.annotators_reaction = settings.annotators_reaction;
                %
                obj.corr_map_W = settings.corr_map_W; 
                %
                obj.ExcelDirectory = settings.ExcelDirectory;
                obj.SGM_MODE = settings.SGM_MODE;
             end
        end
%-------------------------------------------------------------------------%
        function pre_process_ADC(obj,type,~)
               
            if isempty(obj.current_data), return, end;

            obj.current_ADC_pre_processed = zeros(size(obj.current_data.ADC));

            num_ADC_channels = size(obj.current_data.ADC,2);
            hw = waitbar(0,[type ' preprocessing ADC - please wait']);
            for k = 1 : num_ADC_channels
                if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                %
                if strcmp(type,'none')
                    obj.current_ADC_pre_processed(:,k) = obj.current_data.ADC(:,k);
                elseif strcmp(type,'freq.notch-comb + median')
                    %
                    s = obj.current_data.ADC(:,k);
                    %
                    fb = obj.ADC_prprss_notch_comb_fb;
                    bw = obj.ADC_prprss_notch_comb_bw;
                    %
                    for m=1:floor((obj.Fs_ADC/2)/fb)
                        fo = fb*m;
                        wo = fo/(obj.Fs_ADC/2);
                        [b,a] = iirnotch(wo,bw);
                        s = filtfilt(b,a,s);
                    end
                    %
                    if obj.ADC_prprss_Median_Size >= 1
                        mR = fix(obj.ADC_prprss_Median_Size);
                        s = medfilt2(s,[mR 1]);
                    end
                    %
                    obj.current_ADC_pre_processed(:,k) = s;
                end
            end
            if ~isempty(hw), delete(hw), drawnow; end;
            %                        
        end
%-------------------------------------------------------------------------%        
        function segment_ADC(obj,type,~)
            
            if isempty(obj.current_ADC_pre_processed), return, end;

           % extract annotation
            US_data = [];
            %curate
            for k = 1:size(obj.current_data.US,1)
                rec = obj.current_data.US(k,:);
                if ~isempty(cell2mat(rec(1)))
                    US_data = [US_data; rec];
                end
            end

            event_times = US_data(:,1);
            
            % unfortunately "types" index can be 2 or 3..
            event_type_ind = 3;
            tokens = unique(US_data(:,2)); % check if it is 2
            %
            if  ~isempty(intersect(tokens,{'b','g','h','l','s'})) || ...
                ~isempty(intersect(tokens,{'B','G','H','L','S'}))
               event_type_ind = 2;
            end
            
% % % %             event_types = US_data(:,event_type_ind); % 1,2,3,4,5 -> b,g,h,l,s            
% % % %                         
% % % %             obj.current_annotation = zeros(size(event_times,1),1);
% % % %             obj.current_annotation_time = zeros(size(event_times,1),1);            
% % % %             %
% % % %             type_ind = 0; % for now...
% % % %             %
% % % % cnt = zeros(1,6);
% % % %             %
% % % %             for k = 1:length(event_times)
% % % %                 obj.current_annotation_time(k,1) = cell2mat(US_data(k,1));
% % % %                 switch char(event_types(k,1))
% % % %                     case {'b' 'B' 'bb' 'bg'}
% % % %                         type_ind = 1;                        
% % % %                     case {'g' 'G' 'gg' 'gs'}
% % % %                         type_ind = 2;
% % % %                     case {'h' 'H'}
% % % %                         type_ind = 2; % sic!
% % % %                     case {'l' 'L'}
% % % %                         type_ind = 2; % sic!
% % % %                     case {'s' 'S' 'sg'}
% % % %                         type_ind = 5;
% % % %                     otherwise
% % % %                         type_ind = 6; % ehm..
% % % %                         disp(char(event_types(k,1)));
% % % %                 end
% % % %                 obj.current_annotation(k,1) = type_ind;
% % % % cnt(type_ind) = cnt(type_ind)+1;
% % % %             end            
% % % %             %

            % supposed to work only with b,g,s - starts
            event_types = US_data(:,event_type_ind);
            event_times = cell2mat(US_data(:,1));
            %            
            obj.current_annotation = [];
            obj.current_annotation_time = [];
            m_event_times = [];
            m_event_types = [];
            %
            for k = 1:length(event_times)
                type_ind = 0;
                switch char(event_types(k,1))
                    case {'b' 'B' 'bb' 'bg'}
                        type_ind = 1;                        
                    case {'g' 'G' 'gg' 'gs'}
                        type_ind = 2;
                    case {'s' 'S' 'sg'}
                        type_ind = 5;
                    case {'m' 'M'}
                        type_ind = 6;
                    case {'moving' 'MOVING'}                        
                        type_ind = 7;
                    case {'move' 'MOVE'}                        
                        type_ind = 8;                        
                end
                if ~isempty(intersect(type_ind,[1 2 5]))
                    obj.current_annotation = [obj.current_annotation; type_ind];
                    obj.current_annotation_time = [obj.current_annotation_time; event_times(k)];
                elseif type_ind >= 6 % move
                    m_event_times = [m_event_times; event_times(k)];
                    m_event_types = [m_event_types; type_ind];
                end
            end            
            % supposed to work only with b,g,s - end            
            %            
            % extract annotation - ends   

            obj.current_ADC_segmented = zeros(size(obj.current_data.ADC));            
            
            % CREATING M-EVENTS EXCLUSION MASK - START
            %
            % HACKS            
            if strcmp(obj.current_filename,'SUBJECT_009.mat')
                % remove first one - unpaired
                m_event_times = m_event_times(2:numel(m_event_times));
            end
            if strcmp(obj.current_filename,'SUBJECT_012.mat')
                % remove last one - easy fix
                m_event_times = m_event_times(1:numel(m_event_times)-1);
            end            
            % HACKS            
            %
            LMAX = size(obj.current_ADC_segmented,1);
            m_excl_mask = zeros(LMAX,1);
            %
            N_m_events = length(m_event_times);
            if 0 == rem(N_m_events,2) % only if even
                for k=1:2:N_m_events
                    T1=m_event_times(k);
                    T2=m_event_times(k+1);
                    if (T2-T1) < 120 % seconds - this is to treat subj. 8 :)
                        L1 = max(1,round(T1*obj.Fs_ADC));
                        L2 = min(LMAX,round(T2*obj.Fs_ADC));
                        m_excl_mask(L1:L2)=1;
                    end
                end
            end
            % CREATING M-EVENTS EXCLUSION MASK - ENDS
                                          
%debug_h = figure();
            num_ADC_channels = size(obj.current_data.ADC,2);

            hw = waitbar(0,[type ' segmenting ADC - please wait']);            
            if strcmp(type,'moving average subtraction + thresholding') || strcmp(type,'annotator"s b/g/s')            

                for k = 1 : num_ADC_channels
                    if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                    s = obj.current_ADC_pre_processed(:,k);
                        %
                        avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_ADC);
                        %
                        [s,~] = TD_high_pass_filter( s, avr_window );
                        %                    
                        s_ = sqrt(s.*s);
                        %                      
                        t = quantile(s_(:),obj.ADC_sgm_low_signal_quantile); %  take low 10% of signal
                        z = s_(s_<t);
                        %
                            t = obj.ADC_sgm_SN*median(z(:)); % threshold at 100X (or whatever) average noise level
                        %
                        if t < obj.ADC_sgm_signal_cutoff, t = Inf; end; % PRECAUTION AGAINST TOO NOISY SIGNALS                    
                        %
                        z = (s_ > t);
                        %                    
%figure(debug_h);subplot(2,4,k);plot(1:length(s_),s_,'k.-',1:length(s_),t*ones(1,length(s_)),'r-');grid on;xlabel(num2str(t));                    
                        %
                        if 0~=sum(z)
                            min_size = obj.ADC_segm_Minimal_Trail_Duration*obj.Fs_ADC;
                            %
                            % exclude m-events mask - before dilating
                            z = z & (~m_excl_mask);
                            %
                            SE = strel('line',min_size,90);
                            z = imdilate(z,SE);                            
                        end
                        obj.current_ADC_segmented(:,k) = z;
                end
                                               
            elseif strcmp(type,'ICA assisted thresholding') 
                mixedsig = [];        
                for k = 1 : num_ADC_channels
                        s = squeeze(obj.current_ADC_pre_processed(:,k)');
                        if 0~=sum(s);
                            mixedsig = [mixedsig; s];
                        end
                end
                icasig = fastica(mixedsig,'numofic',2);
                    s1 = icasig(1,:)';
                    s2 = icasig(2,:)';
                    
                    % debug
%                          icasig = fastica(mixedsig);
%                          h = figure(22);
%                          [S1,S2] = size(icasig);
%                          for k=1:S1
%                              figure(h);
%                              subplot(S1,1,k);
%                              plot(1:S2,icasig(k,:)','k.-');
%                          end
                    % debug

                        % rest is similar..
                        avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_ADC);
                        %
                        [s1,~] = TD_high_pass_filter( s1, avr_window );
                        [s2,~] = TD_high_pass_filter( s2, avr_window );
                        %                    
                        t1 = quantile(s1(:),obj.ADC_sgm_low_signal_quantile); %  take low 10% of signal
                        t2 = quantile(s1(:),obj.ADC_sgm_low_signal_quantile); %  
                        z1 = s1(s1<t1);
                        z2 = s2(s2<t2);
                        f1 = median(z1(:));
                        f2 = median(z2(:));
                        %
                        s_ = sqrt(abs(s1/f1.*s2/f2));
                        t = quantile(s_(:),obj.ADC_sgm_low_signal_quantile); 
                        z = s_(s_<t);
                        %
                        %t = obj.ADC_sgm_SN*median(z(:)); % threshold at 100X (or whatever) average noise level
                        t = 20*median(z(:)); % better
                        %
                        z = (s_ > t);
                        %                    
%figure(debug_h);subplot(2,4,k);plot(1:length(s_),s_,'k.-',1:length(s_),t*ones(1,length(s_)),'r-');grid on;xlabel(num2str(t));                    
                        %
                        min_size = round(obj.ADC_segm_Minimal_Trail_Duration*obj.Fs_ADC); % sic!
                        %
                        SE = strel('line',min_size,90);
                        z = imdilate(z,SE);
                        for k = 1 : num_ADC_channels% same for all
                            s = squeeze(obj.current_ADC_pre_processed(:,k)');
                            if 0~=sum(s);
                                obj.current_ADC_segmented(:,k) = z;
                            else
                                obj.current_ADC_segmented(:,k) = zeros(size(z));
                            end
                        end   
            end
            %                        
            if strcmp(type,'annotator"s b/g/s') % undo segmentations
                % same for all
                SGM = zeros(1,size(obj.current_ADC_segmented,1));
                                
                for k = 1:length(obj.current_annotation_time)
                    t0 = obj.current_annotation_time(k,1);
                    anno = obj.current_annotation(k,1);
                        lmin = fix((t0-3)*obj.Fs_ADC);
                        lmax = fix((t0+1)*obj.Fs_ADC);
                        if lmin>=1 && lmax<=length(SGM) && 0==sum(SGM(lmin:lmax)) % NB: if not empty on SGM
                            SGM(lmin:lmax)=1;
                        end
                end                
                %
                % the idea is to reject ROI if it doesn't contain enough
                % signal
%                 sgmntd_sgnls = zeros(size(SGM));
%                 for k = 1 : num_ADC_channels% same for all
%                     sgmntd_sgnls = sgmntd_sgnls | obj.current_ADC_segmented(:,k)';
%                 end
%                 %
%                 % does it work properly?
%                 Rt = 0.1;
%                 %
%                 SGM_lab = bwlabel(SGM);
%                 for l=1:max(SGM_lab)
%                     z = (SGM_lab==l);
%                     ratio = sum(z & sgmntd_sgnls)/sum(z);
%                     if ratio < Rt;
%                         SGM = SGM &~ z;
%                     end
%                 end                
                %
                % same for all
                for k = 1 : num_ADC_channels% same for all
                    if 0~=sum(obj.current_ADC_segmented(:,k));
                        obj.current_ADC_segmented(:,k) = SGM;
                    end
                end                                                                 
                % same for all
                %
                % set up available groups
                available_types_ind = unique(obj.current_annotation);
                obj.groups_available = cell(1,length(available_types_ind));
                for k=1:numel(available_types_ind)
                    obj.groups_available(k) = obj.groups_all(available_types_ind(k));
                end   
            end            
            if ~isempty(hw), delete(hw), drawnow; end;
            %            
            
                % one segmentation for all - for debugging
%                 SGM = obj.current_ADC_segmented(:,1);
%                 for k = 2 : num_ADC_channels
%                     SGM = SGM | obj.current_ADC_segmented(:,k);
%                 end
%                 SGM_lab = bwlabel(SGM);
%                 N_segmented = max(SGM_lab);
%                 N_annotated = length(obj.current_annotation);
%                 
%                 str = obj.current_filename;
%                 str = strrep(str,'SUBJECT_','');
%                 str = strrep(str,'.mat','');
%                 subj_num = str2num(str);                                
%                 disp([subj_num N_annotated N_segmented cnt]);                 
                
        end        
%-------------------------------------------------------------------------%        
        function [trails_features, feature_names] = extract_features_current_ADC(obj,~,~)

            feats = []; 
            fnames = [];    
            
            hw = waitbar(0,'Extracting trails features - please wait');
            %
            for subj_ind = 1:length(obj.subj_filenames)                                  
                %
                obj.switch_current_to_subject(char(obj.subj_filenames(subj_ind)));
                %
                num_ADC_channels = size(obj.current_data.ADC,2);        
                %                
                SGM = obj.get_joint_segmentation; % may be OR or AND regime
                %
                norm_meas = inf(1,num_ADC_channels);
                %
                num_ADC_channels = size(obj.current_data.ADC,2);
                
                data_LF_subtracted = zeros(size(obj.current_ADC_pre_processed));
                for k = 1 : num_ADC_channels
                    s = obj.current_ADC_pre_processed(:,k); 
                    % can't normalize otherwise - subtract LF trend
                    avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_ADC);                        
                    [s,~] = TD_high_pass_filter( s, avr_window );                    
                    data_LF_subtracted(:,k) = s;
                end                    
                %
                for k = 1 : num_ADC_channels
                    s = data_LF_subtracted(:,k); 
                    %    
                    if 0~=sum(s) && 0~= sum(squeeze(obj.current_ADC_segmented(:,k)))
                        fl = abs(s(~SGM));
                        norm_meas(k) = median(fl(:)); % norma == the median of segmentation-rejecetd part of record
                    end
                end
                %
                z_lab = bwlabel(SGM);
                %                                                                          
                t =(0:length(SGM)-1)'/obj.Fs_ADC; % seconds
                           
                %%%%%%%%%%%% MATCHING - preparing objects, start
                %
                ANNO_MAP = zeros(size(z_lab));
                % to assign annotations
                ANNO_NAT = zeros(1,max(z_lab)); % "natural"
                ANNO_DOC = zeros(1,max(z_lab)); % "doctor's"
                ANNO_AUT = zeros(1,max(z_lab)); % "automatic"
                %
                anno = obj.current_annotation;
                anno_t = obj.current_annotation_time;
                %
                delay = 0.5; % sec
                for k =1:length(anno)
                    A = anno(k);
                    pos = round((anno_t(k) - delay)*obj.Fs_ADC);
                    pos = max(1,pos);
                    pos = min(length(ANNO_MAP),pos);
                    ANNO_MAP(pos) = A;
                end;
                %
                % 8 sec "dead zone", where ROI is too far from annotations
                deadzone = ~imdilate(ANNO_MAP,strel('line',round(8*obj.Fs_ADC),90));
                %
                ANNO_MAP = dilate_labels(ANNO_MAP); % until 1 (maybe 2?)-gaps only
                %
                % DEBUG
                %figure(22);
                %plot(t,ANNO_MAP,'r-',t,SGM,'b-',t,-deadzone,'g-');                                    
                %
                %%%%%%%%%%%% MATCHING - preparing objects, ends
                
                feats_subj = [];
                fnames_subj = [];
                
                hw1 = waitbar(0,['quantifying ROIs - ' char(obj.subj_filenames(subj_ind))]);                
                for l=1:max(z_lab)
                    if ~isempty(hw1), waitbar(l/max(z_lab),hw1); drawnow, end;
                    %
                    % use "z_lab", "ANNO_MAP", and "deadzone" to do matching                    
                    ROI = (z_lab==l);
                    A = 0; % annotation
                    %
                    ROI_in_deadzone = (0~=sum(ROI&deadzone));
                    %
                    if ~ROI_in_deadzone
                        ROI_annos = ANNO_MAP(z_lab==l);
                        ROI_annos = ROI_annos(0~=ROI_annos);
                        ROI_annos_U = unique(ROI_annos);
                        if 1==length(ROI_annos_U)
                            A = ROI_annos_U;
                        elseif 2==length(ROI_annos_U) % judge by which one is more abundant
                            %
                            A1 = ROI_annos_U(1);
                            A2 = ROI_annos_U(2);
                            p1 = sum(A1==ROI_annos)/sum(ROI_annos(:));
                            if p1 > 0.6 % :) arbitrary
                                A = A1;
                            elseif p1 < 0.4
                                A = A2;
                            end
                            %
                        else
                            disp([l max(z_lab) ROI_annos_U']);
                        end
                    end
                    ANNO_NAT(1,l) = A;
                    % matching - finished
                    %
                    % for each segment.. look for best S/N channel
                    M = -Inf;
                    k_ = [];
                    for k = 1 : num_ADC_channels
                        ref = data_LF_subtracted(:,k);
                        signal = ref(z_lab==l);
                        measure = mean(abs(signal(:)))/norm_meas(k);
                        if measure >= M
                            k_=k;
                            M = measure;
                        end                        
                    end
                    %
                    ref = data_LF_subtracted(:,k_);
                    s_l = ref(z_lab==l);
                    %
                    % rest is more or less clear..
                        p1 = wentropy(s_l,'shannon')/length(s_l);
                        p2 = wentropy(s_l,'log energy')/length(s_l);
                            [C,L] = wavedec(s_l,6,'sym6');
                            [Ea,Ed] = wenergy(C,L); %these are normalized
                        p3 = Ea;
                        p4 = Ed; % this one, - contains 6 numbers
                        %
                        % start time, end time...
                        dt = t(z_lab==l);
                        t1 = min(dt(:));
                        t2 = max(dt(:));
                        %
                        % this is the place to calculate more features ...
                        %                    
                                                
                        [ Omega, psd ] = PSD( s_l,obj.Fs_ADC );
                         
                        I1 = sum(psd(1.1<=Omega&Omega<=1.5)); % Heartbeat
                        I2 = sum(psd(2.2<=Omega&Omega<=3)); % 2x Heartbeat
                        I3 = sum(psd(6<=Omega&Omega<=9)); % 8.3 Hz
                         
                        ILF = sum(psd(0<=Omega&Omega<=16));
                         
                        R1 = I1/ILF;
                        R2 = I2/ILF;
                        R3 = I3/ILF;
                        %                             
                        prm = [k_ l t1 t2 length(s_l) p1 p2 p3 R1 R2 R3 p4];
                        if 0~=sum(isinf(prm)) || 0~=sum(isnan(prm))
                            disp(prm);
                        else
                            feats_subj = [feats_subj; [subj_ind k_ l t1 t2, ...
                                ANNO_NAT(l), ANNO_DOC(l), ANNO_AUT(l),...
                                length(s_l)/obj.Fs_ADC p1 p2 p3 R1 R2 R3 p4]];
                            fnames_subj = [fnames_subj; cellstr(obj.current_filename)];
                        end                                                                                                    
                end % for l=1:max(z_lab)                                                            
                if ~isempty(hw1), delete(hw1), drawnow; end;
                %
                feats = [feats; feats_subj];
                fnames = [fnames; fnames_subj];
                if ~isempty(hw), waitbar(subj_ind/length(obj.subj_filenames),hw); drawnow, end;                                    
            end % for subj_ind = 1:length(obj.subj_filenames)
            if ~isempty(hw), delete(hw), drawnow; end;
            %
            feature_names = {'filename','subject_index','detector_index','trail_index','t1','t2', ...
                'anno_nat','anno_doc','anno_aut', ...
                'trail_length', ...
                'entropy','energy','Ea','R1','R2','R3','Ed1','Ed2','Ed3','Ed4','Ed5','Ed6'};
            obj.ADC_fv_all = {'trail_length', 'entropy','energy','Ea','R1','R2','R3','Ed1','Ed2','Ed3','Ed4','Ed5','Ed6'};
            %        
            trails_features = [fnames num2cell(feats)];            

        end
%-------------------------------------------------------------------------%                
        function pre_process_IMU(obj,type,~)
            
            if isempty(obj.current_data), return, end;
            %
            obj.current_IMU_pre_processed = zeros(size(obj.current_data.IMU));            
            num_IMU_channels = size(obj.current_data.IMU,2);
            for k = 1 : num_IMU_channels
                %
                s = obj.current_data.IMU(:,k);
                    avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_IMU);
                    %
                    [s,~] = TD_high_pass_filter( s, avr_window );
                    %                    
                    obj.current_IMU_pre_processed(:,k) = sqrt(s.*s);                
            end
        end
%-------------------------------------------------------------------------%                
        function segment_IMU(obj,type,~)        
            if isempty(obj.current_data), return, end;

            obj.current_IMU_segmented = zeros(size(obj.current_data.IMU));            
            %            
%debug_h = figure;            
            num_IMU_channels = size(obj.current_data.IMU,2);
            hw = waitbar(0,[type ' segmenting IMU - please wait']);
            for k = 1 : num_IMU_channels
                if ~isempty(hw), waitbar(k/num_IMU_channels,hw); drawnow, end;
                s_ = obj.current_IMU_pre_processed(:,k);
                    %
%                     [N,X] = hist(s_(:),300);
%                     T = find(N==max(N(:)));                    
%                     t = min(s_(:)) + 100*( X(min(T(:))) - min(s_(:)) );
%                     %                                        
%                     z = (s_ > t);        
                    %
                    t = quantile(s_(:),obj.IMU_sgm_low_signal_quantile); %  take low 10% of signal
                    z = s_(s_<t);
                    %
                    t = obj.IMU_sgm_SN*median(z(:)); % threshold at 200X (or whatever) average noise level
                    %
                    if t < obj.IMU_sgm_signal_cutoff, t = Inf; end; % PRECAUTION AGAINST TOO NOISY SIGNALS                    
                    %
                    z = (s_ > t);
                    %                    
%figure(debug_h);subplot(1,num_IMU_channels,k);semilogy(1:length(s_),s_,'k.-',1:length(s_),t*ones(1,length(s_)),'r-');grid on;xlabel(num2str(t));
                    %
                    min_size = round(obj.ADC_segm_Minimal_Trail_Duration*obj.Fs_IMU); % sic!
                    %
                    SE = strel('line',min_size,90);
                    z = imdilate(z,SE); % to expand exclusion area a bit.. 
                    obj.current_IMU_segmented(:,k) = z;
                    %
            end
            if ~isempty(hw), delete(hw), drawnow; end;                                                
        end
%-------------------------------------------------------------------------%                
        function calculate_PSDs(obj,~,~)        
            if isempty(obj.current_data), return, end;

            obj.current_ADC_PSD = [];                
            obj.current_ADC_preprocessed_PSD = [];         
                        
            num_ADC_channels = size(obj.current_data.ADC,2);
            hw = waitbar(0,['Calculating ADC PSDs - please wait']);
            for k = 1 : num_ADC_channels
                if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                s = obj.current_data.ADC(:,k);
                s_prp = obj.current_ADC_pre_processed(:,k);
                [obj.current_omega_ADC, psd] = PSD(s,obj.Fs_ADC);
                [obj.current_omega_ADC, psd_prp] = PSD(s_prp,obj.Fs_ADC);
                if isempty(obj.current_ADC_PSD)
                    L = length(psd);
                    obj.current_ADC_PSD = zeros(L,num_ADC_channels);   
                    obj.current_ADC_preprocessed_PSD = zeros(L,num_ADC_channels);
                end;
                obj.current_ADC_PSD(:,k) = psd;
                obj.current_ADC_preprocessed_PSD(:,k) = psd_prp;                                
            end
            if ~isempty(hw), delete(hw), drawnow; end;
            %
            obj.current_IMU_PSD = [];                
            obj.current_IMU_preprocessed_PSD = [];         
                        
            num_IMU_channels = size(obj.current_data.IMU,2);
            hw = waitbar(0,'Calculating IMU PSDs - please wait');
            for k = 1 : num_IMU_channels
                if ~isempty(hw), waitbar(k/num_IMU_channels,hw); drawnow, end;
                s = obj.current_data.IMU(:,k);
                s_prp = obj.current_IMU_pre_processed(:,k);
                [obj.current_omega_IMU, psd] = PSD(s,obj.Fs_IMU);
                [obj.current_omega_IMU, psd_prp] = PSD(s_prp,obj.Fs_IMU);
                if isempty(obj.current_IMU_PSD)
                    L = length(psd);
                    obj.current_IMU_PSD = zeros(L,num_IMU_channels);   
                    obj.current_IMU_preprocessed_PSD = zeros(L,num_IMU_channels);
                end;
                obj.current_IMU_PSD(:,k) = psd;
                obj.current_IMU_preprocessed_PSD(:,k) = psd_prp;                                
            end
            if ~isempty(hw), delete(hw), drawnow; end;
                        
        end        
%-------------------------------------------------------------------------%  
        function [coords,IDX] = perform_unsupervised_clustering(obj,type,n_clusters,~)
            
            coords = [];
            IDX = [];
            
            if isempty(obj.ADC_trails_features_data), return, end;

            [~,record_length] = size(obj.ADC_trails_features_data);
            
            [data,annos] = obj.get_selected_feature_vector_data;
            
            switch type

                case 'KMeans'
                    try                        
                        IDX = kmeans(data,n_clusters);
                        %[~,~,stats] = manova1(data,IDX); ehm, fails
                        %
                        [coeff,score] = pca(data);
                        coords = score; % score*coeff'; % ??? that gives something weird
                    catch
                        errordlg('kmeans doesnt like these data, no output provided')
                    end
                    
            end %switch type                        
                        
        end
%-------------------------------------------------------------------------%  
        function [coords,IDX] = get_canons_for_supervised_classification(obj,type,~)
            
            coords = [];
            IDX = [];
            
            if isempty(obj.ADC_trails_features_data), return, end;
                                    
            % display annotated data in canonic coords
            switch type                
                case {'annotator"s + segmentation','annotator"s b/g/s'}
                % data composed with selected featue vector but NOT filtered re selected groups    
                [data,IDX] = obj.get_annotators_categorized_data('selected components');                
                case 'auto annotated'
                    [data,IDX] = obj.get_auto_categorized_data('selected components');                    
                otherwise                    
            end
            %
            % fix the data by retaining only wanted groups - starts
            data_reducted = [];
            IDX_reducted = [];
            for k = 1:numel(IDX)
                if ~isempty(intersect(obj.groups_all(IDX(k)),obj.groups_selected))
                    data_reducted = [data_reducted; data(k,:)];
                    IDX_reducted = [IDX_reducted IDX(k)];
                end
            end
            data = data_reducted;
            IDX = IDX_reducted;
            % fix the data by retaining only wanted groups - end
            %
            if ~isempty(IDX)
                [~,score] = pca(data);
                data = score;
                [~,~,stats] = manova1(data,IDX);            
                coords = stats.canon;
            end
                                    
        end        
%-------------------------------------------------------------------------% 
        function [data, annos] = get_selected_feature_vector_data(obj,~,~)
            data = [];
            [~,record_length] = size(obj.ADC_trails_features_data);
            annos = cell2mat(obj.ADC_trails_features_data(:,7)); % annotations column            
            all_data = cell2mat(obj.ADC_trails_features_data(:,10:record_length)); % 10 is offset
            assert(size(all_data,2)==numel(obj.ADC_fv_all));
            for m = 1:numel(obj.ADC_fv_all),
                for k = 1:numel(obj.ADC_fv_selected),
                    if strcmp(char(obj.ADC_fv_all(m)),char(obj.ADC_fv_selected(k)))
                        data = [data all_data(:,m)];
                    end
                end
            end            
        end
%-------------------------------------------------------------------------% 
        function [X1,Y1,X2,Y2,z] = get_pairwise_comparison(obj,feature_vector_name,group1,group2,type,~)
            X1 = [];
            Y1 = [];
            X2 = [];
            Y2 = [];
            
            z = [];
            
            if strcmp(group1,group2), return, end;
            
            if isempty(obj.ADC_trails_features_data), return, end;
                                    
            % display annotated data in canonic coords
            switch type                
                case {'annotator"s + segmentation','annotator"s b/g/s'}
                    % data composed with selected featue vector but NOT filtered re selected groups    
                    [data,IDX] = obj.get_annotators_categorized_data('all components');                
                case 'auto annotated'
                    [data,IDX] = obj.get_auto_categorized_data('all components');                    
                otherwise                    
            end

            %compile 2 corresponding satistical samples and evaluate them
            x1 = [];
            x2 = [];
            %
            fv_index = find(0~=strcmp(obj.ADC_fv_all,feature_vector_name));
            %
            for k = 1:numel(IDX)
                if strcmp(obj.groups_all(IDX(k)),group1)
                    x1 = [x1; data(k,fv_index)];
                elseif strcmp(obj.groups_all(IDX(k)),group2)
                    x2 = [x2; data(k,fv_index)];
                end
            end
            
            %%%%%%%%%%%
            % Cohen's d
            N1 = numel(x1);
            N2 = numel(x2);
            m1 = mean(x1);
            m2 = mean(x2);
            std1 = std(x1);
            std2 = std(x2);
            
            s = sqrt( 1/(N1+N2)*( (N1-1)*var(x1) + (N2-1)*var(x2) ) );
            z.d = abs( m1 - m2 )/s;

            % AUC
            [SEN, SPEC, TH, ACC, AUC,Yi,idx]=roc1(x1,x2);
            %plot(1-SPEC,SEN);

            if AUC < 0.5 AUC = 1-AUC; end;		

            z.AUC = AUC;
%%%%%%%%%%%% should be controlled
nbins = 50;            
a_ks = 0.01;
a_ranksum = 0.01;
%%%%%%%%%%%%
            
            %KS
            [h_ks,p_ks,ks2stat] = kstest2(x1,x2,a_ks);
            [p_rnk,h_rnk] = ranksum(x1,x2, a_ranksum);

            z.h_ks = h_ks;
            z.p_ks = p_ks;
            z.h_rnk = h_rnk;
            z.p_rnk = p_rnk;            
            z.N1=N1;
            z.N2=N2;
            z.m1=m1;
            z.m2=m2;
            z.std1=std1;
            z.std2=std2;
            %

            %            
            [ min1, max1, X1, Y1 ] = histodata(x1,nbins);
            [ min2, max2, X2, Y2 ] = histodata(x2,nbins);
            %
            range = max([ max1 max2 ]) - min( [min1 min2 ]);
            %
            c1 = (max1-min1)/range;
            c2 = (max2-min2)/range;
            Y1=Y1/c1;
            Y2=Y2/c2;
                           
        end
%-------------------------------------------------------------------------%         
        function [data,IDX] = get_annotators_categorized_data(obj,mode,~)
            % data composed with selected featue vector but NOT filtered re selected groups    
            switch mode
                case 'selected components'
                    [fv_data,annos] = obj.get_selected_feature_vector_data;            
                case 'all components'
                    fv_temp = obj.ADC_fv_selected;
                    obj.ADC_fv_selected = obj.ADC_fv_all;
                    [fv_data,annos] = obj.get_selected_feature_vector_data;
                    obj.ADC_fv_selected = fv_temp;
            end
            
            t1 = cell2mat(obj.ADC_trails_features_data(:,5));
            %
            data = [];
            IDX = [];  
                                            
                   for subj = 1:numel(obj.subj_data)
                        %
                        subj_name = obj.subj_data(subj).filename;
                        %
                        hw = waitbar(0,[ subj_name ' annotations - please wait']);
       
%                           % RANDOMIZE!
%                           ANNOS = unique(anno);
%                              for k = 1:length(t1)
%                                  subj_index_k = cell2mat(obj.ADC_trails_features_data(k,2)); % faster via index
%                                  if (subj_index_k == subj)
%                                      L = length(ANNOS);
%                                      IDX = [IDX ANNOS(1+round((L-1)*rand))];
%                                      data = [data; fv_data(k,:)];
%                                  end
%                              end
%                           % RANDOMIZE!                             

                             for k = 1:length(t1)
                                 subj_index_k = cell2mat(obj.ADC_trails_features_data(k,2)); % faster via index
                                 anno = annos(k);
                                 if (subj_index_k == subj) && 0~=anno
                                     IDX = [IDX anno];
                                     data = [data; fv_data(k,:)];                                     
                                 end                                                                                                   
                             end                          
                          
                        if ~isempty(hw), delete(hw), drawnow; end;           
                                                
                   end  
                                     
                    % set up available groups
                    available_types_ind = unique(IDX);
                    obj.groups_available = cell(1,length(available_types_ind));
                    for k=1:numel(available_types_ind)
                        obj.groups_available(k) = obj.groups_all(available_types_ind(k));
                    end   
                    
                    % fix selected groups if needed
                    if isempty(intersect(obj.groups_available,obj.groups_selected))
                        obj.groups_selected = obj.groups_available;
                    end  
        end
%-------------------------------------------------------------------------%
        function data = get_annotators_uncategorized_data(obj,mode,~)
            % data composed with selected featue vector but NOT filtered re selected groups    
            switch mode
                case 'selected components'
                    [fv_data,annos] = obj.get_selected_feature_vector_data;            
                case 'all components'
                    fv_temp = obj.ADC_fv_selected;
                    obj.ADC_fv_selected = obj.ADC_fv_all;
                    [fv_data,annos] = obj.get_selected_feature_vector_data;
                    obj.ADC_fv_selected = fv_temp;
            end
            
            t1 = cell2mat(obj.ADC_trails_features_data(:,5));
            %
            data = [];                                            
                   for subj = 1:numel(obj.subj_data)                        
                        subj_name = obj.subj_data(subj).filename;
                        %
                        hw = waitbar(0,[ subj_name ' picking up the un-annotated data - please wait']);     
                             for k = 1:length(t1)
                                 subj_index_k = cell2mat(obj.ADC_trails_features_data(k,2)); % faster via index
                                 anno = annos(k);
                                 if (subj_index_k == subj) && 0==anno
                                     data = [data; fv_data(k,:)];                                     
                                 end                                                                                                   
                             end                                                    
                        if ~isempty(hw), delete(hw), drawnow; end;                                                           
                   end                                       
        end
%-------------------------------------------------------------------------%
        function possibly_exclude_ADC_findings_with_simultaneous_IMU_response(obj,~,~)
        % NB - this function presumes that both ADC and IMU
        % segmentations are stored in CURRENT, and excludes ADC findings
        % if intersected with IMU
            if obj.exclude_IMU
                    IMU_mask = sum(obj.current_IMU_segmented,2);
                    IMU_mask(0~=IMU_mask)=1;   
                    N2 = size(obj.current_ADC_segmented,1);
                    exclusion_ADC_mask = (imresize(double(IMU_mask),[N2,1])>0.5);
                    %
                    num_ADC_channels = size(obj.current_ADC_segmented,2);
                    hw = waitbar(0,'excluding IMU findings - please wait');
                    for k = 1 : num_ADC_channels  
                        if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                        s = obj.current_ADC_segmented(:,k);
                        %
                        labsegm = bwlabel(s);
                        exclabs = unique(labsegm.*exclusion_ADC_mask);
                        for elb=1:numel(exclabs)
                            L = exclabs(elb);
                            if 0~=L
                               labsegm(labsegm==L)=0;
                            end
                        end                        
                        obj.current_ADC_segmented(:,k) = (labsegm~=0);
                    end
                    if ~isempty(hw), delete(hw), drawnow; end;                                                            
            end;
        end        
%-------------------------------------------------------------------------%
        function [groups,CM] = calculate_confusion_matrix(obj,source_type,method,~)
            
            groups = [];
            CM = [];
            
                switch source_type                
                    case {'annotator"s + segmentation','annotator"s b/g/s'}
                        % data composed with selected featue vector but NOT filtered re selected groups    
                        [data,IDX1] = obj.get_annotators_categorized_data('selected components');                                                                                                
                    case 'auto annotated'                        
                        [data,IDX1] = obj.get_auto_categorized_data('selected components');                        
                    otherwise                    
                end
                
                % data = rand(size(data)); % randomize :)
                
                % fix the data by retaining only wanted groups - starts
                data_reducted = [];
                IDX1_reducted = [];
                for k = 1:numel(IDX1)
                    if ~isempty(intersect(obj.groups_all(IDX1(k)),obj.groups_selected))
                        data_reducted = [data_reducted; data(k,:)];
                        IDX1_reducted = [IDX1_reducted; IDX1(k)];
                    end
                end
                data = data_reducted;
                IDX1 = IDX1_reducted;
                % fix the data by retaining only wanted groups - ends
                
                % PCA data if needed, & choose first 3 components
                truncation = 3;                                
                if ~isempty(strfind(method,'PCA')) && size(data,2)>=truncation
                    V = cov(data);
                    %V = cov(zscore(data)); % ?? no diff..
                    SD = sqrt(diag(V));
                    R = V./(SD*SD');
                    COEFF = pcacov(R);
                    U = data*COEFF;
                    data = U(:,1:truncation);                                
                end
                % PCA data if needed, & choose first 3 compnents
                
                %
                if isempty(data), return, end;
                %
                N = numel(IDX1);
                IDX2 = zeros(size(IDX1));
                %
                groups = sort(unique(IDX1)); 
                %
                        hw = waitbar(0,'"take one out"-ing - please wait');
                        try
                            for i = 1:N    
                                if ~isempty(hw), waitbar(i/N,hw); drawnow, end;
                                %
                                include = setdiff(1:N,i);
                                train_data = data(include,:);      
                                train_group = IDX1(include,:);
                                %
                                switch method                
                                    case {'Linear','Quadratic','PCA->Linear','PCA->Quadratic'}
                                        if ~isempty(strfind(method,'Linear'))
                                            IDX2(i) = classify(data(i,:),train_data,train_group,'Linear');
                                        else
                                            IDX2(i) = classify(data(i,:),train_data,train_group,'Quadratic');
                                        end
                                    case {'kNN','PCA->kNN'}
                                        mdl = fitcknn(train_data,train_group);
                                        IDX2(i) = predict(mdl,data(i,:));
                                end
                                %                            
                            end
                        catch
                            if ~isempty(hw), delete(hw), drawnow; end;                                                        
                                errordlg('classify error - possibly not enough data..');
                                    CM = zeros(numel(groups));
                                        return;
                        end
                        if ~isempty(hw), delete(hw), drawnow; end;
                %                           
                CM = zeros(numel(groups));
                %
                %
                for k=1:N,
                    g1 = IDX1(k);
                    g2 = IDX2(k);
                    i1 = find(groups==g1);
                    i2 = find(groups==g2);
                    CM(i1,i2) = CM(i1,i2)+1;
                end
                %CM = CM/sum(CM(:));                  
        end                                                    
%-------------------------------------------------------------------------%                 
        function [groups,group_nums] = classify_unannotated_data(obj,method,~)
            
            groups = [];
            group_nums = [];
            
            [data,IDX1] = obj.get_annotators_categorized_data('selected components');                                                                                                
                                
            % fix the data by retaining only wanted groups - starts
            data_reducted = [];
            IDX1_reducted = [];
            for k = 1:numel(IDX1)
                if ~isempty(intersect(obj.groups_all(IDX1(k)),obj.groups_selected))
                    data_reducted = [data_reducted; data(k,:)];
                    IDX1_reducted = [IDX1_reducted; IDX1(k)];
                end
            end
            
            TR_data = data_reducted;
            TR_IDX = IDX1_reducted;
            % fix the data by retaining only wanted groups - ends
                
            TST_data = obj.get_annotators_uncategorized_data('selected components');
            
            % PCA data if needed, & choose first 3 components
            truncation = 3;                                
            if ~isempty(strfind(method,'PCA')) && size(TR_data,2)>=truncation
                    V = cov(TR_data);
                    %V = cov(zscore(data)); % ?? no diff..
                    SD = sqrt(diag(V));
                    R = V./(SD*SD');
                    COEFF = pcacov(R);
                    %
                    U = TR_data*COEFF;
                    TR_data = U(:,1:truncation);                                
                    %
                    U = TST_data*COEFF;
                    TST_data = U(:,1:truncation);                                                                        
            end
            % PCA data if needed, & choose first 3 compnents                                             
                    %
                    try                                
                        switch method                
                            case {'Linear','Quadratic','PCA->Linear','PCA->Quadratic'}
                                if ~isempty(strfind(method,'Linear'))
                                    IDX2 = classify(TST_data,TR_data,TR_IDX,'Linear');
                                else
                                    IDX2 = classify(TST_data,TR_data,TR_IDX,'Quadratic');
                                end
                            case {'kNN','PCA->kNN'}
                                mdl = fitcknn(TR_data,TR_IDX);
                                IDX2 = predict(mdl,TST_data);
                        end
                    catch
                        errordlg('classify error - possibly not enough data..');
                        return;
                    end
                    
            groups = sort(unique(TR_IDX));            
            group_nums = zeros(1,length(groups));
            for k = 1:length(group_nums)
                group_nums(k) = sum(groups(k)==IDX2(:));
            end
        end                                                                                    
%-------------------------------------------------------------------------%         
        % "mocking" annotation function for debug
        function [data,IDX] = get_auto_categorized_data(obj,mode,~)
            % data composed with selected featue vector but NOT filtered re selected groups                            
            data = [];
            IDX = [];  

            [Nd,Nrec] = size(obj.ADC_trails_features_data);
            
            hw = waitbar(0,'auto annotating - please wait');
            for k=1:Nd,
                if ~isempty(hw), waitbar(k/Nd,hw); drawnow, end;
                fv = cell2mat(obj.ADC_trails_features_data(k,10:Nrec));
                %
                % 1) multi-threshold and define IDX 
                R3 = fv(7);
                Le = fv(1);
                Ea = fv(4);
                % En = fv(2);
                %
                % groups_all = {'breathe','general','head','limb','startle','other'};
                high_energy = true;
                extra_short = true;
                harmonic = true;
                if R3<0.03 % ?
                    harmonic = false;
                end
                if Ea<95 
                    high_energy = false;
                end
                if Le>513
                    extra_short = false;
                end
                %
                type = 10;       % T4
                if extra_short && high_energy && ~harmonic
                    type = 7;    % T1
                elseif ~extra_short && ~harmonic
                    type = 8;    % T2
                elseif harmonic && high_energy
                    type = 9;    % T3
                end                               

                %                
                % 2) choose fv components                
                switch mode
                    
                    case 'selected components'
                        new_fv = [];
                        for m = 1:numel(obj.ADC_fv_all),            
                            for p = 1:numel(obj.ADC_fv_selected),
                                if strcmp(char(obj.ADC_fv_all(m)),char(obj.ADC_fv_selected(p)))
                                    new_fv = [new_fv fv(:,m)];
                                end                
                            end
                        end
                        data = [data; new_fv];
                                                                                                
                    case 'all components'
                        data = [data; fv];
                end
                %
                IDX = [IDX; type];
            end
                        
            if ~isempty(hw), delete(hw), drawnow; end;

                    % set up available groups
                    available_types_ind = unique(IDX);
                    obj.groups_available = cell(1,length(available_types_ind));
                    for k=1:numel(available_types_ind)
                        obj.groups_available(k) = obj.groups_all(available_types_ind(k));
                    end
                    %
                    % fix selected groups if needed
                    if isempty(intersect(obj.groups_available,obj.groups_selected))
                        obj.groups_selected = obj.groups_available;
                    end            
        end
%-------------------------------------------------------------------------%
function SGM = get_joint_segmentation(obj,~)
                num_ADC_channels = size(obj.current_data.ADC,2); 
                if strcmp(obj.SGM_MODE,'OR') % one segmentation for all
                    SGM = obj.current_ADC_segmented(:,1);
                    for k = 2 : num_ADC_channels
                        SGM = SGM | obj.current_ADC_segmented(:,k);
                    end
                elseif strcmp(obj.SGM_MODE,'AND2')
                    z = obj.current_ADC_segmented(:,1);
                    for k = 2 : num_ADC_channels
                        if 0~=sum(z)                    
                            z = z + obj.current_ADC_segmented(:,k);
                        end
                    end
                    z = (z>1);
                    % dilate
                    min_size = obj.ADC_segm_Minimal_Trail_Duration*obj.Fs_ADC;
                    SE = strel('line',min_size,90);
                    z = imdilate(z,SE);                            
                    % dialte
                    SGM = z;
                elseif strcmp(obj.SGM_MODE,'AND2_best_SN') || strcmp(obj.SGM_MODE,'AND2_larger_total')
                    %
                    norm_meas = inf(1,num_ADC_channels);  % fluctuations
                    sig_meas = zeros(1,num_ADC_channels); % signal
                    sgm_meas = zeros(1,num_ADC_channels); % total area of segmented stuff
                    %
                    num_ADC_channels = size(obj.current_data.ADC,2);

                    data_LF_subtracted = zeros(size(obj.current_ADC_pre_processed));
                    for k = 1 : num_ADC_channels
                        s = obj.current_ADC_pre_processed(:,k); 
                        % can't normalize otherwise - subtract LF trend
                        avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_ADC);                        
                        [s,~] = TD_high_pass_filter( s, avr_window );                    
                        data_LF_subtracted(:,k) = s;
                    end                    
                    %
                    IMU_mask = sum(obj.current_IMU_segmented,2);
                    IMU_mask(0~=IMU_mask)=1;   
                    N2 = size(obj.current_ADC_segmented,1);
                    exclusion_ADC_mask = (imresize(double(IMU_mask),[N2,1])>0.5);                
                    %
                    for k = 1 : num_ADC_channels
                        s = data_LF_subtracted(:,k); 
                        %    
                        SGM_K = squeeze(obj.current_ADC_segmented(:,k));
                        if 0~=sum(s) && 0~= sum(SGM_K)
                            z = s(s~=SGM_K);
                            z = z(z~=exclusion_ADC_mask);
                            fl = abs(z);
                            sg = abs(s(SGM_K>0));
                            norm_meas(k) = median(fl(:));   % norma == the median of segmentation-rejecetd part of record
                            sig_meas(k) = mean(sg(:));      % the median of segmented part of record
                            sgm_meas(k) = sum(SGM_K);
                        end
                    end
                    %                
                    % two modes of qualifier - best S/N or larger total
                    % segmentation
                    if strcmp(obj.SGM_MODE,'AND2_best_SN')
                        qlfr = sig_meas./norm_meas;
                    elseif strcmp(obj.SGM_MODE,'AND2_larger_total')                        
                        qlfr = sgm_meas;
                    end
                    %                
                    % find two "best" detectors, basing on "qualifier"
                    qlfr_srtd = sort(qlfr);
                    % find two "best" detectors
                    D1_meas = qlfr_srtd(length(qlfr_srtd));
                    D2_meas = qlfr_srtd(length(qlfr_srtd)-1);
                    D1_ind = 0;
                    D2_ind = 0;
                    for k=1:length(qlfr)
                        if qlfr(k)==D1_meas
                            D1_ind=k;
                        end
                        if qlfr(k)==D2_meas
                            D2_ind=k;
                        end                    
                    end
                    %
                    disp('-----------------------');
                    disp(char(obj.current_filename));
                    disp([D1_ind D2_ind]);
                    for k = 1 : num_ADC_channels
                        disp([k norm_meas(k) sig_meas(k) qlfr(k)])
                    end
                    disp('-----------------------');
                    %
                    SGM = obj.current_ADC_segmented(:,D1_ind) & obj.current_ADC_segmented(:,D2_ind);                                        
                end                
end
                
%-------------------------------------------------------------------------%         
    end % methods            
end
