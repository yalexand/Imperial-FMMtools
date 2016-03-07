
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
                        
        corr_map_W = 80;
        
    end                    
    
    properties(Transient)
        
        DefaultDirectory = ['C:' filesep];
        RootDirectory = ['C:' filesep];       
                
        current_filename = [];;
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fs_ADC = 1024; %Hz
        Fs_IMU = 64; %Hz

        preprocessing_types = {'freq.notch-comb + median','none'};
        ADC_prprss_Median_Size = 3;
        ADC_prprss_notch_comb_fb = 39.1328; % base frequency, Hz
        ADC_prprss_notch_comb_bw = 0.0022; % bandwidth        
        %
        segmentation_types = {'moving average subtraction + thresholding','...'};
        ADC_segm_Moving_Average_Window = 5; % seconds!
        ADC_segm_Threshold = 0.2; % above noise trail
        ADC_segm_Minimal_Trail_Duration = 3;  %seconds
                        
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
            
         end
%-------------------------------------------------------------------------%        
        function load_single_subject(obj,verbose,~)
                        
            [filename,pathname] = uigetfile({'*.mat','Subject Records Files'},'Select data file',obj.DefaultDirectory);
            if filename == 0, return, end;       

            obj.clear_current_data;                                    
            obj.ADC_trails_features_data = [];
            obj.open_datafile([pathname filesep filename]);
                 
            obj.current_filename = filename;                    
            obj.DefaultDirectory = pathname;
            
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
                settings.ADC_segm_Threshold = obj.ADC_segm_Threshold;
                settings.ADC_segm_Minimal_Trail_Duration = obj.ADC_segm_Minimal_Trail_Duration;                                
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
                obj.ADC_segm_Threshold = settings.ADC_segm_Threshold;
                obj.ADC_segm_Minimal_Trail_Duration = settings.ADC_segm_Minimal_Trail_Duration;                                                
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
                        
        end
%-------------------------------------------------------------------------%        
        function segment_ADC(obj,type,~)
            
            if isempty(obj.current_ADC_pre_processed), return, end;

            obj.current_ADC_segmented = zeros(size(obj.current_data.ADC));
                                             
            num_ADC_channels = size(obj.current_data.ADC,2);
            hw = waitbar(0,[type ' segmenting ADC - please wait']);
            for k = 1 : num_ADC_channels
                if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                s = obj.current_ADC_pre_processed(:,k);
                if strcmp(type,'moving average subtraction + thresholding')
                    %
                    avr_window = round(obj.ADC_segm_Moving_Average_Window*obj.Fs_ADC);
                    %
                    [s,~] = TD_high_pass_filter( s, avr_window );
                    %                    
                    s2 = s.*s;
                    t = quantile(s2(:),obj.ADC_segm_Threshold);
                    z = s2(s2>t);
                    t = median(z(:)) + 3*std(z(:));
                    %                                        
                    z = (s2 > t);
                    min_size = round(obj.ADC_segm_Minimal_Trail_Duration*obj.Fs_ADC);
                    SE = strel('line',min_size,90);
                    z = imclose(z,SE);
                    z = bwareaopen(z,min_size);
                    obj.current_ADC_segmented(:,k) = z;
                    %
                else % shouldn't happen
                    obj.current_ADC_segmented(:,k) = zeros(size(s)); % stupid
                end
            end
            
            % extract annotation

            event_times = obj.current_data.US(:,1);
            event_types = obj.current_data.US(:,3); % 1,2,3,4,5 -> b,g,h,l,s            
            obj.current_annotation = zeros(size(event_times,1),1);
            obj.current_annotation_time = zeros(size(event_times,1),1);            
            if ~isempty(hw), delete(hw), drawnow; end;
            %
            type_ind = 0; % for now...
            %
            for k = 1:length(event_times)
                obj.current_annotation_time(k,1) = cell2mat(obj.current_data.US(k,1));
                switch char(event_types(k,1))
                    case 'b'                    
                        type_ind = 1;                        
                    case 'g'
                        type_ind = 2;
                    case 'h'
                        type_ind = 3;
                    case 'l'
                        type_ind = 4;
                    case 's'    
                        type_ind = 5;                        
                end
                obj.current_annotation(k,1) = type_ind;
            end
            % extract annotation

        end
        
%-------------------------------------------------------------------------%        
        function [trails_features, feature_names] = extract_features_current_ADC(obj,~,~)

            feats = [];
            
            hw = waitbar(0,'Extracting trails features - please wait');
            num_ADC_channels = size(obj.current_data.ADC,2);
            for k = 1 : num_ADC_channels
                if ~isempty(hw), waitbar(k/num_ADC_channels,hw); drawnow, end;
                feats_k = [];
                s = obj.current_ADC_pre_processed(:,k); 
                z = obj.current_ADC_segmented(:,k);
                % first, one needs to normalize the signal, e.g. by dividing it by the
                % intensity of its fluctuation level
                fl = abs(s(~z));
                fl_t = quantile(fl(:),0.75);
                fl = fl(fl<fl_t); % get the weakest 75% part
                s = s/std(fl(:)); % normalize by its std
                %
                z_lab = bwlabel(z);
                for l=1:max(z_lab)
                    s_l = s(z_lab==l);
                    p1 = wentropy(s_l,'shannon')/length(s_l);
                    p2 = wentropy(s_l,'log energy')/length(s_l);
                        [C,L] = wavedec(s_l,4,'sym4');
                        [Ea,Ed] = wenergy(C,L); %these are normalized
                    p3 = Ea;
                    p4 = Ed; % this one, - contains 4 numbers                    
                    %
                    % this is the place to calculate more features ...
                    %                    
                    feats_k = [feats_k; [k l length(s_l) p1 p2 p3 p4]];
                end
                feats = [feats; feats_k];
            end
            if ~isempty(hw), delete(hw), drawnow; end;
            %
            % add new features names here
            feature_names = {'filename','detector_index','trail_index','trail_length', ...
                'entropy','energy','Ea','Ed1','Ed2','Ed3','Ed4'};
            %
            Nrec = size(feats,1);
            trails_features = [cellstr(repmat(obj.current_filename,Nrec,1)) num2cell(feats)];            
                                                
        end
%-------------------------------------------------------------------------%                
        function pre_process_IMU(obj,type,~)
            
            if isempty(obj.current_data), return, end;

            obj.current_IMU_pre_processed = zeros(size(obj.current_data.IMU));            
            %
            % to do
            %
        end
%-------------------------------------------------------------------------%                
        function segment_IMU(obj,type,~)        
            if isempty(obj.current_data), return, end;

            obj.current_IMU_segmented = zeros(size(obj.current_data.IMU));            
            %
            % to do
            %            
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
            
            data = cell2mat(obj.ADC_trails_features_data(:,5:record_length));
            
            switch type
                
                case 'KMeans'
                  IDX = kmeans(data,n_clusters);
                  %[~,~,stats] = manova1(data,IDX); ehm, fails
                  %
                  [coeff,score] = pca(data);
                  coords = score; % score*coeff'; % ??? that gives something weird
                  
            end                        
                        
        end
    end % methods
%-------------------------------------------------------------------------%      
            
end