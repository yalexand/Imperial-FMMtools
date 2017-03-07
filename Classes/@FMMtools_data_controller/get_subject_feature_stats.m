function stats = get_subject_feature_stats(obj,~,~)

stats = [];

if isempty(obj.ADC_trails_features_data), return, end;
            
%%%%%%%%% pat rec
            [ALL_fv_data,ALL_annos] = obj.get_selected_feature_vector_data;

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
%%%%%%%%% pat rec

            names = {'subject','tot time','tot_ROI_time','num_ROIs', 'num_annotations',...
                'N_anno_B','N_anno_G','N_anno_H','N_anno_L','N_anno_S','N_anno_O', ...
                'N_proj_B','N_proj_G','N_proj_H','N_proj_L','N_proj_S','N_proj_O'};
            for k = 1:length(obj.ADC_fv_all)
                names = [names ['mean_' char(obj.ADC_fv_all(k))] ];
                names = [names ['std_' char(obj.ADC_fv_all(k))] ];
                names = [names ['Q25_' char(obj.ADC_fv_all(k))] ];
                names = [names ['median_' char(obj.ADC_fv_all(k))] ];
                names = [names ['Q75_' char(obj.ADC_fv_all(k))] ];                
            end
            
            d2 = obj.annotators_delay;
            d1 = obj.annotators_reaction; % [second] - minimal discernable time between event and annotation
            
            logdata = [];
            
            hw = waitbar(0,'calculating subject stats - please wait');
            for subj_ind = 1:length(obj.subj_filenames)
                if ~isempty(hw), waitbar(subj_ind/length(obj.subj_filenames),hw); drawnow, end;
                obj.switch_current_to_subject(char(obj.subj_filenames(subj_ind)));
                %
                SGM = obj.get_joint_segmentation; % may be OR or AND regime
                %
                anno = obj.subj_data(subj_ind).annotation;
                anno_t = obj.subj_data(subj_ind).annotation_time; 
                N_annotated = length(anno);
                %
                cnt = zeros(1,6);
                cnt_projected = zeros(1,6);
                %
                anno_tot = zeros(size(SGM));
                %
                for k =1:length(anno)
                    A = anno(k);
                    cnt(A) = cnt(A)+1;
                    %
                    T1 = anno_t(k)-3.5; % seconds
                    T2 = anno_t(k)+1.5; % seconds
                    %
                    L = round((T1+T2)/2*obj.Fs_ADC);
                    DL = round(d2/2*obj.Fs_ADC);
                    L1 = L-DL;
                    L2 = L+DL;
                    L1 = max(L1,1);
                    L2 = min(L2,length(SGM));
                    token = zeros(size(SGM));
                    token(L1:L2)=1;
                    if 0~=sum(SGM&token)  
                            cnt_projected(A) = cnt_projected(A) + 1;
                    end
                    %
                    anno_tot = anno_tot|token;
                    %
                end
                %
                z_lab = bwlabel(SGM);
                tot_num_ROIs = max(z_lab);
                tot_time_ROIs = sum(0~=z_lab)/obj.Fs_ADC;
                tot_time = length(z_lab)/obj.Fs_ADC;
                %
                tot_annotators_time = sum(0~=anno_tot)/obj.Fs_ADC;
                tot_annotators_projected_ROIs = 0;
                for l=1:tot_num_ROIs
                    if 0~= sum( (z_lab==l) & anno_tot)
                        tot_annotators_projected_ROIs = tot_annotators_projected_ROIs + 1;
                    end
                end
                %
                % params stats
                params_data = [cellstr(obj.current_filename) num2cell([tot_time tot_time_ROIs tot_num_ROIs N_annotated cnt cnt_projected])];
                subj_indices = cell2mat(obj.ADC_trails_features_data(:,2));
                for k = 10:22
                    param = obj.ADC_trails_features_data(:,k);
                    sample = cell2mat(param(subj_indices==subj_ind));
                    v1 = mean(sample);
                    v2 = std(sample);
                    v3 = quantile(sample,.25);
                    v4 = median(sample);
                    v5 = quantile(sample,.75);
                    %
                    params_data = [params_data num2cell([v1 v2 v3 v4 v5])];
                end
                %                
                % STATS ON LATENCY INTERVALS - STARTS
                z_lab = bwlabel(~SGM);
                STATS = regionprops(z_lab,'Area');
                sample = cat(1, STATS.Area)/obj.Fs_ADC;
                    v1 = mean(sample);
                    v2 = std(sample);
                    v3 = quantile(sample,.25);
                    v4 = median(sample);
                    v5 = quantile(sample,.75);
                params_data = [params_data v1 v2 v3 v4 v5 ];                                    
                % STATS ON LATENCY INTERVALS - ENDS                                                                
                %
                % start %%%%%%% the % of the latent annotated time, also latent according to the detection
                % dil_SGM = imdilate(SGM,strel('line',round(2.5*obj.Fs_ADC),90)); % not sure if dilation is needed
                dil_SGM = SGM;
                prcntg_anno_time_latent_also_latent_by_sensor = sum( (~anno_tot) & (~dil_SGM) )/sum( ~anno_tot);
                % end %%%%%%% the % of the latent annotated time, also latent according to the detection                                                
                %                    
                % number of active channels
                num_ADC_channels = size(obj.current_data.ADC,2);
                n_active_channels = 0;
                n_used_channels = 0;
                for k = 1 : num_ADC_channels                    
                    if 0~=sum(obj.current_ADC_pre_processed(:,k))
                        n_active_channels = n_active_channels + 1;
                    end                        
                    if 0~=sum(obj.current_ADC_segmented(:,k))
                        n_used_channels = n_used_channels + 1;
                    end                                            
                end                                                
                % number of active channels
                %
                
                % number of "physician's annotations according to detection" - estimate
                dtct_annos = bwlabel(imdilate(SGM,strel('line',round(2.5*obj.Fs_ADC),90))); % dilate and count
                num_anno_according_to_detection = max(dtct_annos);
                                                
                % how nmany of these "physician's annotations according to detection" are out of annotations?
                n_intersect = 0;
                for l=1:num_anno_according_to_detection
                    if 0~=sum((dtct_annos==l).*anno_tot)
                        n_intersect = n_intersect+1;
                    end                    
                end
                num_anno_out = num_anno_according_to_detection - n_intersect;                
                %
                %... duration of ROIs only during the doctor’s annotated times
                % SGM and anno_tot
                duration_of_ROIs_only_during_the_doctors_annotated_times = sum(SGM & anno_tot)/obj.Fs_ADC;
                
                params_data = [params_data num2cell(prcntg_anno_time_latent_also_latent_by_sensor)];
                params_data = [params_data num2cell(n_active_channels)];
                params_data = [params_data num2cell(n_used_channels)];
                params_data = [params_data num2cell(tot_annotators_time)];
                params_data = [params_data num2cell(tot_annotators_projected_ROIs)];
                params_data = [params_data num2cell(num_anno_out)];
                params_data = [params_data num2cell(num_anno_according_to_detection)];
                params_data = [params_data num2cell(duration_of_ROIs_only_during_the_doctors_annotated_times)];
                
                num1 = tot_time - tot_time_ROIs;
                denom1 = tot_time;
                ratio1 = num1/denom1;
                %
                num2 = (sum( (~anno_tot) & (~dil_SGM) ))/obj.Fs_ADC;
                denom2 = (sum( ~anno_tot ))/obj.Fs_ADC;
                ratio2 = num2/denom2;
                
                logdata = [logdata; [num1 denom1 ratio1 num2 denom2 ratio2] ];
                
                % per subject pattern recognition for unknowns - starts
                %
                [groups,gr_nums_kNN] = do_per_subject_pattern_recognition_for_unknowns(obj,ALL_fv_data,ALL_annos,TR_data,TR_IDX,obj.current_filename,subj_ind,'PCA->kNN');
                [groups,gr_nums_QDA] = do_per_subject_pattern_recognition_for_unknowns(obj,ALL_fv_data,ALL_annos,TR_data,TR_IDX,obj.current_filename,subj_ind,'PCA->Quadratic');
                [groups,gr_nums_LDA] = do_per_subject_pattern_recognition_for_unknowns(obj,ALL_fv_data,ALL_annos,TR_data,TR_IDX,obj.current_filename,subj_ind,'PCA->Linear');
                if isempty(gr_nums_kNN)
                    gr_nums_kNN = zeros(size(groups));
                end
                if isempty(gr_nums_QDA)
                    gr_nums_QDA = zeros(size(groups));
                end
                if isempty(gr_nums_LDA)
                    gr_nums_LDA = zeros(size(groups));
                end                
                %
                % per subject pattern recognition for unknowns - ends
                params_data = [params_data num2cell([gr_nums_kNN gr_nums_QDA gr_nums_LDA])];
                
                stats = [stats; params_data];
                
            end
            if ~isempty(hw), delete(hw), drawnow; end;            
            
            names = [names {'mean_ltntT','std_ltntT','Q25_ltntT','median_ltntT','Q75_ltntT', ...
                'prcntg_anno_time_latent_also_latent_by_sensor', ...
                'n_actv_dtctrs','n_used_dtctrs', ...
                'tot_anno_time','num_anno_proj_ROIs','num_anno_out','num_anno_according_to_detection', ...
                'duration_of_ROIs_only_during_the_doctors_annotated_times','uncat_kNN_G','uncat_kNN_S','uncat_QDA_G','uncat_QDA_S','uncat_LDA_G','uncat_LDA_S'}];
            
            stats = [names; stats];

            % alternative outputs
            %stats = get_ROI_timing_data(obj);
            %stats = get_Physician_timing_data(obj);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    function [groups,gr_nums] = do_per_subject_pattern_recognition_for_unknowns(obj,ALL_fv_data,ALL_annos,TR_data,TR_IDX,subj_name,subj_ind,method)
        
        gr_nums = [];
        
        disp(subj_name);
        
        t1 = cell2mat(obj.ADC_trails_features_data(:,5));
        TST_data = [];
                             for kk = 1:length(t1)
                                 subj_index_k = cell2mat(obj.ADC_trails_features_data(kk,2)); % faster via index
                                 anno_ = ALL_annos(kk);
                                 if (subj_index_k == subj_ind) && 0==anno_
                                     TST_data = [TST_data; ALL_fv_data(kk,:)];                                     
                                 end                                                                                                   
                             end                                                    
        
        if ~isempty(TST_data)
            
            % PCA data & choose first 3 components
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
            % PCA data & choose first 3 components                                                         
            
                    try                                
                        switch method                
                            case {'PCA->Linear','PCA->Quadratic'}
                                if ~isempty(strfind(method,'Linear'))
                                    IDX2 = classify(TST_data,TR_data,TR_IDX,'Linear');
                                else
                                    IDX2 = classify(TST_data,TR_data,TR_IDX,'Quadratic');
                                end
                            case 'PCA->kNN'
                                mdl = fitcknn(TR_data,TR_IDX);
                                IDX2 = predict(mdl,TST_data);
                        end
                    catch
                        errordlg('classify error - possibly not enough data..');
                        return;
                    end
                    %
                    groups = sort(unique(TR_IDX));
                    gr_nums = zeros(1,length(groups));
                    for kk = 1:length(gr_nums)
                        gr_nums(kk) = sum(groups(kk)==IDX2(:));
                    end                    
                    %
                    disp(groups);
                    disp(gr_nums);
                    %                       
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    function outdata = get_Physician_timing_data(obj,~)
  
        outdata = [];
                    hw = waitbar(0,'calculating subject stats - please wait');
                    for subj_ind = 1:length(obj.subj_filenames)
                        if ~isempty(hw), waitbar(subj_ind/length(obj.subj_filenames),hw); drawnow, end;
                        obj.switch_current_to_subject(char(obj.subj_filenames(subj_ind)));
                        %
                        anno = obj.current_annotation;
                        anno_t = obj.current_annotation_time;
                        disp(numel(anno));
                        disp(numel(anno_t));
                        %
                        fname = char(obj.subj_filenames(subj_ind));
                        for k=1:numel(anno)
                            rec = [ {fname} num2cell(anno(k)) num2cell(anno_t(k)) ];
                            outdata = [outdata; rec];                    
                        end
                    end
                    if ~isempty(hw),delete(hw),end;

                    caption = {'filename','phys_anno','anno_time'};
                    outdata = [caption; outdata];
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    function outdata = get_ROI_timing_data(obj,~)
  
    if isempty(obj.ADC_trails_features_data), return, end;

    [Ndata, ~] = size(obj.ADC_trails_features_data);
    
    [fvD, annos] = obj.get_selected_feature_vector_data();
    
    unannotated_fv_data = [];
    unannotated_index = []; % index in the whole data
    
            training_fv_data = [];
            training_IDX = [];
            for k = 1:numel(annos)
                if 0~=annos(k)
                    training_fv_data = [training_fv_data; fvD(k,:)];
                    training_IDX = [training_IDX; annos(k)];
                else % unannotated
                    unannotated_fv_data = [unannotated_fv_data; fvD(k,:)];
                    unannotated_index = [unannotated_index k];                    
                end
            end
            
            % 3 times apply classifier to the unannotated data
            TR_data = training_fv_data;
            TST_data = unannotated_fv_data; 
            % PCA data & choose first 3 components
            truncation = 3;                                
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
            % PCA data & choose first 3 components                                                         
            
                    IDX_unannotated_LDA = classify(TST_data,TR_data,training_IDX,'Linear');                                
                    IDX_unannotated_QDA = classify(TST_data,TR_data,training_IDX,'Quadratic');                                    
                    mdl = fitcknn(TR_data,training_IDX);
                    IDX_unannotated_kNN = predict(mdl,TST_data);

                    anno_unannotated_LDA = zeros(Ndata,1);
                    anno_unannotated_QDA = zeros(Ndata,1);
                    anno_unannotated_kNN = zeros(Ndata,1);
                    for k=1:numel(IDX_unannotated_kNN)
                        ind = unannotated_index(k);
                        anno_unannotated_LDA(ind)=IDX_unannotated_LDA(k);
                        anno_unannotated_QDA(ind)=IDX_unannotated_QDA(k);
                        anno_unannotated_kNN(ind)=IDX_unannotated_kNN(k);
                    end
                    
                    fnames = obj.ADC_trails_features_data(:,1);
                    t1 = obj.ADC_trails_features_data(:,5);
                    t2 = obj.ADC_trails_features_data(:,6);
                    duration = obj.ADC_trails_features_data(:,10);
                    annotation = obj.ADC_trails_features_data(:,7);
                    
                    outdata = [fnames t1 t2 duration annotation num2cell(anno_unannotated_kNN) num2cell(anno_unannotated_LDA) num2cell(anno_unannotated_QDA)];
                    caption = {'filename','t1','t2','duration','anno_Phys','anno_kNN','anno_LDA','anno_QDA'};
                    %
                    outdata = [caption; outdata];
    end


end



    