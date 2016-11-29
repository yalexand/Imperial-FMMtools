function stats = get_subject_feature_stats(obj,~,~)

stats = [];

if isempty(obj.ADC_trails_features_data), return, end;
            
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

                stats = [stats; params_data];
                
                num1 = tot_time - tot_time_ROIs;
                denom1 = tot_time;
                ratio1 = num1/denom1;
                %
                num2 = (sum( (~anno_tot) & (~dil_SGM) ))/obj.Fs_ADC;
                denom2 = (sum( ~anno_tot))/obj.Fs_ADC;
                ratio2 = num2/denom2;
                
                logdata = [logdata; [num1 denom1 ratio1 num2 denom2 ratio2] ];
                
            end
            if ~isempty(hw), delete(hw), drawnow; end;            
            
            names = [names {'mean_ltntT','std_ltntT','Q25_ltntT','median_ltntT','Q75_ltntT', ...
                'prcntg_anno_time_latent_also_latent_by_sensor', ...
                'n_actv_dtctrs','n_used_dtctrs', ...
                'tot_anno_time','num_anno_proj_ROIs','num_anno_out','num_anno_according_to_detection', ...
                'duration_of_ROIs_only_during_the_doctors_annotated_times'}];            
            
            stats = [names; stats];
end
