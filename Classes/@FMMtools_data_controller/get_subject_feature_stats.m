function stats = get_subject_feature_stats(obj,~,~)

stats = [];

if isempty(obj.ADC_trails_features_data), return, end;
            
            d2 = obj.annotators_delay;
            d1 = obj.annotators_reaction; % [second] - minimal discernable time between event and annotation
            
            hw = waitbar(0,'calculating subject stats - please wait');
            for subj_ind = 1:length(obj.subj_filenames)
                if ~isempty(hw), waitbar(subj_ind/length(obj.subj_filenames),hw); drawnow, end;
                obj.switch_current_to_subject(char(obj.subj_filenames(subj_ind)));

                % one segmentation for all
                num_ADC_channels = size(obj.current_data.ADC,2);
                SGM = obj.current_ADC_segmented(:,1);
                for k = 2 : num_ADC_channels
                    SGM = SGM | obj.current_ADC_segmented(:,k);
                end                
                                
                anno = obj.subj_data(subj_ind).annotation;
                anno_t = obj.subj_data(subj_ind).annotation_time; 
                N_annotated = length(anno);
                %
                cnt = zeros(1,6);
                cnt_projected = zeros(1,6);
                for k =1:length(anno)
                    A = anno(k);
                    cnt(A) = cnt(A)+1;
                    %
                    T1 = anno_t(k) - d1 - d2;
                    T2 = anno_t(k) - d1;
                    L = round((T1+T2)/2*obj.Fs_ADC);
                    if L>=1 && L<=length(SGM) 
                        if 0~= SGM(L)
                            cnt_projected(A) = cnt_projected(A) + 1;
                        end
                    else
                        disp([cellstr(obj.current_filename) num2cell([subj_ind L length(SGM)])]);
                    end
                end                    
                %
                z_lab = bwlabel(SGM);
                tot_num_ROIs = max(z_lab);
                tot_time_ROIs = sum(0~=z_lab)/obj.Fs_ADC;
                tot_time = length(z_lab)/obj.Fs_ADC;
                %
                % params stats
                params_data = [cellstr(obj.current_filename) num2cell([tot_time tot_time_ROIs tot_num_ROIs N_annotated cnt cnt_projected])];
                subj_indices = cell2mat(obj.ADC_trails_features_data(:,2));
                for k = 7:19
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
                stats = [stats; params_data];
            end
            if ~isempty(hw), delete(hw), drawnow; end;

end
