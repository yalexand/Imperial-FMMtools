function stats = get_subject_feature_stats(obj,~,~)

stats = [];

if isempty(obj.ADC_trails_features_data), return, end;
            
            hw = waitbar(0,'calculating subject stats - please wait');
            for subj_ind = 1:length(obj.subj_filenames)
                if ~isempty(hw), waitbar(subj_ind/length(obj.subj_filenames),hw); drawnow, end;
                obj.switch_current_to_subject(char(obj.subj_filenames(subj_ind)));

                num_ADC_channels = size(obj.current_data.ADC,2);        
                %
                % one segmentation for all
                SGM = obj.current_ADC_segmented(:,1);
                for k = 2 : num_ADC_channels
                    SGM = SGM | obj.current_ADC_segmented(:,k);
                end                
                z_lab = bwlabel(SGM);
                tot_num_ROIs = max(z_lab);
                tot_time_ROIs = sum(0~=z_lab)/obj.Fs_ADC;
                %
                % params stats
                params_data = [cellstr(obj.current_filename) num2cell([tot_num_ROIs tot_time_ROIs])];
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
