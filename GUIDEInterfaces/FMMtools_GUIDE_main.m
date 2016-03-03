function varargout = FMMtools_GUIDE_main(varargin)
% FMMTOOLS_GUIDE_MAIN MATLAB code for FMMtools_GUIDE_main.fig
%      FMMTOOLS_GUIDE_MAIN, by itself, creates a new FMMTOOLS_GUIDE_MAIN or raises the existing
%      singleton*.
%
%      H = FMMTOOLS_GUIDE_MAIN returns the handle to a new FMMTOOLS_GUIDE_MAIN or the handle to
%      the existing singleton*.
%
%      FMMTOOLS_GUIDE_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FMMTOOLS_GUIDE_MAIN.M with the given input arguments.
%
%      FMMTOOLS_GUIDE_MAIN('Property','Value',...) creates a new FMMTOOLS_GUIDE_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FMMtools_GUIDE_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FMMtools_GUIDE_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FMMtools_GUIDE_main

% Last Modified by GUIDE v2.5 03-Mar-2016 15:58:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FMMtools_GUIDE_main_OpeningFcn, ...
                   'gui_OutputFcn',  @FMMtools_GUIDE_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FMMtools_GUIDE_main is made visible.
function FMMtools_GUIDE_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FMMtools_GUIDE_main (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

set(handles.pre_processing_type,'String',data_controller.preprocessing_types);
set(handles.segmentation_type,'String',data_controller.segmentation_types);

set(handles.supervised_classification_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.features_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.unsupervised_clustering_pane, 'xticklabel', [], 'yticklabel', []);

% Choose default command line output for FMMtools_GUIDE_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FMMtools_GUIDE_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FMMtools_GUIDE_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in record_type.
function record_type_Callback(hObject, eventdata, handles)
% hObject    handle to record_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_record_pane(handles);

% Hints: contents = cellstr(get(hObject,'String')) returns record_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from record_type


% --- Executes during object creation, after setting all properties.
function record_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in vis_original.
function vis_original_Callback(hObject, eventdata, handles)
% hObject    handle to vis_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vis_original
update_record_pane(handles);

% --- Executes on button press in vis_pre_processed.
function vis_pre_processed_Callback(hObject, eventdata, handles)
% hObject    handle to vis_pre_processed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vis_pre_processed
update_record_pane(handles);

% --- Executes on button press in vis_segmented.
function vis_segmented_Callback(hObject, eventdata, handles)
% hObject    handle to vis_segmented (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vis_segmented
update_record_pane(handles);

% --- Executes on selection change in record_ID.
function record_ID_Callback(hObject, eventdata, handles)
% hObject    handle to record_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_record_pane(handles);

% Hints: contents = cellstr(get(hObject,'String')) returns record_ID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from record_ID


% --- Executes during object creation, after setting all properties.
function record_ID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to record_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in time_frequency.
function time_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to time_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_record_pane(handles);

% Hints: contents = cellstr(get(hObject,'String')) returns time_frequency contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_frequency


% --- Executes during object creation, after setting all properties.
function time_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pre_processing_type.
function pre_processing_type_Callback(hObject, eventdata, handles)
% hObject    handle to pre_processing_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pre_processing_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pre_processing_type


% --- Executes during object creation, after setting all properties.
function pre_processing_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pre_processing_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segmentation_type.
function segmentation_type_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segmentation_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segmentation_type


% --- Executes during object creation, after setting all properties.
function segmentation_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segmentation_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature_extraction_type.
function feature_extraction_type_Callback(hObject, eventdata, handles)
% hObject    handle to feature_extraction_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns feature_extraction_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature_extraction_type


% --- Executes during object creation, after setting all properties.
function feature_extraction_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature_extraction_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in supervised_learning_type.
function supervised_learning_type_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_learning_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns supervised_learning_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from supervised_learning_type


% --- Executes during object creation, after setting all properties.
function supervised_learning_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to supervised_learning_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in unsupervised_clustering_type.
function unsupervised_clustering_type_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns unsupervised_clustering_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unsupervised_clustering_type


% --- Executes during object creation, after setting all properties.
function unsupervised_clustering_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in feature_extraction_go.
function feature_extraction_go_Callback(hObject, eventdata, handles)
% hObject    handle to feature_extraction_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
dc = handles.data_controller;
dc.ADC_trails_features_data = [];
dc.ADC_feature_names = [];
[dc.ADC_trails_features_data,dc.ADC_feature_names] = dc.extract_features_current_ADC;
guidata(hObject, handles);
%
visualize_current_ADC_trails_features_data(handles);
%
set(handles.corrX_chooser,'String',dc.ADC_feature_names(5:length(dc.ADC_feature_names)));
set(handles.corrY_chooser,'String',dc.ADC_feature_names(5:length(dc.ADC_feature_names)));

% --- Executes on button press in supervised_classification_go.
function supervised_classification_go_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in unsupervised_clustering_go.
function unsupervised_clustering_go_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
type_names = get(handles.unsupervised_clustering_type,'String');
type = char(type_names(get(handles.unsupervised_clustering_type,'Value')));

n_clusters = get(handles.n_clusters_chooser,'Value')+1; % should be OK..

[coords,IDX] = dc.perform_unsupervised_clustering(type,n_clusters);
     v1 = coords(:,1);
     v2 = coords(:,2);
     v3 = coords(:,3);
    color_1 = [1 0.2 0.4];
    color_2 = [0.34 0.65 0.87];
    color_3 = [0.5 0.5 0.5];
    cmap = [color_1; color_2; color_3];
    IDX_color = cmap(IDX,:);     
%
% visualize - 3D - ??
%     scatter3(handles.unsupervised_clustering_pane,v1,v2,v3,50,IDX_color,'filled','MarkerEdgeColor','white');
%
% visualize - 2D
    gscatter(v1,v2,IDX_color);
    xlabel(handles.unsupervised_clustering_pane,'C1');
    ylabel(handles.unsupervised_clustering_pane,'C2');
    

% --- Executes on button press in pre_processing_setups.
function pre_processing_setups_Callback(hObject, eventdata, handles)
% hObject    handle to pre_processing_setups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in segmentation_setups.
function segmentation_setups_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_setups (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function pre_processing_menu_Callback(hObject, eventdata, handles)
% hObject    handle to pre_processing_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function feature_extraction_menu_Callback(hObject, eventdata, handles)
% hObject    handle to feature_extraction_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function supervised_classification_menu_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function unsupervised_clustering_menu_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function unsupervised_clustering_set_number_of_clusters_menu_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_set_number_of_clusters_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function unsupervised_clustering_go_menu_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_go_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function supervised_classification_load_training_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_load_training_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function supervised_classification_go_menu_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_go_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function feature_extraction_go_on_curernt_menu_Callback(hObject, eventdata, handles)
% hObject    handle to feature_extraction_go_on_curernt_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function preprocessing_go_on_current_menu_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessing_go_on_current_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_single_subject_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_single_subject_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
data_controller.load_single_subject(true);

pre_processing_type_index = get(handles.pre_processing_type,'Value');
segmentation_type_index = get(handles.segmentation_type,'Value');
pre_processing_type_string = get(handles.pre_processing_type,'String');
segmentation_type_string = get(handles.segmentation_type,'String');
%
dc = handles.data_controller;
dc.pre_process_ADC(pre_processing_type_string{pre_processing_type_index});
dc.segment_ADC(segmentation_type_string{segmentation_type_index});
%
dc.pre_process_IMU([]); % that should be done better (?) - for now returns zeros
dc.segment_IMU([]);
%
dc.calculate_PSDs;
%
update_record_pane(handles);
% clear panes
cla(handles.supervised_classification_pane,'reset');
cla(handles.features_pane,'reset');
cla(handles.unsupervised_clustering_pane,'reset');
%
set(handles.figure1,'Name',['FMMtools : ' num2str(dc.current_filename)]);
guidata(hObject, handles);

% --------------------------------------------------------------------
function save_analysis_state_menu_Callback(hObject, eventdata, handles)
% hObject    handle to save_analysis_state_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%-------------------------------------------------------------------------%    
function update_record_pane(handles)

    dc = handles.data_controller;

    if isempty(dc.current_data), return, end;
    
    rec_type_index = get(handles.record_type,'Value');
    channel_index = get(handles.record_ID,'Value');
    plot_type_index = get(handles.time_frequency,'Value');
        
    % visualizer checkboxes states
    vis_org = get(handles.vis_original,'Value');
    vis_prp = get(handles.vis_pre_processed,'Value');
    vis_sgm = get(handles.vis_segmented,'Value');
        
    num_ADC_channels = size(dc.current_data.ADC,2);
    num_IMU_channels = size(dc.current_data.IMU,2);
        
    if 1==rec_type_index, %ADC
        if channel_index > num_ADC_channels, cla(handles.record_pane); return, end;
        Fs = dc.Fs_ADC;
    elseif 2==rec_type_index, %IMU
        if channel_index > num_IMU_channels, cla(handles.record_pane); return, end;
        Fs = dc.Fs_IMU;
    end

        rec = [];
        prp = [];
        sgm = [];    
    if 1==rec_type_index
        rec = dc.current_data.ADC(:,channel_index);
        prp = dc.current_ADC_pre_processed(:,channel_index);
        sgm = dc.current_ADC_segmented(:,channel_index)*0.05; %!!!                                
        
    elseif 2==rec_type_index
        rec = dc.current_data.IMU(:,channel_index);
        prp = dc.current_IMU_pre_processed(:,channel_index);
        sgm = dc.current_IMU_segmented(:,channel_index)*0.05; %!!!
    end
        
    org_plot_mode = 'k-';
    prp_plot_mode = 'b-';
    sgm_plot_mode = 'r-';
        
    if 1==plot_type_index
        t =(1:length(rec))/Fs;
                                
        XLIMS = [min(t(:)) max(t(:))];
        YLIMS = [min(rec(:)) max(rec(:))];
        
        % 2^3  = 8 combinations
        if vis_org && ~vis_prp && ~vis_sgm || (~vis_org && ~vis_prp && ~vis_sgm)
            plot(handles.record_pane,t,rec,org_plot_mode);
        elseif vis_org && vis_prp && ~vis_sgm
            YLIMS = [min([rec(:); prp(:)]) max([rec(:); prp(:)])];
            plot(handles.record_pane,t,rec,org_plot_mode,t,prp,prp_plot_mode);
        elseif vis_org && vis_prp && vis_sgm
            YLIMS = [min([rec(:); prp(:); sgm(:)]) max([rec(:); prp(:); sgm(:)])];
            plot(handles.record_pane,t,rec,org_plot_mode,t,prp,prp_plot_mode,t,sgm,sgm_plot_mode);
        elseif vis_org && ~vis_prp && vis_sgm
            YLIMS = [min([rec(:); sgm(:)]) max([rec(:); sgm(:)])];
            plot(handles.record_pane,t,rec,org_plot_mode,t,sgm,sgm_plot_mode);            
        elseif ~vis_org && vis_prp && vis_sgm
            YLIMS = [min([prp(:); sgm(:)]) max([prp(:); sgm(:)])];
            plot(handles.record_pane,t,prp,prp_plot_mode,t,sgm,sgm_plot_mode);
        elseif ~vis_org && vis_prp && ~vis_sgm
            YLIMS = [min(prp(:)) max(prp(:))];
            plot(handles.record_pane,t,prp,prp_plot_mode);
        elseif ~vis_org && ~vis_prp && vis_sgm
            YLIMS = [min(sgm(:)) max(sgm(:))];
            plot(handles.record_pane,t,sgm,sgm_plot_mode);            
        end
                    
        grid(handles.record_pane,'on');
        if 0~=sum(rec(:))
            set(handles.record_pane,'XLim',XLIMS);
            if YLIMS(1) < YLIMS(2)
                set(handles.record_pane,'YLim',YLIMS);
            end
        end
        xlabel(handles.record_pane,'[sec]');
        
    elseif 2==plot_type_index % power spectrum density

        if isempty(dc.current_ADC_PSD), return, end;

            psd = [];
            psd_prp = [];
            omega = [];
        
        if 1==rec_type_index
            psd = dc.current_ADC_PSD(:,channel_index);
            psd_prp = dc.current_ADC_preprocessed_PSD(:,channel_index);
            omega = dc.current_omega_ADC;
        elseif 2==rec_type_index
            psd = dc.current_IMU_PSD(:,channel_index);
            psd_prp = dc.current_IMU_preprocessed_PSD(:,channel_index);
            omega = dc.current_omega_IMU;
        end

        XLIMS = [min(omega(:)) max(omega(:))];
        YLIMS = [min(psd(:)) max(psd(:))];
        
        % 2^2  = 4 combinations
        if vis_org && ~vis_prp || (~vis_org && ~vis_prp)
            semilogy(handles.record_pane,omega,psd,org_plot_mode);
        elseif vis_org && vis_prp
            YLIMS = [min([psd(:); psd_prp(:)]) max([psd(:); psd_prp(:)])];
            semilogy(handles.record_pane,omega,psd,org_plot_mode,omega,psd_prp,prp_plot_mode);
        elseif ~vis_org && vis_prp
            YLIMS = [min(psd_prp(:)); max(psd_prp(:))];
            semilogy(handles.record_pane,omega,psd_prp,prp_plot_mode);
        end        

        grid(handles.record_pane,'on');
        if 0~=sum(rec(:))
            set(handles.record_pane,'XLim',XLIMS);
            if YLIMS(1) < YLIMS(2)
                set(handles.record_pane,'YLim',YLIMS);
            end
        end
        xlabel(handles.record_pane,'[Hz]');
                
    end
 %-------------------------------------------------------------------------%       
    function visualize_current_ADC_trails_features_data(handles) 

        dc = handles.data_controller;
        
        if isempty(dc.ADC_trails_features_data), return, end;
                
        piX = get(handles.corrX_chooser,'Value');
        piY = get(handles.corrY_chooser,'Value');        
        offset = 5;                    
        u1 = cell2mat(dc.ADC_trails_features_data(:,offset+piX-1));
        u2 = cell2mat(dc.ADC_trails_features_data(:,offset+piY-1));
                
        corrmap = correlation_map(u1,u2,dc.corr_map_W);
        
        AXES = handles.features_pane;
        imagesc(corrmap,'Parent',AXES);
        daspect(AXES,[1 1 1]);
        set(AXES, 'xticklabel', [], 'yticklabel', []);
        %
        param_names = get(handles.corrX_chooser,'String');
        xlabel(AXES,char(param_names(piX)));
        ylabel(AXES,char(param_names(piY)));

        colorbar('peer',AXES);
        
        

% --- Executes on selection change in corrX_chooser.
function corrX_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to corrX_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns corrX_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from corrX_chooser
visualize_current_ADC_trails_features_data(handles);

% --- Executes during object creation, after setting all properties.
function corrX_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrX_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in corrY_chooser.
function corrY_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to corrY_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns corrY_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from corrY_chooser
visualize_current_ADC_trails_features_data(handles);

% --- Executes during object creation, after setting all properties.
function corrY_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrY_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n_clusters_chooser.
function n_clusters_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to n_clusters_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n_clusters_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n_clusters_chooser
unsupervised_clustering_go_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function n_clusters_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_clusters_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
