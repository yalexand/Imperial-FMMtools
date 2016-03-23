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

% Last Modified by GUIDE v2.5 21-Mar-2016 14:17:54

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
set(handles.supervised_learning_type,'String',data_controller.supervised_learning_types);

set(handles.supervised_classification_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.features_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.unsupervised_clustering_pane, 'xticklabel', [], 'yticklabel', []);
set(handles.pairwise_group_comparison_axes, 'xticklabel', [], 'yticklabel', []);

set(handles.vis_original,'Value',true);
set(handles.vis_pre_processed,'Value',true);
set(handles.vis_segmented,'Value',true);

set(handles.pairwise_group_comparison_feature,'String',data_controller.ADC_fv_all);
set(handles.pairwise_group_comparison_item1,'String',data_controller.groups_selected);
set(handles.pairwise_group_comparison_item2,'String',data_controller.groups_selected);

set(handles.exclude_strong_IMU_checkbox,'Value',data_controller.exclude_IMU);

set(handles.confusion_matrix, 'Data', eye(6));
set(handles.confusion_matrix, 'RowName', data_controller.groups_all);
set(handles.confusion_matrix, 'ColumnName', data_controller.groups_all);

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
set(handles.corrX_chooser,'String',dc.ADC_feature_names(6:length(dc.ADC_feature_names)));
set(handles.corrY_chooser,'String',dc.ADC_feature_names(6:length(dc.ADC_feature_names)));

% --- Executes on button press in supervised_classification_go.
function supervised_classification_go_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
    supervised_learning_type_string = get(handles.supervised_learning_type,'String');
    index = get(handles.supervised_learning_type,'Value');
[sup_coords,sup_IDX] = dc.get_canons_for_supervised_classification(char(supervised_learning_type_string(index)));
handles.sup_coords = sup_coords;
handles.sup_IDX = sup_IDX;    
%
guidata(hObject, handles);
visualize_supervised_clustering(handles);

% --- Executes on button press in unsupervised_clustering_go.
function unsupervised_clustering_go_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
type_names = get(handles.unsupervised_clustering_type,'String');
type = char(type_names(get(handles.unsupervised_clustering_type,'Value')));
%
n_clusters = get(handles.n_clusters_chooser,'Value')+1; % should be OK..
[coords,IDX] = dc.perform_unsupervised_clustering(type,n_clusters);
handles.coords = coords;
handles.IDX = IDX;
guidata(hObject, handles);
visualize_unsupervised_clustering(handles);


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
unsupervised_clustering_go_Callback(hObject, eventdata, handles);

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
supervised_classification_go_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function feature_extraction_go_menu_Callback(hObject, eventdata, handles)
% hObject    handle to feature_extraction_go_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
feature_extraction_go_Callback(hObject, eventdata, handles);

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
pre_processing_type_index = get(handles.pre_processing_type,'Value');
segmentation_type_index = get(handles.segmentation_type,'Value');
pre_processing_type_string = get(handles.pre_processing_type,'String');
segmentation_type_string = get(handles.segmentation_type,'String');
%
dc = handles.data_controller;
res = dc.load_single_subject(pre_processing_type_string{pre_processing_type_index}, ...
    segmentation_type_string{segmentation_type_index});
if ~res, return, end;
%
update_record_pane(handles);
% clear panes
cla(handles.supervised_classification_pane,'reset');
cla(handles.features_pane,'reset');
cla(handles.unsupervised_clustering_pane,'reset');
clear_pairwise_comparison_visuals(handles);
%
set(handles.figure1,'Name',['FMMtools : ' num2str(dc.current_filename)]);
%
if isfield(handles,'coords')
    handles.coords = [];
    handles.sup_coords = [];
end
if isfield(handles,'IDX')
    handles.IDX = [];
    handles.sup_IDX = [];    
end

set(handles.subject_list,'Value',1);
set(handles.subject_list,'String',dc.current_filename);

setup_available_group_items(handles);

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

        org_plot_mode = 'k-';
        prp_plot_mode = 'b-';
        sgm_plot_mode = 'r-';

        % annotations
        b_plot_mode = 'gs';
        g_plot_mode = 'go';
        h_plot_mode = 'gx';
        l_plot_mode = 'gd';
        s_plot_mode = 'g*';     
    
        rec = [];
        prp = [];
        sgm = [];
        base = [];
        cap = [];
        
    if 1==rec_type_index
        
        rec = dc.current_data.ADC(:,channel_index);
        prp = dc.current_ADC_pre_processed(:,channel_index);
        
% not sure how to do it        
%         prp_log_abs = false; 
%         if prp_log_abs
%             avr_window = round(dc.ADC_segm_Moving_Average_Window*dc.Fs_ADC);            
%             [prp,~] = TD_high_pass_filter( prp, avr_window );                        
%             prp = log(abs(prp));
%             vals = prp(~isinf(prp));
%             minval = min(vals(:));
%             prp(isinf(prp)) = minval;
%         end
                
        %sgm = dc.current_ADC_segmented(:,channel_index)*0.05; %!!!                                
        
        % to display segmented signal..
        sgm = dc.current_ADC_segmented(:,channel_index);
        t = quantile(rec,0.1) ;
        z = rec(rec>t);
        base = median(z(:));
        cap = 3*std(z(:));
        sgm = ones(length(sgm),1)*base + sgm*cap;        
        
        %
        annos = dc.current_annotation;
        annos_times = dc.current_annotation_time;

        % 'bghls' - > '12345'
        b = (annos==1);
        t_b = annos_times(b);
        %
        g = (annos==2);
        t_g = annos_times(g);
        %
        h = (annos==3);
        t_h = annos_times(h);
        %
        l = (annos==4);
        t_l = annos_times(l);
        %
        s = (annos==5);
        t_s = annos_times(s);

        b=ones(size(t_b))*base;
        g=ones(size(t_g))*base;
        h=ones(size(t_h))*base;
        l=ones(size(t_l))*base;
        s=ones(size(t_s))*base;
        
    elseif 2==rec_type_index
        
        rec = dc.current_data.IMU(:,channel_index);
        prp = dc.current_IMU_pre_processed(:,channel_index);
                
        b = [];
        t_b = [];
        g = [];
        t_g = [];
        h = [];
        t_h = [];
        l = [];
        t_l = [];
        s = [];
        t_s = [];
                
        sgm = dc.current_IMU_segmented(:,channel_index);
        % to display segmented signal..
        t = quantile(prp,0.1) ;
        z = prp(prp>t);
        base = median(z(:));
        cap = 3*std(z(:));
        sgm = ones(length(sgm),1)*base + sgm*cap;             
    end
                             
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
            plot(handles.record_pane,t,rec,org_plot_mode,t,prp,prp_plot_mode,t,sgm,sgm_plot_mode, ... 
                t_b,b,b_plot_mode, ...
                t_g,g,g_plot_mode, ...
                t_h,h,h_plot_mode, ...
                t_l,l,l_plot_mode, ...
                t_s,s,s_plot_mode);
        elseif vis_org && ~vis_prp && vis_sgm
            YLIMS = [min([rec(:); sgm(:)]) max([rec(:); sgm(:)])];
            plot(handles.record_pane,t,rec,org_plot_mode,t,sgm,sgm_plot_mode, ...
                t_b,b,b_plot_mode, ...
                t_g,g,g_plot_mode, ...
                t_h,h,h_plot_mode, ...
                t_l,l,l_plot_mode, ...
                t_s,s,s_plot_mode);            
        elseif ~vis_org && vis_prp && vis_sgm
            YLIMS = [min([prp(:); sgm(:)]) max([prp(:); sgm(:)])];
            plot(handles.record_pane,t,prp,prp_plot_mode,t,sgm,sgm_plot_mode, ...
                t_b,b,b_plot_mode, ...
                t_g,g,g_plot_mode, ...
                t_h,h,h_plot_mode, ...
                t_l,l,l_plot_mode, ...
                t_s,s,s_plot_mode);
        elseif ~vis_org && vis_prp && ~vis_sgm
            YLIMS = [min(prp(:)) max(prp(:))];
            plot(handles.record_pane,t,prp,prp_plot_mode);
        elseif ~vis_org && ~vis_prp && vis_sgm
            YLIMS = [min(sgm(:)) max(sgm(:))];
            plot(handles.record_pane,t,sgm,sgm_plot_mode, ...
                t_b,b,b_plot_mode, ...
                t_g,g,g_plot_mode, ...
                t_h,h,h_plot_mode, ...
                t_l,l,l_plot_mode, ...
                t_s,s,s_plot_mode);
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
        offset = 6;                    
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
        
%-------------------------------------------------------------------------%        
function visualize_unsupervised_clustering(handles)

    try
        coords = handles.coords;
        IDX = handles.IDX;
    catch
        return;
    end
    if isempty(coords), return, end;
    
    v1 = coords(:,1);
    v2 = coords(:,2);
    try
    v3 = coords(:,3);
    catch
    end
    %
    color_1 = [1 0.2 0.4];
    color_2 = [0.34 0.65 0.87];
    color_3 = [0.5 0.5 0.5];
    color_4 = [0 0 1];
    color_5 = [0 1 0];
    color_6 = [1 1 0];    
    
    cmap = [color_1; color_2; color_3; color_4; color_5; color_6];
    IDX_color = cmap(IDX,:);
    %
    names = get(handles.unsupervised_clustering_vis_mode,'String');
    unsup_vis_mode = char(names(get(handles.unsupervised_clustering_vis_mode,'Value')));
    %        
        if      strcmp('3D',unsup_vis_mode) && exist('v3','var')
            scatter3(handles.unsupervised_clustering_pane,v1,v2,v3,50,IDX_color,'filled','MarkerEdgeColor','white');
        elseif  strcmp('2D',unsup_vis_mode)
            axes(handles.unsupervised_clustering_pane);
            gscatter(v1,v2,IDX_color);
            xlabel(handles.unsupervised_clustering_pane,'                   C1');
            ylabel(handles.unsupervised_clustering_pane,'                   C2');
            legend(handles.unsupervised_clustering_pane,'off');
        end
%-------------------------------------------------------------------------%        
function visualize_supervised_clustering(handles)
        
    try
        coords = handles.sup_coords;
        IDX = handles.sup_IDX;
    catch
        return;
    end
    if isempty(coords), return, end;
    
    v1 = coords(:,1);
    v2 = coords(:,2);
    try
        v3 = coords(:,3);
    catch
    end
    %
    color_1 = [1 0.2 0.4];
    color_2 = [0.34 0.65 0.87];
    color_3 = [0.5 0.5 0.5];
    color_4 = [0 0 1];
    color_5 = [0 1 0];
    color_6 = [1 1 0];    
    
%     cmap = [color_1; color_2; color_3; color_4; color_5; color_6];
%     IDX_color = cmap(IDX,:);
%     %
%     scatter3(handles.supervised_classification_pane,v1,v2,v3,50,IDX_color,'filled','MarkerEdgeColor','white');

    cmap = [color_1; color_2; color_3; color_4; color_5; color_6];
    IDX_color = cmap(IDX,:);
    %
    names = get(handles.supervised_clustering_vis_mode,'String');
    sup_vis_mode = char(names(get(handles.supervised_clustering_vis_mode,'Value')));
    %        
        if      strcmp('3D',sup_vis_mode) && exist('v3','var')
            scatter3(handles.supervised_classification_pane,v1,v2,v3,50,IDX_color,'filled','MarkerEdgeColor','white');
        elseif  strcmp('2D',sup_vis_mode)
            axes(handles.supervised_classification_pane);
            gscatter(v1,v2,IDX_color);
            xlabel(handles.supervised_classification_pane,'C1');
            ylabel(handles.supervised_classification_pane,'C2');
            %legend(handles.supervised_classification_pane,'off');            
        end
        % fix legends
        dc = handles.data_controller;
        group_indices = unique(IDX);
        LEGENDS = [];
        for k=1:numel(group_indices)
            LEGENDS = [LEGENDS dc.groups_all(group_indices(k))];
        end
        legend(handles.supervised_classification_pane,LEGENDS);
                        
%-------------------------------------------------------------------------%  
function setup_available_group_items(handles)
    dc = handles.data_controller;
    set(handles.pairwise_group_comparison_item1,'String',dc.groups_available);
    set(handles.pairwise_group_comparison_item2,'String',dc.groups_available); 
    
    set(handles.confusion_matrix, 'Data', eye(numel(dc.groups_selected)));
    set(handles.confusion_matrix, 'RowName', dc.groups_selected);
    set(handles.confusion_matrix, 'ColumnName', dc.groups_selected);
%-------------------------------------------------------------------------%      
    function clear_pairwise_comparison_visuals(handles)               
    cla(handles.pairwise_group_comparison_axes,'reset');
    grid(handles.pairwise_group_comparison_axes,'on');    
    set(handles.pairwise_comparison_stats,'String', ...
    {'N1 = ', ...
    'N2 = ', ...
    'm1 = ', ...
    'm2 = ', ...
    'std1 = ', ...
    'std2 = ', ...
    'd = ', ...
    'AUC = ', ...
    'p95_KS = ', ...
    'p95_rnk = '});



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


% --- Executes on selection change in unsupervised_clustering_vis_mode.
function unsupervised_clustering_vis_mode_Callback(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_vis_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns unsupervised_clustering_vis_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unsupervised_clustering_vis_mode
visualize_unsupervised_clustering(handles);

% --- Executes during object creation, after setting all properties.
function unsupervised_clustering_vis_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unsupervised_clustering_vis_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_multiple_subjects_Callback(hObject, eventdata, handles)
% hObject    handle to load_multiple_subjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pre_processing_type_index = get(handles.pre_processing_type,'Value');
segmentation_type_index = get(handles.segmentation_type,'Value');
pre_processing_type_string = get(handles.pre_processing_type,'String');
segmentation_type_string = get(handles.segmentation_type,'String');
%
dc = handles.data_controller;
res  = dc.load_multiple_subjects(pre_processing_type_string{pre_processing_type_index}, ...
                            segmentation_type_string{segmentation_type_index},true);
if ~res, return, end;
                        
update_record_pane(handles);
% clear panes
cla(handles.supervised_classification_pane,'reset');
cla(handles.features_pane,'reset');
cla(handles.unsupervised_clustering_pane,'reset');
clear_pairwise_comparison_visuals(handles);
%

%
set(handles.figure1,'Name',['FMMtools : ' num2str(dc.current_filename)]);
%
if isfield(handles,'coords')
    handles.coords = [];
    handles.sup_coords = [];
end
if isfield(handles,'IDX')
    handles.IDX = [];
    handles.sup_IDX = [];    
end

set(handles.subject_list,'String',dc.subj_filenames);

setup_available_group_items(handles);

guidata(hObject, handles);


% --- Executes on selection change in subject_list.
function subject_list_Callback(hObject, eventdata, handles)
% hObject    handle to subject_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subject_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subject_list
subject_list_string = get(handles.subject_list,'String');
data_filename = subject_list_string(get(handles.subject_list,'Value'));
dc = handles.data_controller;
dc.switch_current_to_subject(data_filename);
update_record_pane(handles);
%
set(handles.figure1,'Name',['FMMtools : ' num2str(dc.current_filename)]);
%

% --- Executes during object creation, after setting all properties.
function subject_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subject_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in supervised_clustering_vis_mode.
function supervised_clustering_vis_mode_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_clustering_vis_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns supervised_clustering_vis_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from supervised_clustering_vis_mode
visualize_supervised_clustering(handles);

% --- Executes during object creation, after setting all properties.
function supervised_clustering_vis_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to supervised_clustering_vis_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_fv_components_pushbutton.
function select_fv_components_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_fv_components_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
feature_vector_components_chooser(dc,hObject, eventdata, handles);

% --------------------------------------------------------------------
function save_training_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_training_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function select_feature_vector_Callback(hObject, eventdata, handles)
% hObject    handle to select_feature_vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
select_fv_components_pushbutton_Callback(hObject, eventdata, handles)


% --- Executes on button press in select_groups_pushbutton.
function select_groups_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_groups_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
groups_chooser(dc,hObject, eventdata, handles);





% --- Executes on selection change in pairwise_group_comparison_feature.
function pairwise_group_comparison_feature_Callback(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pairwise_group_comparison_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pairwise_group_comparison_feature
pairwise_group_comparison_go_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function pairwise_group_comparison_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pairwise_group_comparison_go.
function pairwise_group_comparison_go_Callback(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
type = 'annotator"s + segmentation'; % FOR NOW..
%
str = get(handles.pairwise_group_comparison_feature,'String');
feature_vector_name = str(get(handles.pairwise_group_comparison_feature,'Value'));
%
str = get(handles.pairwise_group_comparison_item1,'String');
group1 = str(get(handles.pairwise_group_comparison_item1,'Value'));
%
str = get(handles.pairwise_group_comparison_item2,'String');
group2 = str(get(handles.pairwise_group_comparison_item2,'Value'));
%
[X1,Y1,X2,Y2,z] = dc.get_pairwise_comparison(feature_vector_name,group1,group2,type);
if isempty(z), clear_pairwise_comparison_visuals(handles), return, end;

g1_plot_mode = 'c.-';
g2_plot_mode = 'm.-';
plot(handles.pairwise_group_comparison_axes,X1,Y1,g1_plot_mode,X2,Y2,g2_plot_mode);
grid(handles.pairwise_group_comparison_axes,'on');

set(handles.pairwise_comparison_stats,'String', ...
{['N1 = ', num2str(z.N1)], ...
['N2 = ', num2str(z.N2)], ...
['m1 = ', num2str(z.m1)], ...
['m2 = ', num2str(z.m2)], ...
['std1 = ', num2str(z.std1)], ...
['std2 = ', num2str(z.std2)], ...
['d = ', num2str(z.d)], ...
['AUC = ', num2str(z.AUC)], ...
['p95_KS = ', num2str(z.p_ks)], ...
['p95_rnk = ', num2str(z.p_rnk)]});
legend(handles.pairwise_group_comparison_axes,[group1, group2]);

% --- Executes on selection change in pairwise_group_comparison_item1.
function pairwise_group_comparison_item1_Callback(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_item1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pairwise_group_comparison_item1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pairwise_group_comparison_item1
pairwise_group_comparison_go_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function pairwise_group_comparison_item1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_item1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pairwise_group_comparison_item2.
function pairwise_group_comparison_item2_Callback(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_item2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pairwise_group_comparison_item2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pairwise_group_comparison_item2
pairwise_group_comparison_go_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function pairwise_group_comparison_item2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pairwise_group_comparison_item2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in supervised_classification_type.
function supervised_classification_type_Callback(hObject, eventdata, handles)
% hObject    handle to supervised_classification_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns supervised_classification_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from supervised_classification_type


% --- Executes during object creation, after setting all properties.
function supervised_classification_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to supervised_classification_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in confusion_matrix_go_pushbutton.
function confusion_matrix_go_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to confusion_matrix_go_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
[groups,CM] = dc.calculate_confusion_matrix('annotator"s + segmentation','Linear');
if ~isempty(CM)
    set(handles.confusion_matrix, 'Data',CM);
    set(handles.confusion_matrix, 'RowName', dc.groups_all(groups)); % must be dc_groups_selected
    set(handles.confusion_matrix, 'ColumnName', dc.groups_all(groups));    
end

% --- Executes on selection change in pairwise_comparison_stats.
function pairwise_comparison_stats_Callback(hObject, eventdata, handles)
% hObject    handle to pairwise_comparison_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pairwise_comparison_stats contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pairwise_comparison_stats


% --- Executes during object creation, after setting all properties.
function pairwise_comparison_stats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pairwise_comparison_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
                
end


% --- Executes on button press in exclude_strong_IMU_checkbox.
function exclude_strong_IMU_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to exclude_strong_IMU_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exclude_strong_IMU_checkbox
dc = handles.data_controller;
dc.exclude_IMU = get(handles.exclude_strong_IMU_checkbox,'Value');
