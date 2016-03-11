function varargout = feature_vector_components_chooser(varargin)
% FEATURE_VECTOR_COMPONENTS_CHOOSER MATLAB code for feature_vector_components_chooser.fig
%      FEATURE_VECTOR_COMPONENTS_CHOOSER, by itself, creates a new FEATURE_VECTOR_COMPONENTS_CHOOSER or raises the existing
%      singleton*.
%
%      H = FEATURE_VECTOR_COMPONENTS_CHOOSER returns the handle to a new FEATURE_VECTOR_COMPONENTS_CHOOSER or the handle to
%      the existing singleton*.
%
%      FEATURE_VECTOR_COMPONENTS_CHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEATURE_VECTOR_COMPONENTS_CHOOSER.M with the given input arguments.
%
%      FEATURE_VECTOR_COMPONENTS_CHOOSER('Property','Value',...) creates a new FEATURE_VECTOR_COMPONENTS_CHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before feature_vector_components_chooser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to feature_vector_components_chooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help feature_vector_components_chooser

% Last Modified by GUIDE v2.5 11-Mar-2016 12:33:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @feature_vector_components_chooser_OpeningFcn, ...
                   'gui_OutputFcn',  @feature_vector_components_chooser_OutputFcn, ...
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


% --- Executes just before feature_vector_components_chooser is made visible.
function feature_vector_components_chooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to feature_vector_components_chooser (see VARARGIN)

% Choose default command line output for feature_vector_components_chooser
handles.output = hObject;

dc = varargin{1};
if isempty(dc.ADC_fv_all), 
    errordlg('error while loading the list of available features'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

set(handles.all_features_list,'String',dc.ADC_fv_all);
set(handles.all_features_list,'min',1);
set(handles.all_features_list,'max',length(dc.ADC_fv_all));

set(handles.selected_features_list,'String',dc.ADC_fv_selected);

handles.data_controller = dc;

handles.previous_selection = dc.ADC_fv_selected;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes feature_vector_components_chooser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = feature_vector_components_chooser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in all_features_list.
function all_features_list_Callback(hObject, eventdata, handles)
% hObject    handle to all_features_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns all_features_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from all_features_list
handles.selection = get(handles.all_features_list,'Value');
%disp(get(handles.all_features_list,'Value'));
str = get(handles.all_features_list,'String');
set(handles.selected_features_list,'String',str(handles.selection));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function all_features_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to all_features_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_pushbutton.
function select_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in selected_features_list.
function selected_features_list_Callback(hObject, eventdata, handles)
% hObject    handle to selected_features_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selected_features_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selected_features_list


% --- Executes during object creation, after setting all properties.
function selected_features_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_features_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OK_pushbutton.
function OK_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OK_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;

if isfield(handles,'selection') && numel(handles.selection) < 2
    errordlg('wrong selection, feature vector won"t change');
end

if isfield(handles,'selection') && ~ isempty(handles.selection) && numel(handles.selection) >= 2
    dc.ADC_fv_selected = dc.ADC_fv_all(handles.selection);
else    
    dc.ADC_fv_selected = handles.previous_selection;
end

fh = ancestor(hObject,'figure');     
delete(fh);
