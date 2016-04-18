function varargout = ROI_inspector(varargin)
% ROI_INSPECTOR MATLAB code for ROI_inspector.fig
%      ROI_INSPECTOR, by itself, creates a new ROI_INSPECTOR or raises the existing
%      singleton*.
%
%      H = ROI_INSPECTOR returns the handle to a new ROI_INSPECTOR or the handle to
%      the existing singleton*.
%
%      ROI_INSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_INSPECTOR.M with the given input arguments.
%
%      ROI_INSPECTOR('Property','Value',...) creates a new ROI_INSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROI_inspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROI_inspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_inspector

% Last Modified by GUIDE v2.5 18-Apr-2016 15:23:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROI_inspector_OpeningFcn, ...
                   'gui_OutputFcn',  @ROI_inspector_OutputFcn, ...
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


% --- Executes just before ROI_inspector is made visible.
function ROI_inspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROI_inspector (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

INDEX = 1;
handles.INDEX = INDEX;

set(handles.ROI_index_slider,'Min',1);
set(handles.ROI_index_slider,'Max',size(data_controller.ADC_trails_features_data,1));
set(handles.ROI_index_slider,'Value',1);

display_current_ROI(handles);

% Choose default command line output for ROI_inspector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ROI_inspector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ROI_inspector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in close_button.
function close_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on slider movement.
function ROI_index_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_index_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.INDEX = fix(get(hObject,'Value'));
display_current_ROI(handles);


% --- Executes during object creation, after setting all properties.
function ROI_index_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_index_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function display_current_ROI(handles)
dc = handles.data_controller;
INDEX = handles.INDEX;

data = dc.ADC_trails_features_data(INDEX,:)';
handles.INDEX = INDEX;
%
set(handles.features_table, 'Data',data);
set(handles.features_table, 'RowName', dc.ADC_feature_names); % must be dc_groups_selected
%
subj_ind = cell2mat(data(2));
detector_ind = cell2mat(data(3));
t1 = cell2mat(data(5));
t2 = cell2mat(data(6));
l1 = fix(t1*dc.Fs_ADC); 
    if 0==l1, l1=1; end;
l2 = fix(t2*dc.Fs_ADC);
%
subj_data = dc.subj_data(subj_ind);
%
t = (l1:l2)/dc.Fs_ADC;
minT = min(t(:));
maxT = max(t(:));

cla(handles.signal_axes,'reset');

minY = Inf;
maxY = -Inf;
%
hold(handles.signal_axes,'on');
for k = 1 : 8
        signal = subj_data.ADC_pre_processed(:,k); 
        S = signal(l1:l2)';
        
        if k == detector_ind
            plot(handles.signal_axes,t,S,'k.-','LineWidth',2);
        else
            plot(handles.signal_axes,t,S);
        end
        
        minYk = min(S(:));
        maxYk = max(S(:));
        if minYk < minY && 0~=minYk, minY = minYk; end;
        if maxYk > maxY, maxY = maxYk; end;
end
hold(handles.signal_axes,'off');
%
set(handles.signal_axes,'XLim',[minT maxT]);
set(handles.signal_axes,'YLim',[minY maxY]);        
%
grid(handles.signal_axes,'on');
legend(handles.signal_axes,{'1','2','3','4','5','6','7','8'});
