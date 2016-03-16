function varargout = groups_chooser(varargin)
% GROUPS_CHOOSER MATLAB code for groups_chooser.fig
%      GROUPS_CHOOSER, by itself, creates a new GROUPS_CHOOSER or raises the existing
%      singleton*.
%
%      H = GROUPS_CHOOSER returns the handle to a new GROUPS_CHOOSER or the handle to
%      the existing singleton*.
%
%      GROUPS_CHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROUPS_CHOOSER.M with the given input arguments.
%
%      GROUPS_CHOOSER('Property','Value',...) creates a new GROUPS_CHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before groups_chooser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to groups_chooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help groups_chooser

% Last Modified by GUIDE v2.5 16-Mar-2016 11:57:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @groups_chooser_OpeningFcn, ...
                   'gui_OutputFcn',  @groups_chooser_OutputFcn, ...
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


% --- Executes just before groups_chooser is made visible.
function groups_chooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to groups_chooser (see VARARGIN)

% Choose default command line output for groups_chooser
handles.output = hObject;

dc = varargin{1};
if isempty(dc.groups_available), 
    errordlg('error while loading the list of available groups'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

set(handles.all_groups_list,'String',dc.groups_available);
set(handles.all_groups_list,'min',1);
set(handles.all_groups_list,'max',length(dc.groups_available));

set(handles.selected_groups_list,'String',dc.groups_selected);

handles.data_controller = dc;

handles.previous_selection = dc.groups_selected;


%%%%%%%%%%%%%%
%@(hObject,eventdata)FMMtools_GUIDE_main('load_single_subject_menu_Callback',hObject,eventdata,guidata(hObject))
%%%%%%%%%%%%%%
caller_hObject = varargin{2};
caller_eventdata = varargin{3};
caller_handles = varargin{4};

handles.caller_hObject = caller_hObject;
handles.caller_eventdata = caller_eventdata;
handles.caller_handles = caller_handles;
%@(hObject,eventdata)FMMtools_GUIDE_main('load_single_subject_menu_Callback',hObject,eventdata,guidata(hObject))
%%%%%%%%%%%%%%

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes groups_chooser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = groups_chooser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in all_groups_list.
function all_groups_list_Callback(hObject, eventdata, handles)
% hObject    handle to all_groups_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns all_groups_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from all_groups_list
handles.selection = get(handles.all_groups_list,'Value');
%disp(get(handles.all_groups_list,'Value'));
str = get(handles.all_groups_list,'String');
set(handles.selected_groups_list,'String',str(handles.selection));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function all_groups_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to all_groups_list (see GCBO)
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


% --- Executes on selection change in selected_groups_list.
function selected_groups_list_Callback(hObject, eventdata, handles)
% hObject    handle to selected_groups_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selected_groups_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selected_groups_list


% --- Executes during object creation, after setting all properties.
function selected_groups_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_groups_list (see GCBO)
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
    errordlg('wrong selection, groups won"t change');
end

if isfield(handles,'selection') && ~ isempty(handles.selection) && numel(handles.selection) >= 2
    dc.groups_selected = dc.groups_available(handles.selection);
else    
    dc.groups_selected = handles.previous_selection;
end

caller_hObject = handles.caller_hObject;
caller_eventdata = handles.caller_eventdata;
FMMtools_GUIDE_main('supervised_classification_go_Callback',caller_hObject,caller_eventdata,guidata(caller_hObject));

fh = ancestor(hObject,'figure');     
delete(fh);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

