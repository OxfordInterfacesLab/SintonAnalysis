function varargout = LifetimeAnalyser_3(varargin)
% PLOTLIFETIMES_GUI_3 MATLAB code for PlotLifetimes_GUI_3.fig
%      PLOTLIFETIMES_GUI_3, by itself, creates a new PLOTLIFETIMES_GUI_3 or raises the existing
%      singleton*.
%
%      H = PLOTLIFETIMES_GUI_3 returns the handle to a new PLOTLIFETIMES_GUI_3 or the handle to
%      the existing singleton*.
%
%      PLOTLIFETIMES_GUI_3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTLIFETIMES_GUI_3.M with the given input arguments.
%
%      PLOTLIFETIMES_GUI_3('Property','Value',...) creates a new PLOTLIFETIMES_GUI_3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotLifetimes_GUI_3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotLifetimes_GUI_3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotLifetimes_GUI_3

% Last Modified by GUIDE v2.5 11-Apr-2018 16:39:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotLifetimes_GUI_3_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotLifetimes_GUI_3_OutputFcn, ...
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

% --- Executes just before PlotLifetimes_GUI_3 is made visible.
function PlotLifetimes_GUI_3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotLifetimes_GUI_3 (see VARARGIN)
set(0,'DefaultFigureWindowStyle','normal')
% Choose default command line output for PlotLifetimes_GUI_3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global CurrDirName
global DataCell
global SurfParmCell

DataCell={'SintonLiftime XY','PathName','FileName','SampleName','datenumber', 'thickness', 'user resistivity', 'calculated res', '1 sunsIvoc', 'wafertype', 'Teff at Specif_DelN'};
SurfParmCell={'SampleName','DelN,S0,S1,S2','DelN,J0s-6Mackel,J0s-Kimmerle','Ndop','TauEff at Specif_DelN','JoeAverage','DeltaN for J0s','iVoc at Specif_DelN','Seff0 at Specif_DelN','Seff2 at Specif_DelN','TauSRHparms,Tn,Tp,Et','Fit resnorm'};

% CurrDirName='C:\Users\CORP1959\Google Drive\Research\1. Silicon Surface Passivation\Experimental Data\Lifetime\Single\';
CurrDirName='';
    
% handles.figure1.Resize='off';
handles.figure1.MenuBar='figure';

axes(handles.axes1);

axis([1e14 1e17 1e-5 1e-2]); 
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set( gca, 'DataAspectRatioMode', 'auto' );
xlabel('Injection level \Deltan (cm^{-3})'); 
ylabel('Effective lifetime (s)');
title('Sinton Lifetime Plotter');

lgd=legend('show');
lgd.Interpreter= 'none';
lgd.Location= 'none';
lgd.Position=[0.7 0.7 0.2 0.1];
lgd.FontSize= 10;

hold all

set(gca,'LineStyleOrder',{'-','--',':'});
set(handles.listbox1, 'Value', []);set(handles.listbox1, 'String', {});


% UIWAIT makes PlotLifetimes_GUI_3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotLifetimes_GUI_3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global DataCell

try

    for i=1:length(handles.listbox1.Value)
       
        IndexArray=handles.listbox1.Value(i)+1;
        updatefigure(handles.axes1,IndexArray);
        
    end
    catch err
        errordlg(err.message,' Error');
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% set(handles.listbox1, 'Value', []);


% --- Executes on button press in ButtonLoadData.
function ButtonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ButtonLoadData


% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CurrDirName
global DataCell

try

DataIndex=size(DataCell);

button = questdlg('Importing Data from Sinton WCT120 files','Please select','Import File','Import Folder','Cancel','Cancel');
switch button
    case 'Import File'
        
        [FileName,PathName,FilterIndex] = uigetfile(strcat(CurrDirName,'*.xlsm'));
        DataCell=importsintondata(FileName,PathName,DataCell);   
        CurrDirName=PathName;

        IndexArray=DataIndex(1)+1;
        updatefigure(handles.axes1,IndexArray);
        
        % Add sample to the list
        inputFileNames = get(handles.listbox1,'String');
        inputFileNames{end+1} = FileName(1:end-5);
        set(handles.listbox1,'String',inputFileNames);
        
    case 'Import Folder'
        
        PathName =strcat(uigetdir(CurrDirName),'/');
        listing = dir(PathName);
        FolderSize = length(listing);

        for i=1:FolderSize
            if isempty(strfind(listing(i).name,'.xlsm'))==0 && isempty(strfind(listing(i).name,'._'))==1 && isempty(strfind(listing(i).name,'~$'))==1
                FileName=listing(i).name;
                DataCell=importsintondata(FileName,PathName,DataCell); 
                
                IndexArray=DataIndex(1)+1;
                updatefigure(handles.axes1,IndexArray);

                % Add sample to the list
                inputFileNames = get(handles.listbox1,'String');
                inputFileNames{end+1} = FileName(1:end-5);
                set(handles.listbox1,'String',inputFileNames);
                DataIndex=size(DataCell);
            end
        end
        
        CurrDirName=PathName;
        
    case 'No thank you'
end

catch err
        errordlg(err.message,' Error');
end
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataCell

try

    for i=1:length(handles.listbox1.Value)
        
        IndexArray=handles.listbox1.Value(i)+1;
        SampleName=DataCell{IndexArray,4};
        
        FigHndls=handles.axes1.Children;
        for i=1:length(FigHndls)
            if strcmp(FigHndls(i).DisplayName,SampleName)
                delete(FigHndls(i));
            end
        end
    end
    catch err
        errordlg(err.message,' Error');
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateaxeslims(handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateaxeslims(handles)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateaxeslims(handles)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateaxeslims(handles)
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataCell SurfParmCell

try
    axes(handles.axes1); cla;
    handles.listbox1.Value=[];
    handles.listbox1.String={};
    DataCell={'SintonLiftime XY','PathName','FileName','SampleName','datenumber', 'thickness', 'user resistivity', 'calculated res', '1 sunsIvoc', 'wafertype', 'Teff at Specif_DelN'};
    SurfParmCell={'SampleName','DelN,S0,S1,S2','DelN,J0s-6Mackel,J0s-Kimmerle','Ndop','TauEff at Specif_DelN','JoeAverage','DeltaN for J0s','iVoc at Specif_DelN','Seff0 at Specif_DelN','Seff2 at Specif_DelN','TauSRHparms,Tn,Tp,Et','Fit resnorm'};
    
catch err
    errordlg(err.message,' Error');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global DataCell SurfParmCell
errCase=0;

try

    switch handles.popupmenu2.Value
       
        case 2 % see data
            evalin('base','global DataCell');
            open DataCell
            
        case 3 % print to figure

            Print2Figure(handles);
            
        case 4 % extract recombination parameters
            
            if isempty(handles.listbox1.Value); errCase=1; end
            
            for i=1:length(handles.listbox1.Value)
                IndexArray=handles.listbox1.Value(i)+1;
        
                DeltaN=DataCell{IndexArray,1}(:,1);
                TauEff=DataCell{IndexArray,1}(:,2);
                W=DataCell{IndexArray,6};
                Resistivity=DataCell{IndexArray,7};
                Sitype=DataCell{IndexArray,10};Sitype=Sitype(1);
                DeltaNkey=1e15;

                SampleName=handles.listbox1.String{handles.listbox1.Value(i)};
                
                options.Interpreter = 'tex';
                DlgResults = inputdlg({'Fix value of \tau_{N0} (0 for fit)',...
                    'Fix value of \tau_{P0} (0 for fit)',...
                    'Fix E_t-E_v (0 for fit)',...
                    'Fix value of J_{0S}(0 for fit)',...
                    'Plot (1=active, 0=inactive)','Specify \Deltan'},...
                    'Processing options',[1 50],{'0','0','0.56','0','1','1e15'},options);
                
                SRHfitOptions=[eval(DlgResults{1}),eval(DlgResults{2}),eval(DlgResults{3}),eval(DlgResults{4})];
                
                SurfParmCell=ExtractRecombinationParms(DeltaN,TauEff,W,Resistivity, Sitype, eval(DlgResults{5}),eval(DlgResults{6}),SampleName,SurfParmCell,SRHfitOptions);
                
            end
            
        case 5 % see data
            evalin('base','global SurfParmCell');
            open SurfParmCell
            
        case 6% plot lifetime at specified delta n
            prompt={'Enter \Deltan (cm^{-3})'};
            name = 'Minority carrier injectrion';
            defaultans = {'1e15'};
            options.Interpreter = 'tex';
            SpecifiedMCD = inputdlg(prompt,name,[1 40],defaultans,options);
            Xaxisoption = questdlg('Specify x axis for plot', ...
                        'Independent variable', ...
                        'Date','Absolute Time','Measurement number','Measurement number');
            SpecifiedMCD=eval(SpecifiedMCD{1});
            switch Xaxisoption
                case 'Date'
                    xaxisis = 1;
                case 'Absolute Time'
                    xaxisis = 2;    
                case 'Measurement number'
                    xaxisis = 3;
                otherwise
                    return
            end
            
            plotlifetimemcd(SpecifiedMCD,xaxisis,handles.listbox1.Value);
    end
          

catch err
    switch errCase
        case 1
        errordlg('Please specify one sample to analise',' Error');
        otherwise
        errordlg(err.message,' Error');   
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Down here are all processing functions that do not have to do with GUI


function updateaxeslims(handles)
try
    y2=eval(handles.edit1.String);
    y1=eval(handles.edit2.String);
    x1=eval(handles.edit3.String);
    x2=eval(handles.edit4.String);
    axes(handles.axes1);
    axis([x1 x2 y1 y2]); 
catch err
        errordlg(err.message,' Error');
end

function Print2Figure(handles)

figure;
hdlax=axes;
y2=eval(handles.edit1.String);
y1=eval(handles.edit2.String);
x1=eval(handles.edit3.String);
x2=eval(handles.edit4.String);
axis([x1 x2 y1 y2]);

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlabel('Injection level \Deltan (cm^{-3})'); 
ylabel('Effective lifetime (s)');
lgd=legend('show');
lgd.Interpreter= 'none';
lgd.Location= 'none';
lgd.Position=[0.7 0.7 0.2 0.1];
lgd.FontSize= 10;
hold all;
copyobj(handles.axes1.Children,hdlax)

function updatefigure(axeshandle,IndexArray)
global DataCell

axes(axeshandle);

plot(DataCell{IndexArray,1}(:,1),DataCell{IndexArray,1}(:,2),...
    'DisplayName',DataCell{IndexArray,4},...
    'LineWidth',2);

function plotlifetimemcd(SpecifiedMCD,xaxisis,ix)
global DataCell

figure;

Tau_SMCD=zeros(length(ix),2);
CellTau_SMCD=cell(length(ix),3);
for i=1:length(ix)
    Tau_SMCD(i,1)=DataCell{ix(i)+1,5};
    Tau_SMCD(i,2)=interp1(DataCell{ix(i)+1,1}(:,1),DataCell{ix(i)+1,1}(:,2),SpecifiedMCD,'linear','extrap');
    CellTau_SMCD(i,:)={DataCell{ix(i)+1,4},Tau_SMCD(i,1),Tau_SMCD(i,2)};
end

global CellTau_SMCD
evalin('base','global CellTau_SMCD');
            open CellTau_SMCD

switch xaxisis
    case 1
        plot(Tau_SMCD(:,1),Tau_SMCD(:,2),'o');
        datetick('x', 'dd, mmm, yyyy, HH:MM:SS')
        
    case 2
        plot((Tau_SMCD(:,1)-min(Tau_SMCD(:,1)))*60*24,Tau_SMCD(:,2),'o');
        xlabel('Time (min)'); 
    case 3
        plot(Tau_SMCD(:,2),'o');
        xlabel('Meas Number'); 
     
end

set(gca, 'YScale', 'log');
ylabel('Effective lifetime (s)');
