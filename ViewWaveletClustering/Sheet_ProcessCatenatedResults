function varargout = Sheet_ProcessCatenatedResults_150720(varargin)
% SHEET_PROCESSCATENATEDRESULTS_150720 MATLAB code for Sheet_ProcessCatenatedResults_150720.fig
%      SHEET_PROCESSCATENATEDRESULTS_150720, by itself, creates a new SHEET_PROCESSCATENATEDRESULTS_150720 or raises the existing
%      singleton*.
%
%      H = SHEET_PROCESSCATENATEDRESULTS_150720 returns the handle to a new SHEET_PROCESSCATENATEDRESULTS_150720 or the handle to
%      the existing singleton*.
%
%      SHEET_PROCESSCATENATEDRESULTS_150720('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHEET_PROCESSCATENATEDRESULTS_150720.M with the given input arguments.
%
%      SHEET_PROCESSCATENATEDRESULTS_150720('Property','Value',...) creates a new SHEET_PROCESSCATENATEDRESULTS_150720 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sheet_ProcessCatenatedResults_150720_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sheet_ProcessCatenatedResults_150720_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sheet_ProcessCatenatedResults_150720

% Last Modified by GUIDE v2.5 20-Jul-2015 17:39:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Sheet_ProcessCatenatedResults_150720_OpeningFcn, ...
    'gui_OutputFcn',  @Sheet_ProcessCatenatedResults_150720_OutputFcn, ...
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

% --- Executes just before Sheet_ProcessCatenatedResults_150720 is made visible.
function Sheet_ProcessCatenatedResults_150720_OpeningFcn(hObject, eventdata, handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sheet_ProcessCatenatedResults_150720 (see VARARGIN)

% Choose default command line output for Sheet_ProcessCatenatedResults_150720
handles.output = hObject;

if ispc
    handles.Pref = 'C:/Users/jbilling/Documents/';
else
    handles.Pref = '/data1/shares/brainiac-shared/';
end

% Update handles structure
guidata(hObject, handles);

wlev = 6;
handles.nodes = zeros(wlev+1,1*2^wlev);
for wl = 1:wlev
    handles.nodes(wl+1,1:2^wl) = length(find(handles.nodes))+(1:2^wl);  %Move to previous
end

handles.WvCldata = 0;
handles.ICAdata = 0;
handles.FCdata = 0;
handles.UseIntelDissimClust_radio = 1;
handles.pwdStart = mfilename('fullpath');
[fpath,~,~] = fileparts(handles.pwdStart);
handles.pwdStart = fpath;
% handles.AuxData = load(fullfile(fpath,filesep,'CombinedDataset.mat'));
% handles.AuxData = handles.AuxData.AuxData_EuclideanAverage;

guidata(hObject,handles)

% UIWAIT makes Sheet_ProcessCatenatedResults_150720 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Sheet_ProcessCatenatedResults_150720_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

%--Helper function for plotters
function getdatafun(handles,GUIname,varargin)

WavDepth = handles.d{handles.repCount};
WavPosition = handles.p{handles.repCount};

v = handles.v{handles.repCount};

Ind = handles.nodes(WavDepth+1,WavPosition+1)+1;

Z = handles.ResultsSaves(v).Z{Ind};
score = handles.ResultsSaves(v).score(:);

cutoff = 1-handles.c{handles.repCount};
g = handles.g{handles.repCount};
Inconsist = inconsistent(Z,g);

OrdInc = sort(Inconsist(:,4),'descend');

minOIc = min(OrdInc(OrdInc~=0))+2*eps;
maxOIc = max(OrdInc(OrdInc~=0))-3*eps;

if cutoff>=0 && cutoff<=1
    OIc = (maxOIc-minOIc)*cutoff+minOIc;
else
    error('Inconsistency threshold must be a percentage between 0 and 1')
end

[T,conn] = cluster(Z,'cutoff',OIc,'Depth',g);
handles.conn = conn;
ndxClust = uint16(T);
ClusterNumbers = double(max(ndxClust));
if ClusterNumbers == 1
    error('Threshold would return 1 cluster')
end

handles.ClusterNumbers{handles.repCount} = ClusterNumbers;

maxclst = size(Z,1);

DenLocs = handles.ResultsSaves(v).DendrogramLocs{Ind};

fpath = handles.pwdStart;

dataPath = get(handles.BaseDataFolder_edit,'string');

handles.DendroFilename{handles.repCount} = ...
    fullfile(dataPath,filesep,handles.Datasets(v,:),'Dendrograms',filesep,...
    ['Dendro_D' num2str(WavDepth) 'P' num2str(WavPosition) '.fig']);

handles.Zloc{handles.repCount} = Z(size(Z,1)-ClusterNumbers+2,3);
handles.Zlen{handles.repCount} = size(Z,1);

%Clusters are defined by dendrogram positions
%ndxClust is reorganized to dendrogram positions

ndxClust1 = ndxClust;
handles.ClusterXTickLocs = zeros(2,ClusterNumbers);
for ii = 1:ClusterNumbers
    handles.ClusterXTickLocs(1,ii) = round(mean(DenLocs(ndxClust==ii)));
    ndxClust1(ndxClust==ii) = handles.ClusterXTickLocs(1,ii);
end
[handles.ClusterXTickLocs(1,:)] = sort(handles.ClusterXTickLocs(1,:)); %Build 1:ClusterNumbers list from ClusterNumbers distributed across 1:NVoxels numbers
for ii = 1:ClusterNumbers
    ndxClust(ndxClust1==handles.ClusterXTickLocs(1,ii)) = ii;
    handles.ClusterXTickLocs(2,ii) = ii;
end

handles.ndxClust{handles.repCount} = ndxClust;

data = zeros(32,38,28,1);
data(score==1) = ndxClust(:);

Exclude = str2num(get(handles.Exclude_Edit,'String'));
Include = setdiff(1:ClusterNumbers,Exclude);
if ~isempty(Exclude)
    fun = @(x)find(data==x);
    vec = [];
    for i = 1:length(Exclude)
        vec = [vec fun(Exclude(i))'];
    end
    RemainSet = get(handles.SetRemainClust_edit,'string');
    if ~isempty(RemainSet)
        data(vec) = str2double(RemainSet);
    else
        data(vec) = round(mean(Include)+3*(max(std(Include),mean(Include))));
    end
end

% handles.capt{handles.repCount} = char({['D' num2str(WavDepth) 'P' num2str(WavPosition)];['IF = ' num2str(OIc)];['nC=' num2str(ClusterNumbers)]});
handles.capt{handles.repCount} = char({[''];['IT = ' num2str(OIc)];['|C| = ' num2str(ClusterNumbers)]});
%         data = rand([4,4,4]);
%         data(data<.8) = 0;
%         data = round(data);
sizedata = size(data);
mask = zeros(sizedata);
mask(data>0) = 1;
mask = repmat(mask,[1,1,1,3]);
data_nii = make_nii(data,[5 5 5],[16.4 23.4 11.0]);
save_nii(data_nii, fullfile(handles.pwdStart, ['ProcessedCat' num2str(handles.repCount)]));
data_nii = make_nii(mask,[5 5 5],[16.4 23.4 11.0]);
save_nii(data_nii, fullfile(handles.pwdStart,['ProcessedCatRGBmask' num2str(handles.repCount)]));

zlen = handles.Zlen{handles.repCount};
hlen = 200;
%         vrng = [0.70,0.99];
vrng = [0.60,1];
srng = [0.85,1];

hmap = linspace(0,hlen/240,hlen);
if ClusterNumbers>hlen
    hmap = repmat(hmap,[ceil(zlen/hlen),1]);
    hmap = hmap(:);
    vmap = linspace(vrng(1),vrng(2),ceil(zlen/hlen/2))';
    vmap = [ones(size(vmap));vmap];
    vmap = repmat(vmap,[hlen,1]);
    smap = linspace(srng(1),srng(2),ceil(zlen/hlen/2))';
    smap = [smap;ones(size(smap))];
    smap = repmat(smap,[hlen,1]);
    newcolormap = [[0,0,0];[hmap(1:zlen),vmap(1:zlen),smap(1:zlen)]];
else
    hmap = hmap(:);
    newcolormap = [[0,0,0];[hmap,ones(size(hmap)),ones(size(hmap))]];
end
newcolormap = hsv2rgb(newcolormap);
if get(handles.FixedColormap_checkbox,'value')
    if ClusterNumbers>hlen
        newcolormap = newcolormap([1,handles.ClusterXTickLocs(1,:)],:);
    else
        lc = handles.ClusterXTickLocs(1,:)+1;
        newcolormap = newcolormap([1,ceil((1-(zlen-lc)/zlen)*hlen)],:);
    end
else
    newcolormap = newcolormap([1,round(linspace(2,length(newcolormap),ClusterNumbers))],:);
end
handles.colormap{handles.repCount} = newcolormap;

% I want to build a colorbar that matches the plot

handles.PlotColorbar{handles.repCount} = zeros(1,length(ndxClust),3);
for pb = 1:length(ndxClust);
    handles.PlotColorbar{handles.repCount}(1,DenLocs(pb),:) = newcolormap(ndxClust(pb)+1,:);
end

% db = 0;
% if db == 1
%     figure(40); clf
%     imagesc(handles.PlotColorbar{handles.repCount})
% end

if ~isempty(Exclude)
    data(vec) = ClusterNumbers; %First resets excluded data values to last RGB value
end

dataRGB = reshape(newcolormap(data(:)+1,:),[32,38,28,3]);

if get(handles.StopForData_checkbox,'value')
    dataA = data;
    save('Make3dPlotsData','dataA')
    keyboard;
end

data_nii = make_nii(dataRGB,[5 5 5],[16.4 23.4 11.0]);
save_nii(data_nii, fullfile(handles.pwdStart, ['ProcessedCatRGB' num2str(handles.repCount)]));

if ispc
    handles.name{handles.repCount,1} = [handles.pwdStart,'\','ProcessedCat' num2str(handles.repCount) '.img'];
    handles.name{handles.repCount,2} = [handles.pwdStart,'\','ProcessedCatRGB' num2str(handles.repCount) '.img'];
else
    handles.name{handles.repCount,1} = [handles.pwdStart,'/','ProcessedCat' num2str(handles.repCount) '.img'];
    handles.name{handles.repCount,2} = [handles.pwdStart,'/','ProcessedCatRGB' num2str(handles.repCount) '.img'];
end

handles.WvCldata = 1;

guidata(GUIname,handles)


function WaveletDepth_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to WaveletDepth_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaveletDepth_Edit as text
%        str2double(get(hObject,'String')) returns contents of WaveletDepth_Edit as a double

handles.WvCldata = 0;
handles.FCdata = 0;
handles.ICAdata = 0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function WaveletDepth_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveletDepth_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WaveletPosition_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to WaveletPosition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaveletPosition_Edit as text
%        str2double(get(hObject,'String')) returns contents of WaveletPosition_Edit as a double

handles.WvCldata = 0;
handles.FCdata = 0;
handles.ICAdata = 0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function WaveletPosition_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveletPosition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ClusterNumbers_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterNumbers_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ClusterNumbers_Edit as text
%        str2double(get(hObject,'String')) returns contents of ClusterNumbers_Edit as a double

set(handles.UseClustNums_radio,'value',1)

handles.WvCldata = 0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ClusterNumbers_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClusterNumbers_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [name, capt, mD] = getanatomical(handles,hObject)

v = handles.v{handles.repCount};

dataPath = get(handles.BaseDataFolder_edit,'string');
AnatPath = fullfile(dataPath,filesep,handles.Datasets(v,:), filesep,'Anatomical.nii');
Anat = load_untouch_nii(AnatPath);

capt = char({handles.Datasets(v,1:14),handles.Datasets(v,16:end)});

data = Anat.img;
sizedata = size(Anat.img);
mask = zeros(sizedata);
mask(data>12) = 1;
fmask = repmat(mask,[1,1,1,3]);
Anat.img = mask;
ext = '.nii';
save_untouch_nii(Anat, fullfile(handles.pwdStart, ['RGBmask_Anatomical' handles.Datasets(v,:)  num2str(handles.repCount) ext]));

Anat.img = round(double(data));
Anat.img(data<12) = 0;
mD = max(Anat.img(:));

newcolormap = repmat(linspace(0,mD,mD+1)',3);

% if get(handles.FixedColormap_checkbox,'value')
%     newcolormap = newcolormap([1,handles.ClusterXTickLocs(1,:)+1],:);
% else
%     newcolormap = newcolormap([1,round(linspace(2,length(newcolormap),ClusterNumbers))],:);
% end
% handles.colormap{handles.repCount} = newcolormap;

Anat.img = reshape(newcolormap(Anat.img(:)+1,:),[sizedata(1:3),3]);

Anat.hdr.dime.dim([1,5]) = [4,3];
Anat.hdr.dime.pixdim(5) = 1;

save_untouch_nii(Anat, fullfile(handles.pwdStart, ['RGB_Anatomical' handles.Datasets(v,:) num2str(handles.repCount) ext]));

name{1} = AnatPath;
name{2} = fullfile(handles.pwdStart, ['RGB_Anatomical' handles.Datasets(v,:)  num2str(handles.repCount) ext]);
% name{3} = ['RGBmask_' Locs{v} num2str(handles.repCount) ext];


% --- Executes on button press in PlotClusters_Button.
function PlotClusters_Button_Callback(hObject, eventdata, handles)
% hObject    handle to PlotClusters_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'string','Busy','BackgroundColor',[1 0 0])
pause(.5)

handles = guidata(hObject);
handles.name = [];

if ~isfield(handles,'ResultsSaves') || isempty(handles.ResultsSaves)
    LoadResults(hObject,handles)
    handles = guidata(hObject);
end

gVals = eval(get(handles.Inconsistencyg_Edit,'string'));
InconVals = eval(get(handles.InconsistencyCutoff_Edit,'string'));

lRS = length(handles.ResultsSaves);

handles.WDepths = eval(get(handles.WaveletDepths_edit,'string'));
handles.WPositions = eval(get(handles.WaveletPositions_edit,'string'));
handles.InconCutoff = eval(get(handles.InconsistencyCutoff_Edit,'string'));
handles.Gvalues = eval(get(handles.Inconsistencyg_Edit,'string'));

handles.repCount = 0;

lWD = length(handles.WDepths);
lWP = length(handles.WPositions);
assert(lWD==lWP);

for rs = 1:lRS
    
    newV = 1;
    for wi = 1:lWD
        d = handles.WDepths(wi);
        p = handles.WPositions(wi);
        
        g = handles.Gvalues;
        c = handles.InconCutoff;
        
        assert(isequal(length(g),length(c)),'Values for g and the inconsistency threshold must be of equal length')
        
        for gcInd = 1:length(g)
                handles.repCount = handles.repCount+1;
                rep = handles.repCount;
                
                if newV == 1 && get(handles.IncludeAnatomy_checkbox,'value');
                    handles.v{rep} = rs;
                    handles.d{rep} = d;
                    handles.p{rep} = p;
                    handles.g{rep} = g(gcInd);
                    handles.c{rep} = c(gcInd);
                    [name, handles.capt{handles.repCount}, handles.ClusterNumbers{handles.repCount}] = getanatomical(handles,hObject);
                    handles.name{handles.repCount,1} = name{1};
                    handles.name{handles.repCount,2} = name{2};
                    newV = 0;
                    handles.repCount = handles.repCount+1;
                    rep = handles.repCount;
                end
                
                handles.v{rep} = rs;
                handles.d{rep} = d;
                handles.p{rep} = p;
                handles.g{rep} = g(gcInd);
                handles.c{rep} = c(gcInd);
                
                getdatafun(handles,gcbo,'hObject',hObject);
                
                handles = guidata(hObject);
           
        end
    end
end

sHn = size(handles.name);

name1 = [];

for ii = 1:sHn(1)
    name1{ii} = handles.name{ii,1};
end
name1 = char(name1);

capt = handles.capt;

ClusterNumbers = handles.ClusterNumbers;

spm_check_registration(name1,capt,[],'fontsize',24,'color','w')
spm_orthviews('Interp',0)
spm_orthviews('Reposition',[6.57, -15.33, 29.37])
% spm_orthviews('Reposition',[3.6,-62.9,32])
xhairs = {'off','on'};
xhairs = xhairs{get(handles.Crosshairs_checkbox,'value')+1};
spm_orthviews('Xhairs',xhairs)

for ii = 1:sHn(1)
    spm_orthviews('Window',ii,[1, ClusterNumbers{ii}+1])
    
    colormap(ones(ClusterNumbers{ii}+1,3))
    
    ProcRGBNames = [handles.name{ii,2},',1';handles.name{ii,2},',2';handles.name{ii,2},',3'];
    spm_orthviews('rgb', 'Draw_NoUI', ii, ProcRGBNames )
    %                     spm_orthviews('Window',Wind,[ClusterNumbers-1, ClusterNumbers+1])
    spm_orthviews('Window',ii,[ClusterNumbers{ii}, ClusterNumbers{ii}+1])
end

handles.svDir = get(handles.SaveSheet_edit,'string');
if exist(handles.svDir,'dir')~=7
    handles.svDir = fullfile(handles.pwdStart, 'ImageSaves');
end
svDir = handles.svDir;
set(gcf,'InvertHardCopy','off');

if get(handles.SaveImages_checkbox,'value')
    XhairPositions = spm_orthviews('Pos');
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'color','black');
    hAllAxes = findobj(gcf,'type','axes');
    slcs = {'Sagittal','Coronal','Horozontal'}; % Probably should reorder this to alleviate the following comment...
    mCount = 0; % Actually a bit confusing working to order images to axes handles.
    for ms = (sHn(1)*3):-3:1
        mCount = mCount + 1;
        slCount = 0;
        for ha = 2:-1:0
            slCount = slCount+1;
            F = getframe(hAllAxes(ms-ha));
            Image = frame2im(F);
            set(gca,'color','black');
            set(gcf,'InvertHardCopy','off');
            if isempty(handles.Zlen{mCount})
                imwrite(Image,fullfile(svDir,filesep,[handles.DatasetsIDs{handles.v{mCount}} '_Anatomical_Slices_' slcs{slCount} '_At_' num2str(XhairPositions(slCount)) '.png']),'png')
            else
                imwrite(Image,fullfile(svDir,filesep,[handles.DatasetsIDs{handles.v{mCount}} '_D' num2str(handles.d{mCount}) 'P' num2str(handles.p{mCount}) '_EuAvg_Leaves' num2str(handles.ClusterNumbers{mCount}) '_Slices_' slcs{slCount} '_At_' num2str(XhairPositions(slCount)) '.png']),'png')
            end
        end
    end
end

if get(handles.IncludeDendro_checkbox,'value') && rep == 1;
    figure(40)
    imagesc(handles.PlotColorbar{rep})
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gcf,'position',[231         77        1448          30])
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'color','black');
    set(gca,'color','black');
    if get(handles.SaveImages_checkbox,'value')
        set(gcf,'InvertHardCopy','off');
        saveas(gcf,[svDir filesep 'D' num2str(handles.d{rep}) 'string' 'P' num2str(handles.p{rep})  '_ColorBar_' num2str(handles.ClusterNumbers{rep}) 'Clusts.png'],'png')
        %     saveas(gcf,[svDir 'D' get(handles.WaveletDepth_Edit,'string') 'P' get(handles.WaveletPosition_Edit,'string') '_ColorBar_' num2str(handles.ClusterNumbers) 'Clusts.fig'],'fig')
    end
    if isfield(handles,'CurrentDendroFig')
        if ~strcmp(handles.CurrentDendroFig,handles.DendroFilename)
            try
                close(handles.DendrogramFig)
            end
        end
    end
    handles.DendrogramFig = openfig(handles.DendroFilename{rep},'reuse');
    figure(handles.DendrogramFig);
    handles.CurrentDendroFig = handles.DendroFilename;
    set(gca,'xtick',handles.ClusterXTickLocs(1,:))
    set(gca,'xticklabel',handles.ClusterXTickLocs(2,:))
    hold on
    if get(handles.BuildLines_checkbox,'value')
        DendLineCol = get(handles.DendroLineColor_edit,'string');
        if ~handles.UseIntelDissimClust_radio
            if isfield(handles,'ClustLine') && get(handles.DeleteLines_checkbox,'value')
                try
                    delete(handles.ClustLine)
                end
                handles.ClustLine = plot([1,handles.Zlen],[handles.Zloc,handles.Zloc],DendLineCol,'Linewidth',0.5);
            elseif isfield(handles,'ClustLine') && ~get(handles.DeleteLines_checkbox,'value')
                handles.ClustLine = [handles.ClustLine plot([1,handles.Zlen],[handles.Zloc,handles.Zloc],DendLineCol,'Linewidth',0.5)];
            elseif ~isfield(handles,'ClustLine')
                handles.ClustLine = plot([1,handles.Zlen],[handles.Zloc,handles.Zloc],DendLineCol,'Linewidth',0.5);
            end
        else
            
            if get(handles.DeleteLines_checkbox,'value')
                try
                    delete(handles.ClustLine)
                end
            end
            
            h = findobj(gca,'type','line');
            xdata = cell2mat(get(h,'Xdata'));
            ydata = cell2mat(get(h,'Ydata'));
            v = handles.v{handles.repCount};
            Ind = depo2ind(2,[d,p])+1;
            Z = handles.ResultsSaves(v).Z{Ind};
            conn = handles.conn;
            links = Z(conn==0,3);
            linklocs = ismember(ydata(:,2),links);
            hold on
            if isfield(handles,'ClustLine') && get(handles.DeleteLines_checkbox,'value')
                try
                    delete(handles.ClustLine)
                end
                handles.ClustLine = [];
                for lk = find(linklocs)'
                    %                 [X,Y] = calculateEllipse(mean(xdata(lk,2:3)),ydata(lk,2),diff(xdata(lk,2:3))/2,diff(ydata(lk,1:2))/2,0);
                    %                 handles.ClustLine = [handles.ClustLine, plot(X,Y,DendLineCol,'linewidth',1)];
                    handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),DendLineCol,'linewidth',0.5)];
                    %  handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),'w','linewidth',0.5)];
                end
            elseif isfield(handles,'ClustLine') && ~get(handles.DeleteLines_checkbox,'value')
                for lk = find(linklocs)'
                    %                 [X,Y] = calculateEllipse(mean(xdata(lk,2:3)),ydata(lk,2),diff(xdata(lk,2:3))/2,diff(ydata(lk,1:2))/2,0);
                    %                 handles.ClustLine  = [handles.ClustLine, plot(X,Y,DendLineCol,'linewidth',1)];
                    handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),DendLineCol,'linewidth',0.5)];
                    %  handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),'w','linewidth',0.5)];
                end
            elseif ~isfield(handles,'ClustLine')
                handles.ClustLine = [];
                for lk = find(linklocs)'
                    %                 [X,Y] = calculateEllipse(mean(xdata(lk,2:3)),ydata(lk,2),diff(xdata(lk,2:3))/2,diff(ydata(lk,1:2))/2,0);
                    %                 handles.ClustLine = [handles.ClustLine, plot(X,Y,DendLineCol,'linewidth',1)];
                    handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),DendLineCol,'linewidth',0.5)];
                    %  handles.ClustLine = [handles.ClustLine, plot(xdata(lk,:),ydata(lk,:),'w','linewidth',0.5)];
                end
            end
            hold off
        end
    end
    
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'color','black');
    set(gca,'color','black');
    set(gca,'ycolor','w')
    set(gca,'fontsize',20)
    set(gcf,'position',[229         468        1447         283])
    
    if get(handles.SaveImages_checkbox,'value')
        set(gcf,'InvertHardCopy','off');
        %         saveas(gcf,[svDir 'D' get(handles.WaveletDepth_Edit,'string') 'P' get(handles.WaveletPosition_Edit,'string') '_Dendro_' num2str(handles.ClusterNumbers) 'Clusts.png'],'png')
        saveas(gcf,[svDir filesep 'D' num2str(handles.d{rep}) 'P' num2str(handles.p{rep}) '_Dendro_' num2str(handles.ClusterNumbers{rep}) 'Clusts.png'],'png')
    end
end

set(hObject,'string','Plot','BackgroundColor',[.941 .941 .941])
pause(.5)

guidata(hObject,handles)


function Exclude_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Exclude_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Exclude_Edit as text
%        str2double(get(hObject,'String')) returns contents of Exclude_Edit as a double

handles.WvCldata = 0;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Exclude_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Exclude_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Crosshairs_checkbox.
function Crosshairs_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Crosshairs_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Crosshairs_checkbox


% --- Executes on button press in ResetSaves_button.
function ResetSaves_button_Callback(hObject, eventdata, handles)
% hObject    handle to ResetSaves_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FCdata = 0;
handles.ICAdata = 0;
handles.WvCldata = 0;

try
    close(handles.DendrogramFig)
end

guidata(hObject,handles)

% --- Executes on button press in StopForData_checkbox.
function StopForData_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to StopForData_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StopForData_checkbox

handles.stop = get(handles.StopForData_checkbox,'value');
if handles.stop
    handles.FCdata = 0;
    handles.ICAdata = 0;
    handles.WvCldata = 0;
end
guidata(hObject,handles)

% --- Executes on button press in IncludeDendro_checkbox.
function IncludeDendro_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeDendro_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IncludeDendro_checkbox


% --- Executes on button press in DeleteLines_checkbox.
function DeleteLines_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteLines_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DeleteLines_checkbox


function SetRemainClust_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SetRemainClust_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SetRemainClust_edit as text
%        str2double(get(hObject,'String')) returns contents of SetRemainClust_edit as a double


% --- Executes during object creation, after setting all properties.
function SetRemainClust_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SetRemainClust_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SelectDendrogramNodeforData_PushButton.
function SelectDendrogramNodeforData_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectDendrogramNodeforData_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = openfig(handles.DendroFilename,'reuse');
datacursormode(fig);
dcm_obj = datacursormode(fig);
disp('Hit enter when you have selected the link of interest')
pause
c_info = getCursorInfo(dcm_obj);

WavDepth = str2double(get(handles.WaveletDepth_Edit,'String'));
WavPosition = str2double(get(handles.WaveletPosition_Edit,'String'));
Z = handles.Z;
eval(['score = handles.AuxData.AuxData_' handles.AnalysisType '.score;'])
finds = find(score);

row = find(Z(:,3)>c_info.Position(2),1,'first'); %This is not the actual row selected so the following lines are used for correction
augment = -2:2;
[temp, shift] = min(abs(Z(row+[augment],3)-c_info.Position(2)));
row = row+augment(shift);

tree_indices = [];
tree_indices = getalltree(row,Z,size(Z,1),tree_indices);

load('CatenatedDataSave.mat', ['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)])
data = eval(['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)]); clear(['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)])
data = reshape(data,[],size(data,4));
AveragedSubtree = data(finds(tree_indices),:);
AveragedSubtree = mean(AveragedSubtree,1);
save('AveragedSubtree','AveragedSubtree')

function tree_indices = getalltree(row,Z,szZ,tree_indices)

for i = 1:2
    if Z(row,i) <= (szZ+1)
        tree_indices = [tree_indices, Z(row,i)];
    else
        RowNodeCreated = Z(row,i)-(szZ+1);
        tree_indices = getalltree(RowNodeCreated,Z,szZ,tree_indices);
    end
end

function InconsistencyCutoff_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to InconsistencyCutoff_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InconsistencyCutoff_Edit as text
%        str2double(get(hObject,'String')) returns contents of InconsistencyCutoff_Edit as a double

handles.WvCldata = 0;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function InconsistencyCutoff_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InconsistencyCutoff_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AllFromData_PushButton.
function AllFromData_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to AllFromData_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

WavDepth = str2double(get(handles.WaveletDepth_Edit,'String'));
WavPosition = str2double(get(handles.WaveletPosition_Edit,'String'));
Z = handles.Z;
ndxClust = handles.ndxClust;
eval(['score = handles.AuxData.AuxData_' handles.AnalysisType '.score;'])
finds = find(score);
maxClust = max(ndxClust);

% DataCentroid = zeros(maxClust,3);
%
% for i = 1:maxClust
%     clsts = find(ndxClust==i);
%     [yy,xx,zz] = ind2sub([32,38,28],finds(clsts));
% %     DataCentroid(i,:) = [round(mean(yy)),round(mean(xx)),round(mean(zz))];
%     DataCentroid(i,:) = [mean(yy),mean(xx),mean(zz)];
% end

if get(handles.IndividualData_radio,'value')
    vol = get(handles.WhichIndividual_edit,'string');
    fid = fopen('VolunteerVSaveNames.txt');
    scan = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    VolFName = scan{1}{str2num(vol)};
    if get(handles.GetDataBroadband_checkbox,'value')
        load(VolFName,'VSave1')
        data = zeros(32,38,28,length(VSave1(1).wptCoeff));
        for i = 1:length(VSave1)
            data(VSave1(i).y,VSave1(i).x,VSave1(i).z,:) = VSave1(i).wptCoeff;
        end
    else
        load(VolFName,['VSave' num2str(WavDepth+1)])
        VSave1 = eval(['VSave' num2str(WavDepth+1)]);
        clear(['VSave' num2str(WavDepth+1)]);
        data = zeros(32,38,28,length(VSave1(1).wptCoeff));
        for i = 1:length(VSave1)
            data(VSave1(i).y,VSave1(i).x,VSave1(i).z,:) = VSave1(i).wptCoeff(:,WavPosition+1);
        end
    end
else
    if get(handles.GetDataBroadband_checkbox,'value')
        load('CatenatedDataSave.mat', ['catVsaves_Dep' '1' 'Pos' '1'])
        data = eval(['catVsaves_Dep' '1' 'Pos' '1']);
        clear(['catVsaves_Dep' '1' 'Pos' '1'])
    else
        load('CatenatedDataSave.mat', ['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)])
        data = eval(['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)]); clear(['catVsaves_Dep' num2str(WavDepth+1) 'Pos' num2str(WavPosition+1)])
    end
end

data = reshape(data,[],size(data,4));

rows = max(ndxClust);
AveragedSubtree = cell(double(rows),1);
DataCentroid = zeros(maxClust,3);
% cr0 = [32,38,28];
% cr1 = [82,-117,-55];
% cr2 = [-78,73,85];
for r = 1:rows
    voxels = boolean(ndxClust==r);
    temp = data(finds(voxels),:);
    AveragedSubtree{r} = mean(temp,1);
    % %     val = 0;
    % %     tag = 1;
    % %     for i = 1:size(temp,1)
    % %         target = AveragedSubtree{r}*temp(i,:)';
    % %         if val < target
    % %             val = target;
    % %             tag = i;
    % %         end
    % %     end
    % %     locs = find(voxels);
    % %     gg = locs(tag);
    % %     loc = [VSave1(gg).x,VSave1(gg).y,VSave1(gg).z];
    % % %     DataCentroid(r,:) = [loc(1)/cr0(1)*(cr2(1)-cr1(1))+cr1(1), ...
    % % %         loc(2)/cr0(2)*(cr2(2)-cr1(2))+cr1(2) ...
    % % %         loc(3)/cr0(3)*(cr2(3)-cr1(3))+cr1(3)];
    % %     DataCentroid(r,:) = loc;
end

save('AveragedSubtree','AveragedSubtree','DataCentroid')


% --- Executes on button press in FixedColormap_checkbox.
function FixedColormap_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to FixedColormap_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixedColormap_checkbox
handles.WvCldata = 0;
guidata(hObject,handles)

% --- Executes on button press in SaveImages_checkbox.
function SaveImages_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to SaveImages_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SaveImages_checkbox


% --- Executes on button press in BuildLines_checkbox.
function BuildLines_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to BuildLines_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BuildLines_checkbox


function DendroLineColor_edit_Callback(hObject, eventdata, handles)
% hObject    handle to DendroLineColor_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DendroLineColor_edit as text
%        str2double(get(hObject,'String')) returns contents of DendroLineColor_edit as a double


% --- Executes during object creation, after setting all properties.
function DendroLineColor_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DendroLineColor_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WaveletDepths_edit_Callback(hObject, eventdata, handles)
% hObject    handle to WaveletDepths_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaveletDepths_edit as text
%        str2double(get(hObject,'String')) returns contents of WaveletDepths_edit as a double


% --- Executes during object creation, after setting all properties.
function WaveletDepths_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveletDepths_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WaveletPositions_edit_Callback(hObject, eventdata, handles)
% hObject    handle to WaveletPositions_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WaveletPositions_edit as text
%        str2double(get(hObject,'String')) returns contents of WaveletPositions_edit as a double


% --- Executes during object creation, after setting all properties.
function WaveletPositions_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveletPositions_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IncludeAnatomy_checkbox.
function IncludeAnatomy_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeAnatomy_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IncludeAnatomy_checkbox


function SaveSheet_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSheet_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaveSheet_edit as text
%        str2double(get(hObject,'String')) returns contents of SaveSheet_edit as a double

[pname] = fileparts(get(hObject,'string'));
handles.SaveSheetPname = pname;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function SaveSheet_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveSheet_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SaveSheetBrowse_pushbutton.
function SaveSheetBrowse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSheetBrowse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pname] = uigetdir;
handles.SaveSheetPname = pname;
set(handles.SaveSheet_edit,'string',fullfile(pname))
guidata(hObject,handles)

function DescriptiveTitle_edit_Callback(hObject, eventdata, handles)
% hObject    handle to DescriptiveTitle_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DescriptiveTitle_edit as text
%        str2double(get(hObject,'String')) returns contents of DescriptiveTitle_edit as a double


% --- Executes during object creation, after setting all properties.
function DescriptiveTitle_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DescriptiveTitle_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PlotCrossWPT_pushbutton.
function PlotCrossWPT_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCrossWPT_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Inconsistencyg_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Inconsistencyg_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inconsistencyg_Edit as text
%        str2double(get(hObject,'String')) returns contents of Inconsistencyg_Edit as a double

handles.WvCldata = 0;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Inconsistencyg_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inconsistencyg_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HowManyVolunteers_edit_Callback(hObject, eventdata, handles)
% hObject    handle to HowManyVolunteers_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HowManyVolunteers_edit as text
%        str2double(get(hObject,'String')) returns contents of HowManyVolunteers_edit as a double
if isfield(handles,'ResultsSaves')
    handles = rmfield(handles,'ResultsSaves');
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function HowManyVolunteers_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HowManyVolunteers_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function WhichGroup_edit_Callback(hObject, eventdata, handles)
% hObject    handle to WhichGroup_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhichGroup_edit as text
%        str2double(get(hObject,'String')) returns contents of WhichGroup_edit as a double

if isfield(handles,'ResultsSaves')
    handles = rmfield(handles,'ResultsSaves');
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function WhichGroup_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhichGroup_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotEntropies_pushbutton.
function PlotEntropies_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotEntropies_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'ResultsSaves') || isempty(handles.ResultsSaves)
    LoadResults(hObject,handles)
    handles = guidata(hObject);
end

lRS = length(handles.ResultsSaves);

for rs = 1:lRS
    
    depths = 7;
    positions = (2^(depths-1));
    EntMat = zeros(depths,positions);
    pMat = zeros(depths,positions,100);
    vertlines = zeros(sum(2.^(1:depths-1)-1),4);
    Ent = handles.ResultsSaves(rs).Ent2;
    
    for i = 1:depths
        
        count = 0;
        jvec = 0:positions/(2^(i-1)):positions;
        
        for j = 1:(2^(i-1))
            
            EntMat(i,jvec(j)+1:jvec(j+1)) = Ent(depo2ind(2,[i-1,j-1])+1); % See Notes for Entropy type
            
            if i==1; begcap = 0; endcap = 0;
            else
                endcap = 2^(i-1)-1;
                vertlines(begcap+1:begcap+endcap,3:4) = repmat([i-0.5,i+0.5],endcap,1);
            end
            if count > 0
                vertlines(begcap+count,1:2) = [jvec(j)+0.5,jvec(j)+0.5];
            end
            
            count = count+1;
        end
        begcap = begcap+endcap;
    end
    
    h = figure(9);
    % colormap(h,[cool(256);hot(256)])
    spR = [1,2,2,2,3,3,4,4,4,4,4,4,5,5,5];
    spC = [1,1,2,2,2,2,2,2,3,3,3,3,3,3,3];
    
    minEM = min(EntMat(:));
    maxEM = max(EntMat(:));
    
    EMpos = EntMat;
    EMpos(EMpos<0) = nan;
    EMposA = logspace(log10(min(EMpos(:))),log10(maxEM),112);
    cmpos = flipud(autumn(112));
    EMneg = EntMat;
    EMneg(EMneg>0) = nan;
    EMnegA = logspace(log10(-max(EMneg(:))),log10(-minEM),112);
    cmneg = winter(112);
    
    colormap([cmneg;repmat([0,0,0],32,1);cmpos])
    
    EMx = zeros(length(EntMat(:)),1);
    EMy = zeros(length(EntMat(:)),1);
    EMc = zeros(7,64,3);
    
    for i = 1:length(EntMat(:))
        [EMx(i),EMy(i)] = ind2sub(size(EntMat),i);
        if EntMat(i)>0
            targ = abs(EMposA-EntMat(i))==min(abs(EMposA-EntMat(i)));
            EMc(EMx(i),EMy(i),:) = cmpos(targ,:);
            %         EMc(EMx(i),EMy(i),:) = cmpos(ceil((EntMat(i)/maxEM)*128),:);
        elseif EntMat(i)<0
            targ = abs(EMnegA+EntMat(i))==min(abs(EMnegA+EntMat(i)));
            EMc(EMx(i),EMy(i),:) = cmneg(fliplr(targ),:);
            %         EMc(EMx(i),EMy(i),:) = cmneg(ceil((1-(EntMat(i)/minEM))*128+eps),:);
        elseif EntMat(i)==0
            EMc(EMx(i),EMy(i),:) = [0,0,0];
        end
    end
    
    subplot(spR(lRS),spC(lRS),rs)
    
    image(EMc)
    set(gca,'yticklabel',str2num(get(gca,'yticklabel'))-1)
    freqs = [1/(0.645*900)*2,1/0.645/2];
    freqs = linspace(freqs(1),freqs(2),64);
    xt = get(gca,'xtick');
    set(gca,'xticklabel',sprintf('%3.2g | ',freqs(xt)))
    cbh = colorbar;
    set(cbh,'ytickmode','manual');
    set(cbh,'yticklabelmode','manual');
    set(gca,'fontsize',13)
    % rlgsp = round(logspace(log10(1),log10(112),3));
    rlgsp = round(linspace(1,112,5));
    ytkA = [113-fliplr(rlgsp),rlgsp+112+32];
    set(cbh,'ytick',ytkA)
    minspan = linspace(minEM,0,3);
    maxspan = linspace(maxEM,0,3);
    ylbls = cell(2*length(rlgsp),1);
    for rp = 1:length(rlgsp)*2;
        if rp<=length(rlgsp)
            ylbls{rp} = sprintf('%4.2e',-1*EMnegA(ytkA(length(rlgsp)-rp+1)));
        else
            ylbls{rp} = sprintf('%4.2e',EMposA(ytkA(rp)-112-32));
        end
    end
    % set(cbh,'yticklabel',{sprintf('%4.2e',minspan(1)),sprintf('%4.2e',minspan(2)),sprintf('%4.2e',0),sprintf('%4.2e',maxspan(2)),sprintf('%4.2e',maxspan(1))})
    set(cbh,'yticklabel',ylbls)
    ylabel('Wavelet packet depth')
    xlabel('Approximate Fourier frequency (Hz)')
    title('')
    set(gca,'fontsize',13)
    
    for i = 1:size(vertlines,1); line(vertlines(i,1:2),vertlines(i,3:4),'color','k','linestyle','-'); end
    
    EntMatMask = ones(depths,positions);
    count = 0;
    for i = 2:depths
        
        jvec = 0:positions/(2^(i-1)):positions;
        
        for j = 1:2:(length(jvec)-1)
            % Zero packets based upon cost function
            if ~(sum([EntMat(i,jvec(j)+1),EntMat(i,jvec(j+1)+1)])<(EntMat(i-1,jvec(j)+1)))
                count = count+1;
%                 ln = jvec(j)+0.5:0.5:jvec(j+2)+0.5;
%                 xline(1:length(ln),count) = ln;
%                 yline(1:length(ln),count) = repmat(i-0.5,length(ln),1);
                EntMatMask(i,jvec(j)+1:jvec(j+2)) = 0;
                % Ensure subsequent packets are also zeroed
                for i2 = i+1:depths
                    EntMatMask(i2,jvec(j)+1:jvec(j+2)) = 0;
                    
                end
            end
        end
    end
    
%     vmLines = diff(EntMatMask,1,1);
%     hmLines = diff(EntMatMask,1,2);
    if get(handles.BuildLines_checkbox,'value')
        aa = bwboundaries(EntMatMask);
        line(aa{1}(:,2),aa{1}(:,1)+0.5,'linewidth',4)
    end
end

function LoadResults(hObject,handles)

HMV = eval(get(handles.HowManyVolunteers_edit,'string'))';
WG = eval(get(handles.WhichGroup_edit,'string'))';
assert(isequal(size(HMV,1),size(WG,1)),'Length of inputs "How Many Volunteers" and "Which Groups" must be equal')

for i = 1:length(HMV)
    HMVs(i,:) = sprintf('%03d', HMV(i));
    WGs(i,:) = sprintf('%03d', WG(i));
end

if isfield(handles,'Datasets')
    handles = rmfield(handles,'Datasets');
end

for i = 1:size(WG,1)
    handles.Datasets(i,:) = ['Volunteers_' HMVs(i,1:3) '_Group_' WGs(i,1:3)];
end

dataPath = get(handles.BaseDataFolder_edit,'string');

for i = 1:size(handles.Datasets,1)
    h = load(fullfile(dataPath,filesep,handles.Datasets(i,:), filesep,'ResultsSave.mat'));
    hf = fieldnames(h);
    handles.ResultsSaves(i) = h.(hf{1});
end
guidata(hObject,handles)

function BaseDataFolder_edit_Callback(hObject, eventdata, handles)
% hObject    handle to BaseDataFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseDataFolder_edit as text
%        str2double(get(hObject,'String')) returns contents of BaseDataFolder_edit as a double

% --- Executes during object creation, after setting all properties.
function BaseDataFolder_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaseDataFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

[pname] = fileparts(get(hObject,'string'));
handles.GetDataPname = pname;
guidata(hObject,handles)


% --- Executes on button press in BaseDataFolderBrowse_pushbutton.
function BaseDataFolderBrowse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to BaseDataFolderBrowse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pname] = uigetdir;
handles.GetDataPname = pname;
set(handles.SaveSheet_edit,'string',fullfile(pname))
guidata(hObject,handles)


% --- Executes on button press in CrossWPT_Dendro.
function CrossWPT_Dendro_Callback(hObject, eventdata, handles)
% hObject    handle to CrossWPT_Dendro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlotCrossWPT_v2_pushbutton.
function PlotCrossWPT_v2_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCrossWPT_v2_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'ResultsSaves') || isempty(handles.ResultsSaves)
    LoadResults(hObject,handles)
    handles = guidata(hObject);
end

gVals = eval(get(handles.Inconsistencyg_Edit,'string'));
InconVals = eval(get(handles.InconsistencyCutoff_Edit,'string'));

lRS = length(handles.ResultsSaves);

depths = 7;
positions = (2^(depths-1));
indices = positions*2-1;

spR = [1,2,2,2,3,3,4,4,4,4,4,4,5,5,5];
spC = [1,1,2,2,2,2,2,2,3,3,3,3,3,3,3];

for rs = 1:lRS
    
    MIZ = handles.ResultsSaves(rs).MIZ;
    
    Y = zeros(1,((positions*2-1)*(positions*2-2))/2);
    
    count = 0;
    
    for i = 1:(positions*2-2)
        for j = i+1:(positions*2-1)
            count = count+1;
            Y(count) = MIZ(i,j);
        end
    end
    
    Z = linkage(Y,'average');
    
    count=0;
    
    for i = 0:positions*2-2
        count = count+1;
        [dep, pos] = ind2depo(2,i);
        names{count} = ['D' sprintf('%01d',dep) 'P' sprintf('%02d',pos)];
    end
    
    %Plot Dendro
    
    figure(2); clf;
    subplot(spR(lRS),spC(lRS),rs)
    [aaa,bbb,ccc] = dendrogram(Z,indices);
    % set(gcf,'position',[256        -270        1438         501])
    % 478          43        3097         500])
    set(gca,'xticklabel',names(ccc))
    
    rotateticklabel(gca,90);
    texts = findall(2,'type','text');
    gcaPos = get(gca,'position');
    set(gca,'position',gcaPos+[0,0.1,0,-0.1])
    gcaYLim = get(gca,'ylim');
    for i = 1:2:indices; set(texts(i),'position',[indices-(i-1), gcaYLim(1)-diff(gcaYLim)*0.02, 0],'fontsize',10);end
    for i = 2:2:indices; set(texts(i),'position',[indices-(i-1), gcaYLim(1)-diff(gcaYLim)*0.16, 0],'fontsize',10);end
    set(gca,'fontsize',13)
    
    try
        gval = gVals(rs);
        InconVal = InconVals(rs);
    catch
        gval = gVals(1);
        InconVal = InconVals(1);
    end

    Inconsist = inconsistent(Z,gval);
    OrdInc = sort(Inconsist(:,4),'descend');
    minOIc = min(OrdInc(OrdInc~=0))+2*eps;
    maxOIc = max(OrdInc(OrdInc~=0))-3*eps;
    
    InconVal = 1-InconVal;
    
    if InconVal>=0 && InconVal<=1
        OIc = (maxOIc-minOIc)*InconVal+minOIc;
    else
        error('Inconsistency threshold must be a fraction of min/max inconsistency values indexed between 0 and 1')
    end
%     OIc = 1.6584250;
    [T,conn] = cluster(Z,'cutoff',OIc-eps*2,'Depth',gval);
    
    hold on
    h = findobj(gca,'type','line');
    xdata = cell2mat(get(h,'Xdata'));
    ydata = cell2mat(get(h,'Ydata'));
    links = Z(conn==0,3);
    linklocs = ismember(ydata(:,2),links);
    
    if ~get(handles.DeleteLines_checkbox,'value') && isfield(handles,'xdata') && isequal(handles.xdata,xdata)
        handles.xdata = xdata;
        handles.ydata = ydata;
        DendLineCol = get(handles.DendroLineColor_edit,'string');
        rep = length(handles.DendLineCol)+1;
        handles.DendLineCol{rep} = DendLineCol;
        handles.linklocs{rep} = linklocs;
        for lk = 1:size(ydata,1)
            plot(xdata(lk,:),ydata(lk,:),'c','linewidth',1);
        end
        for r = 1:rep
            for lk = find(handles.linklocs{r})'
                plot(xdata(lk,:),ydata(lk,:),handles.DendLineCol{r},'linewidth',1);
            end
        end
    else
        DendLineCol = get(handles.DendroLineColor_edit,'string');
        handles.DendLineCol{1} = DendLineCol;
        for lk = 1:size(ydata,1)
            plot(xdata(lk,:),ydata(lk,:),'c','linewidth',1);
        end
        for lk = find(linklocs)'
            plot(xdata(lk,:),ydata(lk,:),DendLineCol,'linewidth',1);
        end
        handles.xdata = xdata;
        handles.ydata = ydata;
        handles.linklocs{1} = linklocs;
    end
    guidata(hObject,handles)
    hold off
    set(gca,'color',[0 0 0]);
    ndxClust = uint16(T);
    
    ClusterNumbers = double(max(ndxClust));
    if ClusterNumbers == 1
        error('Threshold would return 1 cluster')
    end
    
    maxclst = size(Z,1)+1;
    ddd = zeros(1,maxclst);
    for ii = 1:(maxclst)
        ddd(ii) = find(ccc==bbb(ii));
    end
    DenLocs = ddd;
    
    zloc = Z(size(Z,1)-ClusterNumbers+2,3);
    zlen = size(Z,1);
    
    %Clusters are defined by dendrogram positions
    %ndxClust is reorganized to dendrogram positions
    
    ndxClust1 = ndxClust;
    ClusterXTickLocs = zeros(2,ClusterNumbers);
    for ii = 1:ClusterNumbers
        ClusterXTickLocs(1,ii) = round(mean(DenLocs(ndxClust==ii))); %Tells me where clusters are on the Dendrogram
        ndxClust1(ndxClust==ii) = ClusterXTickLocs(1,ii); % Replaces cluster number with location on Dendrogram
    end
    [ClusterXTickLocs(1,:)] = sort(ClusterXTickLocs(1,:)); %Build 1:ClusterNumbers list from ClusterNumbers distributed across 1:NVoxels numbers
    for ii = 1:ClusterNumbers
        ndxClust(ndxClust1==ClusterXTickLocs(1,ii)) = ii;
        ClusterXTickLocs(2,ii) = ii;
    end
    
    hlen = 200;
    offset = 16;
    %         vrng = [0.70,0.99];
    vrng = [0.60,1];
    srng = [0.85,1];
    
    hmap = linspace(offset/240,(hlen+offset)/240,hlen);
    if ClusterNumbers>hlen
        hmap = repmat(hmap,[ceil(zlen/hlen),1]);
        hmap = hmap(:);
        vmap = linspace(vrng(1),vrng(2),ceil(zlen/hlen/2))';
        vmap = [ones(size(vmap));vmap];
        vmap = repmat(vmap,[hlen,1]);
        smap = linspace(srng(1),srng(2),ceil(zlen/hlen/2))';
        smap = [smap;ones(size(smap))];
        smap = repmat(smap,[hlen,1]);
        newcolormap = [[hmap(1:zlen),vmap(1:zlen),smap(1:zlen)]];
    else
        hmap = hmap(:);
        newcolormap = [[hmap,ones(size(hmap)),ones(size(hmap))]];
    end
    newcolormap = hsv2rgb(newcolormap);
    newcolormap = colorcube(length(newcolormap));
    if get(handles.FixedColormap_checkbox,'value')
        if ClusterNumbers>hlen
            newcolormap = newcolormap([1,ClusterXTickLocs(1,:)+1],:);
        else
            lc = ClusterXTickLocs(1,:);
            if min(lc)>1
                lc = lc-1;
            end
            newcolormap = newcolormap([1,ceil((1-(zlen-lc)/zlen)*hlen)],:);
        end
    else
        newcolormap = newcolormap([1,round(linspace(2,length(newcolormap),ClusterNumbers))],:);
    end
    
    % Plot clusters in pyramidal map
    
    MIZMat = zeros(depths,positions/2);
    vertlines = zeros(sum(2.^(1:depths-1)-1),4);
    ind = 0;
    
    for i = 1:depths
        count = 0;
        jvec = 0:positions/(2^(i-1)):positions;
        for j = 1:(2^(i-1))
            ind = ind+1;
            MIZMat(i,jvec(j)+1:jvec(j+1)) = ndxClust(ind); % Calcs unnormalized entropy
            if i==1; begcap = 0; endcap = 0;
            else
                endcap = 2^(i-1)-1;
                vertlines(begcap+1:begcap+endcap,3:4) = repmat([i-0.5,i+0.5],endcap,1);
            end
            if count > 0
                vertlines(begcap+count,1:2) = [jvec(j)+0.5,jvec(j)+0.5];
            end
            count = count+1;
        end
        begcap = begcap+endcap;
    end
    figure(3); clf;
    subplot(spR(lRS),spC(lRS),rs)
    imagesc(MIZMat);
    colormap(newcolormap)
    colorbar
    
    for i = 1:size(vertlines,1); line(vertlines(i,1:2),vertlines(i,3:4),'color','k','linestyle','-'); end
    % set(gcf,'position',[ 256        -270        1437         438])
    % 721   644   817   353])
    set(gca,'yticklabel',0:6)
    % set(gca,'xtick',get(gca,'xtick')+1)
    % set(gca,'xticklabel',get(gca,'xtick')-1)
    freqs = [1/(0.645*900)*2,1/0.645/2];
    freqs = linspace(freqs(1),freqs(2),64);
    xt = get(gca,'xtick');
    set(gca,'xticklabel',sprintf('%3.2g | ',freqs(xt)))
    % rectangle('position',[0.5,1.5,32,5],'linewidth',5,'linestyle','--')
    
    % line([12.5,12.5,12.5,16.5,16.5,16.5,16.5,32.5,32.5,32.5,32.5,64.5],[7.5,6.5,6.5,6.5,6.5,3.5,3.5,3.5,3.5,2.5,2.5,2.5],'color',[1 0 0],'linestyle','--','linewidth',3)
    set(gca,'fontsize',13)
    
    ylabel('Wavelet packet depth')
    xlabel('Approximate Fourier frequency (Hz)')
    title(num2str(OIc-eps*2))
end
