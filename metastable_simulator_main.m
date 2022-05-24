function [] = metastable_simulator_main(varargin)
% metastable_simulator is a simple investigation tool and does not
% calculate the metastability using our calibration method. 
%
% Pierre Lanari (SWOT 24.05.2022)
%

close all
clc

version = '1.0';

% -- Save results in a separate folder (LastResults)
if isequal(exist([cd,'/LastResults']),7)
    Reply = questdlg('Do you want to replace LastResults?','Modeling','Yes','No (cancel)','Yes');
    switch Reply
        case 'Yes'
            rmdir([cd,'/LastResults'],'s'); 
        otherwise
            return
    end
end
[Success,Message,MessageID] = mkdir('LastResults');
% --

clc, close all
disp(' ')
disp(['    ------------------------------------------'])
disp(['   |           METASTABILITY SIMULATOR        |']);
disp(['   |                    (',version,')                 |']);
disp(['    ------------------------------------------'])
disp(' ')

% -- Activate Fortran libraries
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin')
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')
% --

diary LastResults/output.txt

if isempty(varargin)
    InputFile = 'MetastabilitySimulatorIN.txt';
else
    InputFile = varargin;
end

Job = readInputFile(InputFile);

% Check for database in working directory
if ~exist(fullfile(cd,Job.Database))
    disp([' *** WARNING *** Database ',Job.Database,' not found'])
    disp(' '),disp(' '),disp(' ')
    disp('        -------------- END -------------- ')
    disp(' ')
    return
end

LastStable.ID = [];
LastStable.Minerals = {};
LastStable.Elem = {};
LastStable.MOLES = [];
LastStable.Gsys = [];


for iStep = 1:size(Job.PT,1)
    
    disp(['Calculating step #',num2str(iStep),'/',num2str(size(Job.PT,1))])

    % (1) Equilibrium calculation (minimum G)
    dlmwrite('THERIN',char( ['    ',char(num2str(Job.PT(iStep,1))),'     ',char(num2str(Job.PT(iStep,2)))],['1    ',Job.Bulk,'   * '] ),'delimiter','');

    dlmwrite('XBIN',char(Job.Database,'no'),'delimiter','');

    [wum,yum]=system([Job.PathTher,'   XBIN   THERIN']);
    
    [WorkVariMod] = Core_ReadResTheriak(yum,'');
    
    if isequal(Job.PT(iStep,3),1)
        LastStable.ID = iStep;
        LastStable.Minerals = WorkVariMod.Names4Moles;
        LastStable.Elem = WorkVariMod.Els;
        LastStable.MOLES = WorkVariMod.MOLES;
        LastStable.Gsys = WorkVariMod.Gsys;
        
        GsytMeta(iStep) = WorkVariMod.Gsys;
        GsysEqui(iStep) = WorkVariMod.Gsys;

    else
        
        % Recalculate the G of the metastable system one-by-one...
        GminMeta = zeros(size(LastStable.Minerals));
        for i = 1:length(LastStable.Minerals)-1

            TempBulk = GenerateBulkForMetastablePhase(LastStable.Elem,LastStable.MOLES(i,:));
            
            dlmwrite('THERIN',char( ['    ',char(num2str(Job.PT(iStep,1))),'     ',char(num2str(Job.PT(iStep,2)))],['1    ',TempBulk,'   * '] ),'delimiter','');
            dlmwrite('XBIN',char([Job.Database,'   ',LastStable.Minerals{i}],'no'),'delimiter','');

            [wum,yum]=system([Job.PathTher,'   XBIN   THERIN']);

            [WorkVariMod_META] = Core_ReadResTheriak(yum,'');

            GminMeta(i) = WorkVariMod_META.Gsys;
        end
        
        GsytMeta(iStep) = sum(GminMeta);
        GsysEqui(iStep) = WorkVariMod.Gsys;
    end

end

figure, 

p1 = nexttile; hold on
plot(Job.PT(:,1),GsysEqui/1e3,'o-k')
plot(Job.PT(:,1),GsytMeta/1e3,'o-r')
xlabel('T (°C)')
ylabel('G_s_y_s (kJ/mol)')
legend({'G_s_y_s at equilibrium','Gsys metastable'})

nexttile
plot(Job.PT(:,1),GsytMeta-GsysEqui,'o-b')
xlabel('T (°C)')
ylabel('A (J/mol)')

nexttile
plot(Job.PT(:,1),(GsysEqui(1)-GsysEqui)./1e3,'o-b')
xlabel('T (°C)')
ylabel('G (kJ/mol)')







keyboard










keyboard



end

function [TempBulk] = GenerateBulkForMetastablePhase(Elem,Moles)

TempBulk = '';
for i = 1:length(Elem)
    if Moles(i) > 0
        TempBulk = [TempBulk,char(Elem{i}),'(',num2str(Moles(i),'%.6f'),')'];
    end
end

end



function [WorkVariMod] = Core_ReadResTheriak(OutputTheriakd,ListRefMiner)
% 
% New version (21.03.19)
%
% -> This version reads the chemical potential from the Theriak's output
% and does not required to define the components in the Job.Database.
%

TestInput = strread(OutputTheriakd,'%s','delimiter','\n');

% To display Theriak output: 
%OutputTheriakd

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (0) Gsystem from theriak
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
WhereGsystem = find(ismember(TestInput,'equilibrium assemblage:'))+7;

TheStr = strread(TestInput{WhereGsystem},'%s');

WorkVariMod(1).Gsys = str2num(TheStr{6});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (1) Elements and test for error with the Job.Database ...
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%
WhereElOrder = find(ismember(TestInput,'elements in stable phases:'))+3;

if isempty(WhereElOrder)
    OutputTheriakd
    error('ERROR IN READING THE THERIAK OUTPUT FOR ELEMENTS (see theriak output above); the table containing the elements in stable phases is not available');
    return
end

El1 = strread(char(TestInput(WhereElOrder)),'%s')';
if length(El1) == 10
    SecondRow = strread(char(TestInput(WhereElOrder+1)),'%s')';
    if length(SecondRow) < length(El1) % Could be a second row...
        El1 = [El1,SecondRow];
    end
end

WorkVariMod(1).Els = El1;
WorkVariMod(1).NbEl = length(El1);

if WorkVariMod(1).NbEl > 10
    Shift = 1;
else
    Shift = 0;
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (2a) Elements in stable phases
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

WhereCompositions = find(ismember(TestInput,'elements per formula unit:'))+1;

Compt = 1; Indice = 1;
while 1
    ShiftRows = (Compt-1)*Shift;
    Temp = strread(char(TestInput(WhereCompositions+Compt+ShiftRows)),'%s')';
    if isempty(Temp)
        break          % OK
    end
    if Shift
        Temp2 = strread(char(TestInput(WhereCompositions+Compt+ShiftRows+1)),'%s')';
    end
    
    NbElem1 = length(str2double(Temp(2:end)));
    ASS_Indice2(Indice) = Compt;
    ASS_Names2{Indice} = Temp{1};
    ASS_COMP2(Indice,1:NbElem1) = str2double(Temp(2:end));
    if Shift
        ASS_COMP2(Indice,NbElem1+1:NbElem1+length(Temp2)) = str2double(Temp2);
    end
    Indice = Indice+1;    
        
    Compt = Compt+1;
end


WhereCompositions = find(ismember(TestInput,'elements in stable phases:'))+3+Shift;

Compt = 1; Indice = 1;
while 1
    ShiftRows = (Compt-1)*Shift;
    Temp = strread(char(TestInput(WhereCompositions+Compt+ShiftRows)),'%s')';
    if isempty(Temp)
        break          % OK
    end
    if Shift
        Temp2 = strread(char(TestInput(WhereCompositions+Compt+ShiftRows+1)),'%s')';
    end
    
    NbElem1 = length(str2double(Temp(2:end)));
    ASS_Indice3(Indice) = Compt;
    ASS_Names3{Indice} = Temp{1};
    ASS_COMP3(Indice,1:NbElem1) = str2double(Temp(2:end));
    if Shift
        ASS_COMP3(Indice,NbElem1+1:NbElem1+length(Temp2)) = str2double(Temp2);
    end
    Indice = Indice+1;    
        
    Compt = Compt+1;
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (2b) Volume and densities of solids
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

WhereVol = find(ismember(TestInput,'volumes and densities of stable phases:'))+4;

Compt = 1; Indice = 1; Skip = 0;
while 1
    Temp = strread(char(TestInput(WhereVol+Compt+Skip)),'%s')';
    
    if ~isequal(length(Temp),11)
        break
    end

    FractPer = str2double(Temp(5));

    if FractPer > 0.0001
        VOL_VolCCM(Compt) = str2double(Temp(4));
        VOL_VolFrac(Compt) = FractPer/100;
        VOL_Dens(Compt) = str2double(Temp(11)); 

        VOL_Names{Compt} = Temp{1};

        Compt = Compt+1;
    else
        Skip = Skip+1;
    end
end


% Check for Gases     *** New 1.4.1 (PL - 10.02.2018) 
WhereGf = WhereVol+Compt+Skip+4;

Temp = strread(char(TestInput(WhereGf)),'%s')';

Compt2 = 1;
WhereGf = WhereGf + 1;

if isequal(char(Temp{1}),'gases') && isequal(char(Temp{2}),'and') && isequal(char(Temp{3}),'fluids')    
    while 1
        Temp = strread(char(TestInput(WhereGf+Compt2)),'%s')';

        if ~isequal(length(Temp),7)
            break
        end
        
        VOLgf_VolCCM(Compt2) = str2double(Temp(4));
        VOLgf_VolFrac(Compt2) = 0;
        VOLgf_Dens(Compt2) = str2double(Temp(7));
        
        VOLgf_Names{Compt2} = Temp{1};
        
        Compt2 = Compt2+1;        
    end
    
else
    VOLgf_VolCCM = [];
    VOLgf_VolFrac = [];
    VOLgf_Dens = [];
    VOLgf_Names = '';
end

% Check if a liquid phase should be added to the solids
DoWeUdpateVol = 0;
VOL_Names = '';
for i=1:length(VOLgf_VolCCM)
    TempName = VOLgf_Names{i};
    TempNameStr = strread(TempName,'%s','delimiter','_');
    if iscell(TempNameStr)
        NameGf = TempNameStr{1};
    else
        NameGf = TempName;
    end
    Ok = find(ismember(ListRefMiner,NameGf));
    
    if Ok
        DoWeUdpateVol = 1;
        %Compt = Compt +1;
        VOL_VolCCM(Compt) = VOLgf_VolCCM(i);
        VOL_VolFrac(Compt) = VOLgf_VolFrac(i);
        VOL_Dens(Compt) = VOLgf_Dens(i); 
    
        VOL_Names{Compt} = TempName;
    end
end

if DoWeUdpateVol
    VOL_VolFrac = VOL_VolCCM/sum(VOL_VolCCM);
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (3) SYNCHRONIZATION (SOLIDS + FLUIDS)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%keyboard

for i=1:length(VOL_Names)
    
    [OK,WhereC] = ismember(VOL_Names{i},ASS_Names2);
    
    if ~OK
        OutputTheriakd
        disp(' ')
        disp(' ')
        disp(['IT SEEMS THAT THE COMPOSITION OF PHASE ',char(VOL_Names{i}),' IS NOT AVAILABLE']) 
        disp('please send the theriak output above with error code (ERR2048) to pierre.lanari@geo.unibe.ch  ...')
        keyboard
    end
    
    WorkVariMod(1).Indice(i) = i;
   	WorkVariMod(1).Names{i} = VOL_Names{i};
    WorkVariMod(1).COMP(i,:) = ASS_COMP2(WhereC,:);
    WorkVariMod(1).VolFrac(i) = VOL_VolFrac(i);
    WorkVariMod(1).Dens(i) = VOL_Dens(i);
end

if ~length(VOL_Names)
    WorkVariMod(1).Indice = [];
    WorkVariMod(1).Names = '';
    WorkVariMod(1).COMP = [];
    WorkVariMod(1).VolFrac = [];
    WorkVariMod(1).Dens = [];
end

% Addition of MOLES --------------------------------------------
for i = 1:length(ASS_Names3)
    WhereUnd = find(ismember(char(ASS_Names3{i}),'_'));
    if isequal(length(WhereUnd),1)
        NameTemp = ASS_Names3{i}(1:WhereUnd(1)-1);
    elseif isequal(length(WhereUnd),2)
        NameTemp = ASS_Names3{i}(1:WhereUnd(2)-1);
        NameTemp(WhereUnd(1)) = ' ';
    else
        NameTemp = char(ASS_Names3{i});
    end
    Names4Moles{i} = NameTemp;
end

WorkVariMod(1).Names4Moles = Names4Moles;
WorkVariMod(1).MOLES = ASS_COMP3;
% ---------------------------------------------------------------

WorkVariMod(1).NbPhases = length(WorkVariMod(1).Names);

% FLUID PL - 07.03.2019
if Compt2 > 1 % then there is at least one fluid phase
    WorkVariMod(1).FluidPhases = Compt2-1;
    WorkVariMod(1).FluidDens = VOLgf_Dens;
    
    for i= 1:Compt2-1
        [OK,WhereC] = ismember(VOLgf_Names{i},ASS_Names2);
        
        % added (PL 25.07.2019):
        TheNameFromTheriak = VOLgf_Names{i};
        WereD = ismember(TheNameFromTheriak,'_');
        switch sum(WereD)
            case 1
                Where = find(WereD);
                WorkVariMod(1).FluidNames{i} = char(TheNameFromTheriak(1:Where-1));
                
            case 2
                % we have to delete the first one (Compatibility with the
                % MELT model of DOUG (considered as fluid)
                Where = find(WereD);
                NameTemp = TheNameFromTheriak(1:Where(1)-1);
                if isequal(NameTemp,'LIQtc6')
                    WorkVariMod(1).FluidNames{i} = NameTemp;
                else
                    % we delete the second one ...
                   WorkVariMod(1).FluidNames{i} = TheNameFromTheriak(1:Where(2)-1);
                end
                
            otherwise
                WorkVariMod(1).FluidNames{i} = VOLgf_Names{i};
        end
        
        WorkVariMod(1).FluidCOMP(i,:) = ASS_COMP2(WhereC,:);
    end
    
    % Check for melt (PL 25.07.2019)
    
    %keyboard
    
    % ---
    
else
    WorkVariMod(1).FluidPhases = 0;
    WorkVariMod(1).FluidDens = [];
    WorkVariMod(1).FluidCOMP = [];
end

% ------------------------------------------------------------------------
% Check Theriak Names and remove EM abbreviations
% ------------------------------------------------------------------------
for i=1:length(WorkVariMod(1).Names) 
    TheNameFromTheriak = WorkVariMod(1).Names{i};
    if ~ismember(TheNameFromTheriak,ListRefMiner)
        WereD = ismember(TheNameFromTheriak,'_');
        %keyboard
        switch sum(WereD)
            case 1
                % we delete it...
                Where = find(WereD);
                WorkVariMod(1).Names{i} = TheNameFromTheriak(1:Where-1);
            case 2
                % we have to delete the first one (Compatibility with the
                % MELT model of DOUG (considered as solid)
                Where = find(WereD);
                NameTemp = TheNameFromTheriak(1:Where(1)-1);
                if isequal(NameTemp,'LIQtc6')
                    WorkVariMod(1).Names{i} = NameTemp;
                else
                    % we delete the second one ...
                    WorkVariMod(1).Names{i} = TheNameFromTheriak(1:Where(2)-1);
                end
            case 3
                disp('Oups, too many underscores in this name contact Pierre Lanari')
                keyboard
        end
    end
end


% ------------------------------------------------------------------------
% Check for phase demixion (not identified in ListRefMiner). 
% ------------------------------------------------------------------------
% Note this is typically caused by flat G function in complex solid 
% solution models from Roger Powell (amphiboles).  
%
% in this case we select the phase with the higher volume fraction for
% comparison with the observation and rename the other phases.
%
%                                                  Pierre Lanari (24.10.16)

for i=1:length(WorkVariMod(1).Names)
    TheName = WorkVariMod(1).Names{i};
    Ind = find(ismember(WorkVariMod(1).Names,TheName));
    
    if length(Ind)>1
        if isequal(TheName,'OMPH') || isequal(TheName,'CPXo')
            % Manual selection of OMPH by Lanari (SWOT 2022) for selecting
            % the highest Jad content... 
            WhereNa = find(ismember(WorkVariMod(1).Els,'NA'));
            Na = WorkVariMod.COMP(Ind,WhereNa);
            [Val,IndSort] = sort(Na',2,'descend');
            for j=2:length(IndSort)
                WorkVariMod(1).Names{Ind(IndSort(j))} = [WorkVariMod(1).Names{Ind(IndSort(j))},num2str(j)];
            end
        elseif isequal(TheName,'ClAMP') || isequal(TheName,'ClAMg')
            % Manual selection of Glaucophase by Lanari (SWOT 2022)
            WhereNa = find(ismember(WorkVariMod(1).Els,'NA'));
            Na = WorkVariMod.COMP(Ind,WhereNa);
            [Val,IndSort] = sort(Na',2,'descend');
            for j=2:length(IndSort)
                WorkVariMod(1).Names{Ind(IndSort(j))} = [WorkVariMod(1).Names{Ind(IndSort(j))},num2str(j)];
            end
        else
            Vols = WorkVariMod.VolFrac(Ind);
            [Val,IndSort] = sort(Vols,2,'descend');

            for j=2:length(IndSort)
                WorkVariMod(1).Names{Ind(IndSort(j))} = [WorkVariMod(1).Names{Ind(IndSort(j))},num2str(j)];
            end
        end
          
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (4) CHEMICAL POTENTIAL OF OXIDES
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% WhereChemCo = find(ismember(TestInput,'chemical potentials of components:'))+3;
% 
% Temp = strread(char(TestInput(WhereChemCo)),'%s')';
% NbComponents = str2num(Temp{end});
% 
% Message = char(TestInput(WhereChemCo+7));
% 
% if isequal(Message,'oxydes probably buffered')
%     WeStoreChemPot = 1;
% else
%     WeStoreChemPot = 0;
% end
% 
% WhereChemRead = find(ismember(TestInput,'component   chem.pot.'))+2;
% 
% for i = 1:NbComponents
%     Temp = strread(char(TestInput(WhereChemRead+i-1)),'%s')';
%     % Sometimes there is n components and n-1 chemical potential displayed
%     if ~isempty(Temp)
%         Oxide = char(Temp{1});
%         Oxide = upper(Oxide(2:end-1));
%         WorkVariMod.ChemComp{i} = Oxide;
%         if WeStoreChemPot
%             WorkVariMod.ChemPot(i) = str2num(Temp{2});
%         else
%             WorkVariMod.ChemPot(i) = NaN;
%         end
%     end
% end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% (5) CHEMICAL POTENTIAL OF ELEMENTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

WhereChemEl = find(ismember(TestInput,'--------------------------------------------------------------------'))+1;

CheckEl = zeros(size(WorkVariMod(1).Els));

WorkVariMod(1).El4ChemPot = WorkVariMod(1).Els;
WorkVariMod(1).ChemPotEl = zeros(size(CheckEl));

Compt3 = 0;
while 1
    Temp = strread(char(TestInput(WhereChemEl+Compt3)),'%s')';

    if ~isequal(length(Temp),6)
        break
    end
    if isequal(char(Temp{3}(1)),'"')
        ElTemp = Temp{3}(2:end-1);
        WhereEl = find(ismember(WorkVariMod(1).Els,ElTemp));
        if WhereEl
            WorkVariMod(1).El4ChemPot{WhereEl} = ElTemp;
            WorkVariMod(1).ChemPotEl(WhereEl) = str2num(Temp{5});
        end
    end
    Compt3 = Compt3+1;
end

WorkVariMod(1).Gsys2 = sum(WorkVariMod(1).ChemPotEl.*WorkVariMod(1).MOLES(end,:));

end



function [Job] = readInputFile(FileName)
%

fid = fopen(FileName,'r');

TheL = fgetl(fid);
TheS = textscan(TheL,'%s');
Job.MetastableSimulatorVersion = char(TheS{1}(2));

TheL = fgetl(fid);
TheS = textscan(TheL,'%s');
Job.PathTher = char(TheS{1}(2));

TheL = fgetl(fid);
TheS = textscan(TheL,'%s');
Job.Database = char(TheS{1}(2));

TheL = fgetl(fid);
TheS = textscan(TheL,'%s');
Job.MetaCalc = char(TheS{1}(2));

TheL = fgetl(fid);
TheS = textscan(TheL,'%s');
Job.Bulk = char(TheS{1}(2));

while 1
    TheLine = fgetl(fid);
    if isequal(TheLine,-1)
        break
    end
    if length(TheLine) > 1
        if isequal(TheLine(1),'>')
            % Read PT path
            CountPT = 0;
            while 1
                TheLine = fgetl(fid);
                if isequal(TheLine,-1)
                    break
                end
                if length(TheLine) > 1
                    CountPT = CountPT + 1;
                    Job.PT(CountPT,:) = strread(TheLine,'%f');
                end            
            end
        end
    end
end     
fid = fclose(fid);

end
