% The MIT License (MIT)
% Copyright © 2023 Christian Stenz
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software”), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


%% Initialize
clearvars

% -------------------------------
% LOAD YOUR XYZ- & CIF-FILE HERE:
% -------------------------------
filename  = 'CIFs\GST225_HT_cub.xyz';
CIFname   = [filename(1:end-4), '.cif'];  % CIF is also required in the same directory
indstring = max(strfind(filename, '\'));
Name      = filename(indstring+1:end-4);

startRow   = 3;
formatSpec = '%2s%12f%12f%f%[^\n\r]';
fileID     = fopen(filename,'r');
% Read columns of data according to the format.
dataArray    = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray{1} = strtrim(dataArray{1}); % Remove white space around all cell columns.
fclose(fileID);
% Create output variable
data = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','VarName2','VarName3','VarName4'});
clearvars filename startRow formatSpec fileID dataArray ans; % Clear temporary variables

UnitCell   = table2array(data(:,2:end));
NumAtoms   = size(data,1);
Species    = table2array(data(:,1));
AllSpecies = unique(Species,'stable');
NumSpecies = size(AllSpecies,1);

if NumSpecies == 1
    S1 = {};
    for i = 1:NumAtoms
        S1 = [S1, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
    end
    S1           = cell2mat(S1');
    UnitCellSpec = {S1};
elseif NumSpecies == 2
    S1 = {};
    S2 = {};
    for i = 1:NumAtoms
        if Species(i) == AllSpecies(1)
            S1 = [S1, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
        elseif Species(i) == AllSpecies(2)
            S2 = [S2, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
        end
    end
    S1           = cell2mat(S1');
    S2           = cell2mat(S2');
    UnitCellSpec = {S1, S2};
elseif NumSpecies == 3
    S1 = {};
    S2 = {};
    S3 = {};
    for i = 1:NumAtoms
        if Species(i) == AllSpecies(1)
            S1 = [S1, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
        elseif Species(i) == AllSpecies(2)
            S2 = [S2, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
        elseif Species(i) == AllSpecies(3)
            S3 = [S3, [UnitCell(i,1),UnitCell(i,2),UnitCell(i,3)]];
        end
    end
    S1           = cell2mat(S1');
    S2           = cell2mat(S2');
    S3           = cell2mat(S3');
    UnitCellSpec = {S1, S2, S3};
end

%% Find basis vectors from CIF
opts = delimitedTextImportOptions("NumVariables", 1);
% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "";
% Specify column names and types
opts.VariableNames = "VarName1";
opts.VariableTypes = "string";
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
% Specify variable properties
opts = setvaropts(opts, "VarName1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName1", "EmptyFieldRule", "auto");
% Import the data
CIF = readmatrix(CIFname, opts);
clear opts

t = CIF(:,:);

cellOutput = 1;
inda = strfind(CIF(:,:),'_cell_length_a','ForceCellOutput',cellOutput);
indb = strfind(CIF(:,:),'_cell_length_b','ForceCellOutput',cellOutput);
indc = strfind(CIF(:,:),'_cell_length_c','ForceCellOutput',cellOutput);
indalpha = strfind(CIF(:,:),'_cell_angle_alpha','ForceCellOutput',cellOutput);
indbeta  = strfind(CIF(:,:),'_cell_angle_beta','ForceCellOutput',cellOutput);
indgamma = strfind(CIF(:,:),'_cell_angle_gamma','ForceCellOutput',cellOutput);

a = 0;
b = 0;
c = 0;
alpha = 0;
beta = 0;
gamma = 0;
for i = 1:size(inda(:,1),1)
    if isempty(inda{i,1}) == 0 && a == 0
        a = i;
    end
    if isempty(indb{i,1}) == 0 && b == 0
        b = i;
    end
    if isempty(indc{i,1}) == 0 && c == 0
        c = i;
    end
    if isempty(indalpha{i,1}) == 0 && alpha == 0
        alpha = i;
    end
    if isempty(indbeta{i,1}) == 0 && beta == 0
        beta = i;
    end
    if isempty(indgamma{i,1}) == 0 && gamma == 0
        gamma = i;
    end
end

Nums = cell(6,1);
vars = {a, b, c, alpha, beta, gamma};
for ix = 1:length(vars)
    
    B = regexp(CIF(vars{ix},:),'\d*','Match');
    for ii = 1:size(B,2)
        if ~isempty(B{ii})
            if ii == 1
                Nums{ix} = str2double(B{ii});
            elseif ii == 2
                Nums{ix} = Nums{ix} + str2double(['0.',B{ii}]);
            end
        else
            Nums{ix} = NaN;
        end
    end
    
end

a = Nums{1,1};
b = Nums{2,1};
c = Nums{3,1};
alpha = Nums{4,1};
beta  = Nums{5,1};
gamma = Nums{6,1};


basevecs = cell(1,3);
if alpha == 90 && beta == 90 && gamma == 90
    basevecs{1,1} = a*[1,0,0];
    basevecs{1,2} = b*[0,1,0];
    basevecs{1,3} = c*[0,0,1];
elseif gamma ~= 90 && alpha == 90 && beta == 90
    basevecs{1,1} = a*[1,0,0];
    basevecs{1,2} = b*[cosd(gamma),sind(gamma),0];
    basevecs{1,3} = c*[0,0,1];
elseif gamma == 90 && alpha == 90 && beta ~= 90
    basevecs{1,1} = a*[1,0,0];
    basevecs{1,2} = b*[0,1,0];
    basevecs{1,3} = c*[cosd(beta),0,sind(beta)];
else
    basevecs{1,1} = a*[1,0,0];
    basevecs{1,2} = b*[cosd(gamma),sind(gamma),0];
    
    Roty = [cosd(-beta)     0           sind(-beta);
            0               1           0;
            -sind(-beta)    0           cosd(-beta)];
    cprime = c*Roty*[1,0,0]';
    xn = cosd(alpha)*b*c-basevecs{1,2}(1)*cprime(1);
    an = basevecs{1,2}(2)*cprime(2)+basevecs{1,2}(3)*cprime(3);
    bn = basevecs{1,2}(2)*cprime(3)-basevecs{1,2}(3)*cprime(2);
    delta = 2*atand((bn-sqrt(an^2+bn^2-xn^2))/(an+xn));
    Rotx = [1   0               0;
            0   cosd(-delta)    -sind(-delta);
            0   sind(-delta)    cosd(-delta)];
    basevecs{1,3} = (Rotx*cprime)';
end

%% Plot Unit Cell

MaxXY = max(max(UnitCell(:,1)),max(UnitCell(:,2)))*1.3;
MaxZ = max(UnitCell(:,3))*1.1;

fig1 = figure('Name','Unit Cell','units','normalized','outerposition',[0.1 0 0.7 1]);
for i = 1:NumSpecies
    scatter3(UnitCellSpec{i}(:,1), UnitCellSpec{i}(:,2), UnitCellSpec{i}(:,3), 50,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[((i-0.9)/(NumSpecies-0.9))*0.9 0.75 (1-(i-0.9)/(NumSpecies-0.9))],...
        'DisplayName',AllSpecies(i));
    hold on
end
grid on
plot3([0,0],[0,0],[0 MaxZ],'-k','HandleVisibility','off','PickableParts','none')
plot3([0,0],[-MaxXY MaxXY],[0,0],'-k','HandleVisibility','off','PickableParts','none')
plot3([-MaxXY MaxXY],[0,0],[0,0],'-k','HandleVisibility','off','PickableParts','none')
xlabel('$x$ (\AA)', 'Interpreter', 'Latex')
ylabel('$y$ (\AA)', 'Interpreter', 'Latex')
zlabel('$z$ (\AA)', 'Interpreter', 'Latex')
legend('show','AutoUpdate','on');
ax = gca;
axis([-MaxXY MaxXY -MaxXY MaxXY 0 MaxZ]);
axis equal
ax.Toolbar.Visible = 'off';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

for i = 1:3
quiver3(0,0,0,basevecs{1,i}(1),basevecs{1,i}(2),basevecs{1,i}(3),'r',...
         'HandleVisibility','off','LineWidth',1.75)
end
saveas(fig1, ['Reduced_Unit_Cells\', Name, '_Unreduced.fig'])

%% Find reduced unit cell

translations = [basevecs{1};basevecs{2};basevecs{3};...
    basevecs{1}+basevecs{2};basevecs{1}+basevecs{3};basevecs{2}+basevecs{3};...
    basevecs{1}-basevecs{2};basevecs{1}-basevecs{3};basevecs{2}-basevecs{3};...
    basevecs{1}+basevecs{2}+basevecs{3};...
    basevecs{1}+basevecs{2}-basevecs{3};...
    basevecs{1}-basevecs{2}+basevecs{3};...
    -basevecs{1}+basevecs{2}+basevecs{3}];

% set under basis vectors translational invariant atoms to [Inf, Inf, Inf]
ReducedUnitCellSpec = UnitCellSpec;
for n = 1:NumSpecies
    for i = 1:size(ReducedUnitCellSpec{1,n},1)
        for j = 1:size(translations,1)
            for sign = [+1,-1]
                index = [];
                index = [index, find(ismembertol(ReducedUnitCellSpec{1,n},...
                    ReducedUnitCellSpec{1,n}(i,:)+sign*translations(j,:), 1e-5, 'ByRows', true)==1)];
                ReducedUnitCellSpec{1,n}(unique(index),1) = [inf];
                ReducedUnitCellSpec{1,n}(unique(index),2) = [inf];
                ReducedUnitCellSpec{1,n}(unique(index),3) = [inf];
            end
        end
    end
end

% delete unneccesary entries ([Inf, Inf, Inf])
for n = 1:NumSpecies
    [r,l] = find(ismember(ReducedUnitCellSpec{1,n},Inf));
    ReducedUnitCellSpec{1,n}(r,:) = [];
end

%% Plot reduced unit cell
fig2 = figure('Name','Reduced Unit Cell','units','normalized','outerposition',[0.1 0 0.7 1]);

title('Reduced Unit Cell')
for i = 1:NumSpecies
    scatter3(ReducedUnitCellSpec{i}(:,1), ReducedUnitCellSpec{i}(:,2), ReducedUnitCellSpec{i}(:,3), 50,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[((i-0.9)/(NumSpecies-0.9))*0.9 0.75 (1-(i-0.9)/(NumSpecies-0.9))],...
        'DisplayName',AllSpecies(i));
    hold on
end
for i = 1:3
    quiver3(0,0,0,basevecs{i}(1),basevecs{i}(2),basevecs{i}(3),'r',...
        'DisplayName',sprintf('Basis Vector %u',i),'LineWidth',1.75,'PickableParts','none')
end
grid on
plot3([0,0],[0,0],[0 MaxZ],'-k','HandleVisibility','off','PickableParts','none')
plot3([0,0],[-MaxXY MaxXY],[0,0],'-k','HandleVisibility','off','PickableParts','none')
plot3([-MaxXY MaxXY],[0,0],[0,0],'-k','HandleVisibility','off','PickableParts','none')
xlabel('$x$ (\AA)', 'Interpreter', 'Latex')
ylabel('$y$ (\AA)', 'Interpreter', 'Latex')
zlabel('$z$ (\AA)', 'Interpreter', 'Latex')
axis equal;
legend('show');

saveas(fig2, ['Reduced_Unit_Cells\', Name, '.fig'])

%% Construct crystal lattice
MaxR = 20; % in Angstroem
MeasStep = 0.01;

% Find number of unit cells to construct lattice
minUCRange = min([norm(cell2mat(basevecs(1,1))),...
    norm(cell2mat(basevecs(1,2))),...
    norm(cell2mat(basevecs(1,3))),...
    norm(cell2mat(basevecs(1,1))+cell2mat(basevecs(1,2))),...
    norm(cell2mat(basevecs(1,2))+cell2mat(basevecs(1,3))),...
    norm(cell2mat(basevecs(1,1))+cell2mat(basevecs(1,3))),...
    norm(cell2mat(basevecs(1,1))+cell2mat(basevecs(1,2))+cell2mat(basevecs(1,3)))]);

maxUCRange = max([norm(cell2mat(basevecs(1,1))),...
    norm(cell2mat(basevecs(1,2))),...
    norm(cell2mat(basevecs(1,3)))]);

IndmaxUCR = find([norm(cell2mat(basevecs(1,1))),...
    norm(cell2mat(basevecs(1,2))),...
    norm(cell2mat(basevecs(1,3)))] == maxUCRange);

CellNo1 = ceil(MaxR/minUCRange)+1;
CellNo2 = ceil(MaxR/minUCRange)+1;
CellNo3 = ceil(MaxR/minUCRange)+1;

switch IndmaxUCR(1,1)
    case 1
        CellNo1 = ceil(MaxR/maxUCRange)+1;
    case 2
        CellNo2 = ceil(MaxR/maxUCRange)+1;
    case 3
        CellNo3 = ceil(MaxR/maxUCRange)+1;
end

empty = cell(NumSpecies,1);
for bv1 = -CellNo1:CellNo1
    for bv2 = -CellNo2:CellNo2
        for bv3 = -CellNo3:CellNo3
            basevec = bv1*basevecs{1} + bv2*basevecs{2} + bv3*basevecs{3};
            for i = 1:NumSpecies
                empty{i} = [empty{i}, (ReducedUnitCellSpec{i}+basevec)'];
            end
        end
    end
    disp(['Step 1/3: Constructing Lattice: ', num2str(round((bv1+CellNo1+1)/(2*CellNo1+1)*100,2)), char(9), '%']);
end

%% Calculate atomic distances with respect to each reference atom in unit cell

% number of correlation functions for NumSpecies atom species:
Ntot = 0;
for n = 1:NumSpecies
    Ntot = Ntot + n;
end

Lattice  = cell(NumSpecies,1);
DistSelf = cell(NumSpecies,1);

% calculate distances between atoms of same type (e.g.: g(TeTe) & g(SbSb))
for n = 1:NumSpecies
    Lattice{n,1} = empty{n,1}';
    NumRefAtoms = size(ReducedUnitCellSpec{1,n},1);
    for k = 1:NumRefAtoms
        for j = 1:size(Lattice{n,1},1)
            DistSelf{n,1}(k,j) = norm(Lattice{n,1}(j,:)-ReducedUnitCellSpec{1,n}(k,:));
        end
        disp(['Step 2/3: Calculating Atomic Distances (', convertStringsToChars(AllSpecies(n)), '-', convertStringsToChars(AllSpecies(n)), '): ', num2str(k/NumRefAtoms*100), char(9), '%']);
    end
end

% calculate distances between atoms of different type (e.g.: g(SbTe))
if NumSpecies > 1
    
    DistMix    = cell(Ntot-NumSpecies,1);
    CombMixPDF = combnk(1:NumSpecies,2);
    
    for nc = 1:size(CombMixPDF,1)
        n1 = CombMixPDF(nc,1);
        n2 = CombMixPDF(nc,2);
        Lattice{n1,1} = empty{n1,1}';
        NumRefAtoms = size(ReducedUnitCellSpec{1,n2},1);
        for k = 1:NumRefAtoms
            for j = 1:size(Lattice{n1,1},1)
                DistMix{nc,1}(k,j) = norm(Lattice{n1,1}(j,:)-ReducedUnitCellSpec{1,n2}(k,:));
            end
            disp(['Step 2/3: Calculating Atomic Distances (', convertStringsToChars(AllSpecies(n1)), '-', convertStringsToChars(AllSpecies(n2)), '): ', num2str(k/NumRefAtoms*100), char(9), '%']);
        end
    end
end

Dist                     = cell(Ntot,1);
Dist(1:NumSpecies,1)     = DistSelf;
if NumSpecies > 1; Dist(NumSpecies+1:end,1) = DistMix; end

%% Plot lattice

fig3 = figure('Name','Lattice','units','normalized','outerposition',[0.1 0 0.7 1]);
for n = 1:NumSpecies
    scatter3(Lattice{n,1}(:,1), Lattice{n,1}(:,2), Lattice{n,1}(:,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[((n-0.9)/(NumSpecies-0.9))*0.9 0.75 (1-(n-0.9)/(NumSpecies-0.9))],...
        'DisplayName',AllSpecies(n),'MarkerFaceAlpha',1.0);
    hold on
end
for i = 1:3
    quiver3(0,0,0,basevecs{i}(1),basevecs{i}(2),basevecs{i}(3),'r',...
        'DisplayName',sprintf('Basis Vector %u',i),'LineWidth',1.75)
end
grid on
plot3([0,0],[0,0],[-MaxZ MaxZ],'-k','HandleVisibility','off')
plot3([0,0],[-MaxXY MaxXY],[0,0],'-k','HandleVisibility','off')
plot3([-MaxXY MaxXY],[0,0],[0,0],'-k','HandleVisibility','off')
xlabel('$x$ (\AA)', 'Interpreter', 'Latex')
ylabel('$y$ (\AA)', 'Interpreter', 'Latex')
zlabel('$z$ (\AA)', 'Interpreter', 'Latex')
axis equal
[X,Y,Z] = sphere();
X=X*MaxR;
Y=Y*MaxR;
Z=Z*MaxR;
surf(X,Y,Z,'DisplayName','Integration Radius')
legend('show');

saveas(fig3, ['Integrated_Lattices\', Name, '.fig'])

%% Count atoms in spherical shells
disp(['Step 3/3: Counting Atoms (1/', num2str(Ntot), '): 0  %'])

g       = cell(Ntot,1);
gscaled = cell(Ntot,1);
gtotal  = cell(Ntot,1);

for n = 1:Ntot
    NumRefAtoms = size(Dist{n,1},1);
    for k = 1:NumRefAtoms
        nri = 1;
        for ri = 1:MeasStep:MaxR
            count = 0;
            for j = 1:size(Dist{n},2)
                if ri <= Dist{n}(k,j) && Dist{n}(k,j) < ri+MeasStep
                    count = count + 1;
                end
            end
            if n <= NumSpecies
                g{n,1}(k,nri) = count;
                gscaled{n,1}(k,nri) = count/(4*pi*ri.^2);
            elseif n > NumSpecies
                g{n,1}(k,nri) = 2*count;
                gscaled{n,1}(k,nri) = 2*count/(4*pi*ri.^2);
            end
            nri = nri + 1;
        end
        disp(['Step 3/3: Counting Atoms (', num2str(n), '/', num2str(Ntot), '): ', num2str(k/NumRefAtoms*100), char(9), '%']);
    end
    
    dummy = 0.*gscaled{n,1}(1,:);
    for k = 1:NumRefAtoms
        dummy = dummy+gscaled{n,1}(k,:)/NumRefAtoms;
    end
    gtotal{n,1} = dummy;
end

%% Save PDF to txt file

r = 1:MeasStep:MaxR;
r = r';
varNames   = {'r'};
tableInput = {r};
for n = 1:Ntot
    if n <= NumSpecies
        varNames = [varNames, [AllSpecies{n}, '-', AllSpecies{n}]];
    else
        nc = n - NumSpecies;
        n1 = CombMixPDF(nc,1);
        n2 = CombMixPDF(nc,2);
        varNames = [varNames, [AllSpecies{n1}, '-', AllSpecies{n2}]];
    end
    tableInput = [tableInput, gtotal{n,1}'];
end
T = table(tableInput{:}, 'VariableNames', varNames);
writetable(T, ['PDFs\', Name, '.txt'])

%% Plot computed PDF
fig4 = figure('Name','Pair Distribution Function','units','normalized','outerposition',[0.1 0 0.7 1]);

Gr = cell2mat(gtotal)';

% plot
b = bar(r, Gr, 'stacked', 'barwidth', 3.0);

% change plot colors and legend labels
for k = 1:NumSpecies
    b(k).FaceColor   = element_color(AllSpecies{k});
    b(k).EdgeColor   = 'k';
    b(k).LineWidth   = 0.1;
    b(k).DisplayName = ['\textsf{', varNames{k+1}, '}'];
end
if Ntot > NumSpecies
    for nc = 1:size(CombMixPDF,1)
        k  = nc + NumSpecies;
        n1 = CombMixPDF(nc,1);
        n2 = CombMixPDF(nc,2);
        b(k).FaceColor   = element_color(AllSpecies{n1});
        b(k).EdgeColor   = element_color(AllSpecies{n2});
        b(k).LineWidth   = 1.3;
        b(k).LineStyle   = '--';
        
        b(k).DisplayName = ['\textsf{', varNames{k+1}, '}'];
    end
end
xlabel('$r$ [\AA]', 'Interpreter', 'LaTex');
ylabel('$G(r)$', 'Interpreter', 'LaTex');
xlim([0 10]);
lgd = legend('show', 'Interpreter', 'Latex');
    
saveas(fig4, ['PDFs\PDF_', Name, '.fig'])


function color = element_color(PSE_Sym)
%% Import VESTA atom colors
filename = 'elements.txt';
formatSpec = '%3f%3s%*6*s%*6*s%*7f%11f%11f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
dataArray{2} = strtrim(dataArray{2});
fclose(fileID);
elements = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','H','VarName6','VarName7','VarName8'});
clearvars filename formatSpec fileID dataArray ans;

Ar  = table2array(elements);
Z   = str2double(Ar(:,1));
El  = Ar(:,2);
RGB = str2double(Ar(:,3:5));

color = [1, 1, 1];
i = 1;
while color == [1, 1, 1]
    if i <= 97 && string(PSE_Sym) == El(i)
        color = [RGB(i,1), RGB(i,2), RGB(i,3)];
    elseif i == 97
        color = [0.1,0.1,0.1];
    end
    i = i + 1;
end
end