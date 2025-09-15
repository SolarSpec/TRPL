%% Input description
% The workspace should contain the following:
% A row vector of RawData, which is your TRPL decay data
% A row vector of RawTime, which serves as the X axis info. The units are assumed to be in nanoseconds
% A row vector of RawIRF, which is the TRPL decay of the IRF.

% The X axis for Data and IRF need to be the same. All vectors need to be of the same length. Only fit one trace at a time.

Time = RawTime; % Define from raw to avoid issues of different vector lengths when just changing one of them and there are trailling zeros.
IRF = RawIRF;
Data = RawData;

%% Parameters that the user might want to change

% Toggle figure export on or off
FigExport = false;

% Define fit window. Editable from caller via WindowLB, WindowUB
if ~exist('WindowLB','var') || isempty(WindowLB)
    WindowLB = 48;
end
if ~exist('WindowUB','var') || isempty(WindowUB)
    WindowUB = 80;
end

% Define area to get baseline level guess. Editable from caller via WindowBaseline
if ~exist('WindowBaseline','var') || isempty(WindowBaseline)
    WindowBaseline = 900:1100;
end
%% Choose the fit model (now handled by app)
% Check if model is defined by the app
if ~exist('FitModel', 'var') || isempty(FitModel)
    % Fallback to old dialog for standalone script usage
    FitList = {'1 exponential', '1 power law', '1 exponential, 1 power law',...
        '2 exponentials','2 exponentials, 1 power law',...
        '1 exponential, 1 second order','2 exponentials, 1 second order'};
    [ListChoice,tf] = listdlg('ListString',FitList,'SelectionMode','single',...
        'ListSize',[200,100],'Name','Choose fitting model');

    if tf==0
        error("Didn't choose a fit model. Exiting...")
    end

    switch ListChoice
        case 1
            FitModel = '1exp';
            NumExpComponents = 1; NumPowerComponents = 0; NumSecondComponents = 0;
        case 2
            FitModel = '1power';
            NumExpComponents = 0; NumPowerComponents = 1; NumSecondComponents = 0;
        case 3
            FitModel = '1exp1power';
            NumExpComponents = 1; NumPowerComponents = 1; NumSecondComponents = 0;
        case 4
            FitModel = '2exp';
            NumExpComponents = 2; NumPowerComponents = 0; NumSecondComponents = 0;
        case 5
            FitModel = '2exp1power';
            NumExpComponents = 2; NumPowerComponents = 1; NumSecondComponents = 0;
        case 6
            FitModel = '1exp1second';
            NumExpComponents = 1; NumPowerComponents = 0; NumSecondComponents = 1;
        case 7
            FitModel = '2exp1second';
            NumExpComponents = 2; NumPowerComponents = 0; NumSecondComponents = 1;
    end
end

disp(['Fit model chosen: ', FitModel])

%% Start the parallel pool
gcp;


%% Get rid of the trailling 0 points
while Data(end) == 0
    Data(end) = [];
    IRF(end) = [];
    Time(end) = [];
end

%% Adjust WindowBaseline after trimming data
% Ensure WindowBaseline indices are still valid after removing trailing zeros
if exist('WindowBaseline','var') && ~isempty(WindowBaseline)
    dataLength = length(Data);
    % Remove any baseline indices that are now out of bounds
    WindowBaseline = WindowBaseline(WindowBaseline <= dataLength);
    
    % If we lost too many baseline points, create a new safe baseline
    if length(WindowBaseline) < 10
        % Use last 10% of remaining data for baseline
        startIdx = max(1, round(dataLength * 0.9));
        endIdx = dataLength;
        WindowBaseline = startIdx:endIdx;
        fprintf('Warning: Baseline range adjusted due to data trimming. New range: %d:%d\n', startIdx, endIdx);
    end
end

%% Prepare for fit data at the different IRF shifts

% Calculate number of parameters dynamically
% Y offset (1) + Exp components (2 each) + Power components (3 each) + Second components (2 each) + IRF shift (1)
numParam = 1 + NumExpComponents*2 + NumPowerComponents*3 + NumSecondComponents*2 + 1;

%% Restrict to fit window

TimeZero = Time(max(Data)==Data);
if ~exist('ConvPad', 'var') || isempty(ConvPad)
    ConvPad = 15;
end

%Editable. Define how much time to add before and after the window to get rid of the edge effects of the convolution.

WindowData = Data(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);
WindowIRF = IRF(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);
WindowTime = Time(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);

%% Setup initial parameters. 
%Editable for initial guesses and bounds.

% Check if parameters are provided by the app
if exist('InitialGuess', 'var') && exist('LowerBound', 'var') && exist('UpperBound', 'var')
    % Use parameters from app - data-dependent values are already set correctly
    % Just use the provided parameters directly
    
    % Use generalized fitting function
    TRPLfitfun = @(param)Fit_IRF_General(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad,NumExpComponents,NumPowerComponents,NumSecondComponents);
    
else
    % Fallback to original hardcoded parameters for standalone usage
    % First ensure baseline is safe
    if exist('WindowBaseline','var') && ~isempty(WindowBaseline)
        maxIdx = max(WindowBaseline);
        if maxIdx > length(Data)
            dataLength = length(Data);
            startIdx = max(1, round(dataLength * 0.9));
            endIdx = dataLength;
            WindowBaseline = startIdx:endIdx;
            fprintf('Warning: Baseline indices adjusted for fallback. New range: %d:%d\n', startIdx, endIdx);
        end
        baselineValue = mean(Data(WindowBaseline));
    else
        dataLength = length(Data);
        startIdx = max(1, round(dataLength * 0.9));
        endIdx = dataLength;
        baselineValue = mean(Data(startIdx:endIdx));
        fprintf('Warning: No baseline defined for fallback. Using last 10%% of data.\n');
    end
    
    switch FitModel
        case '1exp'         % (1 exponential parameters)
            InitialGuess = [baselineValue; max(Data); 5; 0];
            LowerBound = [0; 0; 0; -3];
            UpperBound = [baselineValue*5; inf; inf; 3];
            TRPLfitfun = @(param)Fit_IRF_1exp(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '1power'       % (1 power law parameters)
            InitialGuess = [baselineValue; max(Data)/2; 0.1; 0.5; 0];
            LowerBound = [0; 0; 0; 0; -3];
            UpperBound = [baselineValue*5; inf; inf; 5; 3];
            TRPLfitfun = @(param)Fit_IRF_1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '1exp1power'       % (1 exp, 1 power law parameters)
            InitialGuess = [baselineValue; max(Data)/2; 2 ; max(Data)/2; 0.765; 1.75; 0];
            LowerBound = [0; 0; 0; 0; 0.765; 1.75; -3];
            UpperBound = [baselineValue*5; inf; 2*max(Time) ; inf; 0.765; 1.75; 3];
            TRPLfitfun = @(param)Fit_IRF_1exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '2exp'         % (2 exponential parameters)
            InitialGuess = [baselineValue; max(Data)/2; 0.5 ; max(Data)/2; 50; 0];
            LowerBound = [0; 0; 0; 0; 0; -3];
            UpperBound = [baselineValue*5; inf; max(Time)/2; inf; 2*max(Time); 3];
            TRPLfitfun = @(param)Fit_IRF_2exp(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '2exp1power'   % (2 exp, 1 power law parameters)
            InitialGuess = [baselineValue; max(Data)/2; 0.5 ; max(Data)/2; 2; max(Data); 0.88; 0.5; 0];
            LowerBound = [0; 0; 0; 0; 0; 0; 0.86; 1.74; -3];
            UpperBound = [baselineValue*5; inf; max(Time)/2; inf; 2*max(Time); inf; 0.89; 1.76; 3];
            TRPLfitfun = @(param)Fit_IRF_2exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '1exp1second'  % (1 exp, 1 second order components)
            InitialGuess = [baselineValue; max(Data)/2; 2 ; max(Data); 0.1; 0];
            LowerBound = [0; 0; 0.2; 0; 0.001; -3];
            UpperBound = [baselineValue*5; inf; 2*max(Time); inf; inf; 3];
            TRPLfitfun = @(param)Fit_IRF_1exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

        case '2exp1second'    % (2 exp, 1 power law parameters)
            InitialGuess = [baselineValue; max(Data)/2; 0.4 ; max(Data)/2; 0.8; max(Data); 0.01; 0];
            LowerBound = [0; 0; 0.2; 0; 0.2; 0; 0.001; -3];
            UpperBound = [baselineValue*5; inf; 8; inf; 8; inf; 0.03; 3];
            TRPLfitfun = @(param)Fit_IRF_2exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);
    end
end

%% Global minimization

options = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','objective',TRPLfitfun,'x0',InitialGuess,...
    'lb',LowerBound,'ub',UpperBound);
problem.options.Display = 'none';
% %problem.options.PlotFcns = @optimplotfval;
% %problem.options.PlotFcns = @optimplotx;

%[x,fval,eflag,output] = fmincon(problem); %Local minimization for testing problem structure.

ms = MultiStart;
ms.UseParallel = true;
ms.StartPointsToRun = 'all';
ms.XTolerance = 1E-8;   
ms.FunctionTolerance = 1E-8;
numMultiStartRuns = 50;

%gs = GlobalSearch(ms);
%rng(14,'twister')   % for reproducibility
tic                 % for timing
%[paramgs,fvalgs] = run(gs,problem);
[paramgs,fvalgs] = run(ms,problem,numMultiStartRuns);
toc


%% Generate best fit data

% Check if using generalized function (from app)
if exist('NumExpComponents', 'var') && (NumExpComponents + NumPowerComponents + NumSecondComponents > 0) && ...
   exist('InitialGuess', 'var') && exist('LowerBound', 'var') && exist('UpperBound', 'var')
    % Use generalized decay function
    [FinalFit, componentData] = Decay_IRF_General(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad,NumExpComponents,NumPowerComponents,NumSecondComponents);
    
    % Extract individual components for plotting
    compIdx = 1;
    if NumExpComponents >= 1
        FitYdataExp = componentData{compIdx};
        compIdx = compIdx + 1;
    end
    if NumExpComponents >= 2
        FitYdataExp2 = componentData{compIdx};
        compIdx = compIdx + 1;
    end
    if NumPowerComponents >= 1
        FitYdataPower = componentData{compIdx};
        compIdx = compIdx + 1;
    end
    if NumSecondComponents >= 1
        FitYdataSecond = componentData{compIdx};
        compIdx = compIdx + 1;
    end
    
else
    % Use original specific functions for backward compatibility
    switch FitModel
        case '1exp'
            [FinalFit,FitYdataExp] = Decay_IRF_1exp(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '1power'
            [FinalFit, FitYdataPower] = Decay_IRF_1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '1exp1power'
            [FinalFit, FitYdataExp, FitYdataPower] = Decay_IRF_1exp1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '2exp'
            [FinalFit, FitYdataExp, FitYdataExp2] = Decay_IRF_2exp(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '2exp1power'
            [FinalFit, FitYdataExp, FitYdataExp2, FitYdataPower] = Decay_IRF_2exp1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '1exp1second'
            [FinalFit, FitYdataExp, FitYdataSecond] = Decay_IRF_1exp1second(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
        case '2exp1second'
            [FinalFit, FitYdataExp, FitYdataExp2, FitYdataSecond] = Decay_IRF_2exp1second(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    end
end

% Calculate the Chi2 from within the fit window, which is what is outputted
% by the Decay functions.Keep in mind that the WindowTime, WindowIRF, and
% WindowData have 'padding' defined by ConvPad.

% 1) build the window mask
winMask = (WindowTime > WindowTime(1)+ConvPad) & (WindowTime < WindowTime(end)-ConvPad);

% 2) extract the vectors you need
timeWin = WindowTime(winMask);
dataWin = WindowData(winMask);
fitWin  = FinalFit;   % FinalFit is already filtered inside your Decay function

% 3) sanity-check their lengths
if numel(timeWin) ~= numel(fitWin)
    error('WindowTime (%d) and FinalFit (%d) lengths do not match.', numel(timeWin), numel(fitWin));
end

% 4) drop any zeros in the fit
keep   = fitWin ~= 0;
CleanTime     = timeWin(keep);
CleanData     = dataWin(keep);
CleanFinalFit = fitWin(keep);

FitDelta = CleanData - CleanFinalFit;
FitRes = FitDelta./sqrt(CleanData);
ChiSq = sum(FitDelta.^2./CleanFinalFit);
RedChiSq = ChiSq/(length(CleanData)-length(InitialGuess)-1);        % Length Initial Guess gives you the number of fitting parameters.


%% Plot ouputs
subplot(10,1,1:6)

semilogy(Time,Data)
hold on
semilogy(CleanTime,CleanFinalFit,'LineWidth',0.75)
yline(paramgs(1),'--k')
ylabel('Counts')
fontsize(gca,12,"points")

% Dynamic plotting based on components
legendEntries = {'Data','Total Fit','Baseline'};

% Check if using generalized function (from app)
if exist('NumExpComponents', 'var') && (NumExpComponents + NumPowerComponents + NumSecondComponents > 0) && ...
   exist('InitialGuess', 'var') && exist('LowerBound', 'var') && exist('UpperBound', 'var')
    
    % Plot exponential components
    for i = 1:NumExpComponents
        if i == 1 && exist('FitYdataExp', 'var')
            plot(CleanTime,FitYdataExp)
            if NumExpComponents == 1
                legendEntries{end+1} = 'Exponential comp.';
            else
                legendEntries{end+1} = 'Exponential comp. 1';
            end
        elseif i == 2 && exist('FitYdataExp2', 'var')
            plot(CleanTime,FitYdataExp2)
            legendEntries{end+1} = 'Exponential comp. 2';
        end
    end
    
    % Plot power law components
    for i = 1:NumPowerComponents
        if i == 1 && exist('FitYdataPower', 'var')
            plot(CleanTime,FitYdataPower)
            if NumPowerComponents == 1
                legendEntries{end+1} = 'Power law comp.';
            else
                legendEntries{end+1} = 'Power law comp. 1';
            end
        end
    end
    
    % Plot second order components
    for i = 1:NumSecondComponents
        if i == 1 && exist('FitYdataSecond', 'var')
            plot(CleanTime,FitYdataSecond)
            if NumSecondComponents == 1
                legendEntries{end+1} = 'Second ord. comp.';
            else
                legendEntries{end+1} = 'Second ord. comp. 1';
            end
        end
    end
    
else
    % Use original plotting for backward compatibility
    switch FitModel
        case '1exp'
            plot(CleanTime,FitYdataExp)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp.'};

        case '1power'
            plot(CleanTime,FitYdataPower)
            legendEntries = {'Data','Total Fit','Baseline','Power law comp.'};

        case '1exp1power'
            plot(CleanTime,FitYdataExp)
            plot(CleanTime,FitYdataPower)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp.','Power law comp.'};

        case '2exp'
            plot(CleanTime,FitYdataExp)
            plot(CleanTime,FitYdataExp2)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp. 1','Exponential comp. 2'};
            
        case '2exp1power'
            plot(CleanTime,FitYdataExp)
            plot(CleanTime,FitYdataExp2)
            plot(CleanTime,FitYdataPower)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp. 1','Exponential comp. 2','Power law comp.'};

        case '1exp1second'
            plot(CleanTime,FitYdataExp)
            plot(CleanTime,FitYdataSecond)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp.','Second ord. comp.'};

        case '2exp1second'
            plot(CleanTime,FitYdataExp)
            semilogy(CleanTime,FitYdataExp2)
            semilogy(CleanTime,FitYdataSecond)
            legendEntries = {'Data','Total Fit','Baseline','Exponential comp. 1','Exponential comp. 2','Second ord. comp.'};
    end
end

legend(legendEntries,'Box','off','FontSize',8,'NumColumns',2)
xlim([WindowLB WindowUB])

hold off
subplot(10,1,8:10)
plot(CleanTime,FitRes)
yline(0)
xlim([WindowLB WindowUB])
xlabel('Time (ns)')
ylabel('Residuals')
fontsize(gca,12,"points")

if FigExport == false
    title(['IRF shift of ' num2str(paramgs(end),3), '. Chi-squared = ' num2str(ChiSq,3), '. Reduced Chi-squared = ' num2str(RedChiSq,3)])
else
    [Figfile,Figpath] = uiputfile('.emf','Export figure to...');
    exportgraphics(gcf,fullfile(Figpath,Figfile))
end

%% Output additional parameters from the fit

% Find t50% from fit
SearchTime = CleanTime(CleanTime>TimeZero);
SearchData = CleanFinalFit(CleanTime>TimeZero);

HalfIntTimes = SearchTime(min(abs(SearchData-max(SearchData)/2))==abs(SearchData-max(SearchData)/2));
RoughHalfIntTime = mean(HalfIntTimes);              %In case there are multiple values

DataSpacing = CleanTime(end)-CleanTime(end-1);
FineTime = RoughHalfIntTime-2*DataSpacing:DataSpacing/50:RoughHalfIntTime+2*DataSpacing;
FineData = interp1(CleanTime,CleanTime,FineTime);
FinalHalfIntTime = FineTime(min(abs(FineData-max(SearchData)/2))==abs(FineData-max(SearchData)/2));

t50 = FinalHalfIntTime(1) - TimeZero;              %TimeZero is defined as the time where the Ydata is max


% Find relative emission from the different components
switch FitModel
    case '1exp'
        RelEmExp = 1;           % Trivial since there is only 1 component.

    case '1power'
        RelEmPower = 1;         % Trivial since there is only 1 component.

    case '1exp1power' 
        ExpCompIntegral = trapz(CleanTime,FitYdataExp);
        PowerCompIntegral = trapz(CleanTime,FitYdataPower);

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(ExpCompIntegral+PowerCompIntegral);

    case '2exp'
        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral);

    case '2exp1power'
        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);
        PowerCompIntegral = trapz(CleanTime,FitYdataPower);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);

    case '1exp1second'
        ExpCompIntegral = trapz(CleanTime,FitYdataExp);
        SecondCompIntegral = trapz(CleanTime,FitYdataSecond);

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(ExpCompIntegral+SecondCompIntegral);

    case '2exp1second'

        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);
        SecondCompIntegral = trapz(CleanTime,FitYdataSecond);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);

end

%% Setup fit function (1 exponential)

function RedChiSq = Fit_IRF_1exp(param,Time,Data,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
IRFshift = param(4);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 power law)

function RedChiSq = Fit_IRF_1power(param,Time,Data,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
PowerAmp = param(2);        %Amplitude of the Power law component
PowerOnset = param(3);      %Power law component onset parameter
PowerAlpha = param(4);      %Power law exponent parameter (alpha)
IRFshift = param(5);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataPL;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 exp, 1 power law)

function RedChiSq = Fit_IRF_1exp1power(param,Time,Data,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
PowerAmp = param(4);        %Amplitude of the Power law component
PowerOnset = param(5);      %Power law component onset parameter
PowerAlpha = param(6);      %Power law exponent parameter (alpha)
IRFshift = param(7);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataExp + YdataPL;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exponentials)

function RedChiSq = Fit_IRF_2exp(param,Time,Data,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
IRFshift = param(6);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exp, 1 power law)

function RedChiSq = Fit_IRF_2exp1power(param,Time,Data,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
PowerAmp = param(6);        %Amplitude of the Power law component
PowerOnset = param(7);      %Power law component onset parameter
PowerAlpha = param(8);      %Power law exponent parameter (alpha)
IRFshift = param(9);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2 + YdataPL;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 exp, 1 second order)

function RedChiSq = Fit_IRF_1exp1second(param,Time,Data,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
SecondAmp = param(4);       %Amplitude of the second order component
SecondKin = param(5);       %Second order kinetic parameter
IRFshift = param(6);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for Second Order Component
TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataSecond;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exp, 1 second order)

function RedChiSq = Fit_IRF_2exp1second(param,Time,Data,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
SecondAmp = param(6);       %Amplitude of the second order component
SecondKin = param(7);       %Second order kinetic parameter
IRFshift = param(8);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Second Order Component
TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2 + YdataSecond;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup decay function for getting the convolved fitted data (1 exponential)

function [TotalOutput,ExpComp] = Decay_IRF_1exp(param,Time,IRF,TimeZero,ConvPad)
% Monoexponential fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
IRFshift = param(4);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Calculate output

TotalOutput = Y0 + ExpComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 power)

function [TotalOutput,PowerComp] = Decay_IRF_1power(param,Time,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
PowerAmp = param(2);        %Amplitude of the Power law component
PowerOnset = param(3);      %Power law component onset parameter
PowerAlpha = param(4);      %Power law exponent parameter (alpha)
IRFshift = param(5);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + PowerComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 power)

function [TotalOutput,ExpComp,PowerComp] = Decay_IRF_1exp1power(param,Time,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
PowerAmp = param(4);        %Amplitude of the Power law component
PowerOnset = param(5);      %Power law component onset parameter
PowerAlpha = param(6);      %Power law exponent parameter (alpha)
IRFshift = param(7);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + PowerComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exponentials)

function [TotalOutput,ExpComp,ExpComp2] = Decay_IRF_2exp(param,Time,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
IRFshift = param(6);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 power)

function [TotalOutput,ExpComp,ExpComp2,PowerComp] = Decay_IRF_2exp1power(param,Time,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
PowerAmp = param(6);        %Amplitude of the Power law component
PowerOnset = param(7);      %Power law component onset parameter
PowerAlpha = param(8);      %Power law exponent parameter (alpha)
IRFshift = param(9);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2 + PowerComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 second)

function [TotalOutput,ExpComp,SecondComp] = Decay_IRF_1exp1second(param,Time,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
SecondAmp = param(4);       %Amplitude of the second ordercomponent
SecondKin = param(5);       %Second order component kinetic parameter
IRFshift = param(6);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for second order component
TempYSecond = SecondAmp./(1+(t-TimeZero)/SecondKin);
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
SecondComp = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + SecondComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
SecondComp = SecondComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 second)

function [TotalOutput,ExpComp,ExpComp2,SecondComp] = Decay_IRF_2exp1second(param,Time,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
SecondAmp = param(6);       %Amplitude of the second ordercomponent
SecondKin = param(7);       %Second order component kinetic parameter
IRFshift = param(8);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for second order component
TempYSecond = SecondAmp./(1+(t-TimeZero)/SecondKin);
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
SecondComp = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2 + SecondComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
SecondComp = SecondComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Generalized fitting function for arbitrary component combinations

function RedChiSq = Fit_IRF_General(param,Time,Data,IRF,TimeZero,ConvPad,NumExp,NumPower,NumSecond)
% Generalized fit function with experimental IRF convolution
% Supports arbitrary numbers of exponential, power law, and second order components

% Set up parameters
t = Time;
paramIdx = 1;

% Y offset (always first parameter)
Y0 = param(paramIdx);
paramIdx = paramIdx + 1;

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

% IRF shift (always last parameter)
IRFshift = param(end);

%Get the IRF values at the shifted time
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Initialize total signal
GuessKin = Y0;

% Process exponential components
for i = 1:NumExp
    ExpAmp = param(paramIdx);
    ExpTau = param(paramIdx+1);
    paramIdx = paramIdx + 2;
    
    % Generate exponential component
    TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau);
    TempYExp(t<TimeZero) = 0;
    ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
    YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);
    
    GuessKin = GuessKin + YdataExp;
end

% Process power law components
for i = 1:NumPower
    PowerAmp = param(paramIdx);
    PowerOnset = param(paramIdx+1);
    PowerAlpha = param(paramIdx+2);
    paramIdx = paramIdx + 3;
    
    % Generate power law component
    TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
    TempYPL(t<TimeZero) = 0;
    ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
    YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);
    
    GuessKin = GuessKin + YdataPL;
end

% Process second order components
for i = 1:NumSecond
    SecondAmp = param(paramIdx);
    SecondKin = param(paramIdx+1);
    paramIdx = paramIdx + 2;
    
    % Generate second order component
    TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
    TempYSecond(t<TimeZero) = 0;
    ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
    YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);
    
    GuessKin = GuessKin + YdataSecond;
end

% Calculate Chi squared taking into consideration the padding
GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Generalized decay function for getting the convolved fitted data

function [TotalOutput, componentData] = Decay_IRF_General(param,Time,IRF,TimeZero,ConvPad,NumExp,NumPower,NumSecond)
% Generalized decay function with experimental IRF convolution
% Returns total fit and individual components

% Set up parameters
t = Time;
paramIdx = 1;

% Y offset (always first parameter)
Y0 = param(paramIdx);
paramIdx = paramIdx + 1;

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

% IRF shift (always last parameter)
IRFshift = param(end);

%Get the IRF values at the shifted time
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Initialize total output (full length first)
TotalOutput = Y0;
componentData = {};

% Process exponential components
for i = 1:NumExp
    ExpAmp = param(paramIdx);
    ExpTau = param(paramIdx+1);
    paramIdx = paramIdx + 2;
    
    % Generate exponential component
    TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau);
    TempYExp(t<TimeZero) = 0;
    ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
    YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);
    
    % Add to total (before filtering)
    TotalOutput = TotalOutput + YdataExp;
    
    % Apply window filtering for component storage
    YdataExpFiltered = YdataExp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
    componentData{end+1} = YdataExpFiltered;
end

% Process power law components
for i = 1:NumPower
    PowerAmp = param(paramIdx);
    PowerOnset = param(paramIdx+1);
    PowerAlpha = param(paramIdx+2);
    paramIdx = paramIdx + 3;
    
    % Generate power law component
    TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
    TempYPL(t<TimeZero) = 0;
    ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
    YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);
    
    % Add to total (before filtering)
    TotalOutput = TotalOutput + YdataPL;
    
    % Apply window filtering for component storage
    YdataPLFiltered = YdataPL(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
    componentData{end+1} = YdataPLFiltered;
end

% Process second order components
for i = 1:NumSecond
    SecondAmp = param(paramIdx);
    SecondKin = param(paramIdx+1);
    paramIdx = paramIdx + 2;
    
    % Generate second order component
    TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
    TempYSecond(t<TimeZero) = 0;
    ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
    YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);
    
    % Add to total (before filtering)
    TotalOutput = TotalOutput + YdataSecond;
    
    % Apply window filtering for component storage
    YdataSecondFiltered = YdataSecond(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
    componentData{end+1} = YdataSecondFiltered;
end

% Apply window filtering to total output (do this LAST)
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end