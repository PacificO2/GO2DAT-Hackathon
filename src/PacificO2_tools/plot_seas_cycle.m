%% Header

%% Initialize
clear variables; close all; clc
%% User defined parameters

start_trend = 2015;             % Starting year to calculate trends.
end_trend   = 2021;             % Year to end of the trend period.

% List of variables
% Temperatura del Aire [?C] -- Velocidad Viento [m s^-1] -- Direcci?n Viento [?]
% Presi?n atmosf?rica [hPa] -- Saturaci?n de Ox?geno [%] -- Saturaci?n de Ox?geno [mL L^-1]
% Ox?geno disuelto [mL L^-1] -- Salinidad de agua [psu] -- Conductividad [S m^-1]
% Temperatura agua [?C] -- Turbidez [NTU] -- Clorofila [mg m^-3]
varname  = {'Air temperature [?C]','Wind velocity [m s^-1]','Wind direction [?]',...
    'Atm pressure [hPa]','Oxygen saturation [%]','Oxygen saturation [mL L^-1]',...
    'Dissolved oxygen [mL L^-1]','Salinity [psu]','Conductivity [S m^-1]',...
    'Temperature [?C]','[NTU]','Chlorophil [mg m^-3]'};
savename = {'air_temp','wind_vel','wind_dir','atm_press','o2sat','o2satml','diss_oxygen'...
    ,'salt','cond','water_temp','turb','chla'};

var2plot = 10;  %

save_plot = false;              % Set to false to not save any plot

select_season = 'all_year';     % Select the season to analyse (winter (djf), spring (mam), summer (jja), fall (son))

yaxis2plot = 'pressure';        % Plot the profiles versus 'pressure' or 'sigma0'

plot_anom = 'false';            % Set to true to plo the anomaly of the seasonal cycle instead of the absolute values.

plot_timeseries = false;        % Plot the timeseries figures and standard deviations (heavy and takes time)
save_timeseries_plots = false;  % Save time series plots (only enters here if plot_timeseries set to true)

%% Load file
load DATATGYDAILY2014202202.mat
%% === Treat the data ===
ydata = TGYdaily(:,var2plot);
tdata = datedTGY;

% Remove nan's from the data and time timeseries to apply fminsearch
tdata = tdata(~isnan(ydata)); tdata = tdata';
ydata = ydata(~isnan(ydata)); ydata = ydata';

% Subsample for the selected start and end years.
start_idx = find(year(tdata) >= start_trend,1,'first');
end_idx = find(year(tdata) <= end_trend,1,'last');

ydata = ydata(start_idx:end_idx);
tdata = tdata(start_idx:end_idx); clear start_idx end_idx

tfit = (tdata(1):1:tdata(end))';    % Evaluate function at these dates
%% === Prepare input for fitgrp ===
% 1. Define function handle
h1 = @(u) [ones(size(u)), u , cos(2*pi*u/(365.25)), sin(2*pi*u/(365.25))];
% For subsurface, isopycnal space
%h1 = @(u) [ones(size(u)), u , cos(2*pi*u/(18.6*365.25)), sin(2*pi*u/(18.6*365.25))];

% 2. Provide initial guess for parameters. Generate them
% in relation to the data in this order: constant, slope, amplitude of the
% cosine and amplitude of the sine terms.
beta01 = [mean(ydata) (max(ydata)-min(ydata))/(max(tdata)-min(tdata)) std(ydata) std(ydata)];

% 3. Calculate the gpm. The 'kernel parameters' are 300 = correlation length
% or link parameter.
% 0.001 = amplitude of the kernel paremeter, has to be many orders of
% magnitude smaller than Y, so that it does not "inform" the posterior
% distribution.
gpm1 = fitrgp(tdata,ydata','basisfunction',h1,'beta',beta01,'kernelfunction','squaredexponential','KernelParameters',[300 0.001]);

% 4. Read beta and sigma
beta1 = gpm1.Beta;
sy1 = gpm1.Sigma;

% H is the set of basis functions evaluated at all training points (i.e., at all available times)
H1 = h1(tdata)';

% 5. Calculate mean and covariance matrices
% This has to do with how to obtain the CI. Betabar1 and beta1 have to be
% the same. According to Adam, this proves that the theory (the formula in
% betabar1) agrees with what matlab calculates.
% betacov1 is the matrix of covariances. The diagonal is the variance of
% each parameter, and it is used to calculate the standard deviation and
% generate random sets of paramters based on the mean and covariance
% matrix.
betabar1 = inv(H1*H1')*H1*ydata';
betacov1 = (sy1^2)*inv(H1*H1');

% 6. Evaluate function with the mean parameters (beta1) to generate the mean
% response and slope
yfit = beta1(1) + beta1(2)*tfit + beta1(3)*cos(2*pi*tfit/(365.25)) + beta1(4)*sin(2*pi*tfit/(365.25));
mean_slope = beta1(1) + beta1(2)*tfit;

% 7. Evaluate function only for the available dates, to match the observed data and calculate residuals
yfit_sub = beta1(1) + beta1(2)*tdata + beta1(3)*cos(2*pi*tdata/(365.25)) + beta1(4)*sin(2*pi*tdata/(365.25));
%yfit_sub = A(l)*(tdata-t0) + B(l) + C(l)*(sin(2*pi*(tdata - t0 - D(l))/365));

% 8. Evaluate the function without trend
tfit_seas = datenum(['01-jan-' num2str(start_trend)]):day(1):datenum(['31-dec-' num2str(end_trend)]);
yfit_seas = beta1(1) + beta1(2)*tfit(1) + beta1(3)*cos(2*pi*tfit_seas/(365.25)) + beta1(4)*sin(2*pi*tfit_seas/(365.25));
%yfit_seas = B(l) + C(l)*(sin(2*pi*(tfit_seas - t0 - D(l))/365));

% gpc = figure('Renderer','painters','Position',[100 100 1000 500]);
% plot(tfit_seas,yfit_seas(1:length(tfit_seas)),'LineWidth',5);
% datetick('x',10)
% ylabel(varname{var2plot},'FontSize',15,'FontWeight','bold')
% xlabel('Time [days]','FontSize',15,'FontWeight','bold')
% title('Fitted seasonal cycle','FontSize',20,'FontWeight','bold')
% print([savename{var2plot},'fit_seas','.png'],'-dpng','-r300')
%% Calculate statistics for each parameter

% 1. Calculate the means of the parameters
b = beta1(1)+beta1(2)*tfit(1);      % Mean intercept
A = sqrt(beta1(3)^2 + beta1(4)^2);  % Mean amplitude
s = beta1(2)*365.25;                % Mean annual slope

% The variance of each parameter is the diagonal of the betacov1 matrix.
% 2. Get the standard deviations
% Intercept
b_std = sqrt(betacov1(1,1));                        % Std dev intercept
A_std = sqrt(betacov1(3,3) + betacov1(4,4));  % Std dev amplitude
% How to calculate the standard deviation of the amplitude?
% I think is the squared root of the sum of the variances of the cos and
% sin component.
s_std = sqrt(betacov1(2,2))*365.25;                 % Std dev annual slope

% 3. Generate 1000 sets of (a0,a1,a2,a3) or (a0,a1) based on the
% multivariate sampling distribution of these parameters.
n = 1000;
R = mvnrnd(beta1,betacov1,n);

gpc = figure('Renderer','painters','Position',[100 100 1000 500]);
for r = 1:n
    % Whole fit
    yfit_tmp(r,:) = R(r,1) + R(r,2)*tfit + R(r,3)*cos(2*pi*tfit/(365.25)) + R(r,4)*sin(2*pi*tfit/(365.25));
    plot(tfit',yfit_tmp(r,:),'color',[0.8,0.8,0.8])
    hold on
    % Only slope
    slope_tmp(r,:) =  R(r,1) + R(r,2)*tfit;
    plot(tfit',slope_tmp(r,:),'color',[0.8,0.8,0.8])
end
axis tight
datetick('x',10)
ylabel(varname{var2plot},'FontSize',15,'FontWeight','bold')
xlabel('Time [days]','FontSize',15,'FontWeight','bold')
print([savename{var2plot},'fit_seas_uncertainty','.png'],'-dpng','-r300')
hold on
% plot the data and the fit
plot(tfit,yfit,'r','linewidth',3);          % Add the mean function
hold on
plot(tfit,mean_slope,'r','linewidth',3);    % Add the mean slope
hold on
plot(tdata,ydata,'bo')    % Plot the original data
hold off
print([savename{var2plot},'fit_seas_uncertainty_and_rawdata','.png'],'-dpng','-r500')
%%