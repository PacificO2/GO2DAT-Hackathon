% This script is based on "plot_seas_cycle_fminsearch_lineP.m" and
% "get_surface_timeseries_lineP.m". It extracts time series from vertically
% interpolated profiles and fits a seasonal cycle using fitgrp instead of
% fminsearch, taking into account any secular trend. 
% Uses output from: "get_interpolated_profiles_lineP.m"
%
% Input needed: 
% Station
% Variable
% Depth to average
% initial and final year
% 
% 
% History
% December 2020: created by Ana C. Franco (afranco@eoas.ubc.ca), based on
% plot_seas_cycle_fminsearch_lineP.m


%% === Initialize ===
clear all
close all
clc

% === User defined ===

sref = 33;
remin_fact = 117/170;

start_trend = 1990;         % Starting year to calculate trends.
end_trend = 2019;           % Year to end of the trend period.

major_st = [26,20,16,12,4];                      % Select the station to analyze

%variable = [{'dic'}, {'ndic'}, {'arag_2mlr'}, {pH_2mlr}, {temp}, {salt}, {sigma0},{pco2_2mlr},{alk_2mlr},{nalk_2mlr}];
var2plot = 'temp'; % pH_2mlr % revelle_2mlr, 'O2sat'
remove_heatwaves = false;        % Set this to true to remove the springs of 2016 and 2019

save_plot = false;           % Set to false to not save any plot

select_season = 'all_year';   % Select the season to analyse (winter (djf), spring (mam), summer (jja), fall (son))

yaxis2plot = 'pressure';            % Plot the profiles versus 'pressure' or 'sigma0'

plot_anom = 'false';           % Set to true to plo the anomaly of the seasonal cycle instead of the absolute values.

plot_timeseries = false;        % Plot the timeseries figures and standard deviations (heavy and takes time)
save_timeseries_plots = false;  % Save time series plots (only enters here if plot_timeseries set to true)


layers = [10:10:100];

datafrom = 'observations'; % Choose model or observations to do the analysis on the data from line P (interpolated profile) or simulated pactcs30


%% ==== Load the file ===
f1 = figure(1);
for st = 1:length(major_st)
    if strcmp(datafrom,'observations')
        
        datadir = ['/Users/ana/Data/LineP_data/interpolated_profiles/st',num2str(major_st(st)),'/pacifica/QCaug2020/'];
        if strcmp(var2plot,'ndic')
            filename = [datadir,'interpolated_profiles_lineP_1982-2019_',var2plot,num2str(sref),'_st',num2str(major_st(st)),'_interp_pchip_on_',yaxis2plot,'_levels.mat'];
        else
            filename = [datadir,'interpolated_profiles_lineP_1982-2019_',var2plot,'_st',num2str(major_st(st)),'_interp_pchip_on_',yaxis2plot,'_levels.mat'];
        end
        load(filename);
        
    elseif strcmp(datafrom,'model')
        
        datadir = ['/Users/ana/Data/hindcast_pactcs30/colocated_mod_output/'];
        filename = [datadir,'colocated_',var2plot,'_profiles_hc_newsrc_hc05_monthly_pactcs30_lineP.mat'];
        load(filename);
        
    end
    
    
    %% === Treat the data ===
    for l = 1:length(layers)
        
        if strcmp(datafrom,'observations')
            layer_idx = find(interp.interpolated_press == layers(l));
            ydata = interp.interpolated_press_var_profiles(layer_idx,:);
            tdata = floor(interp.datetime2plot);
            
            % Remove nan's from the data and time timeseries to apply fminsearch
            tdata = tdata(~isnan(ydata)); tdata = tdata';
            ydata = ydata(~isnan(ydata)); ydata = ydata';
            
            % Subsample for the selected start and end years.
            start_idx = find(year(tdata) >= start_trend,1,'first');
            end_idx = find(year(tdata) <= end_trend,1,'last');
            
            ydata = ydata(start_idx:end_idx);
            tdata = tdata(start_idx:end_idx); clear start_idx end_idx
            
            
        elseif strcmp(datafrom,'model')
            layer_idx = find(mod_var.depth == layers(l));
            eval(['ydata = squeeze(mod_var.st',num2str(major_st(st)),'(layer_idx,:,:));']);
            [SM,SY] = size(ydata);
            ydata = reshape(ydata,[SM*SY,1]);
            tdata = squeeze(mod_var.dates);
            tdata = datenum(reshape(tdata,[SM*SY,1])); clear SM SY
            
            % Subsample for the selected start and end years.
            start_idx = find(year(tdata) >= start_trend,1,'first');
            end_idx = find(year(tdata) <= end_trend,1,'last');
            
            ydata = ydata(start_idx:end_idx);
            tdata = tdata(start_idx:end_idx); clear start_idx end_idx
            
        end
        
        tfit = [tdata(1):1:tdata(end)]';    % Evaluate function at these dates
        
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
        gpm1 = fitrgp(tdata,ydata,'basisfunction',h1,'beta',beta01,'kernelfunction','squaredexponential','KernelParameters',[300 0.001]);
        
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
        betabar1 = inv(H1*H1')*H1*ydata;
        betacov1 = (sy1^2)*inv(H1*H1');
        
        % 6. Evaluate function with the mean parameters (beta1) to generate the mean
        % response and slope
        yfit = beta1(1) + beta1(2)*tfit + beta1(3)*cos(2*pi*tfit/(365.25)) + beta1(4)*sin(2*pi*tfit/(365.25));
        mean_slope = beta1(1) + beta1(2)*tfit;
        
        % 7. Evaluate function only for the available dates, to match the observed data and calculate residuals
        yfit_sub = beta1(1) + beta1(2)*tdata + beta1(3)*cos(2*pi*tdata/(365.25)) + beta1(4)*sin(2*pi*tdata/(365.25));
        %yfit_sub = A(l)*(tdata-t0) + B(l) + C(l)*(sin(2*pi*(tdata - t0 - D(l))/365));
        
        % 8. Evaluate the function only for one year and without trend
        tfit_seas = datenum('01-jan-1990'):1:datenum('31-jun-1991');
        yfit_seas = beta1(1) + beta1(2)*tfit(1) + beta1(3)*cos(2*pi*tfit_seas/(365.25)) + beta1(4)*sin(2*pi*tfit_seas/(365.25));
        %yfit_seas = B(l) + C(l)*(sin(2*pi*(tfit_seas - t0 - D(l))/365));
        
        
        % 9. Plot seasonal cycle
        f1 = subplot(2,5,st);
        if l<=7
            plot(tfit_seas,yfit_seas(1:length(tfit_seas)),'LineWidth',5);
        else
            plot(tfit_seas,yfit_seas(1:length(tfit_seas)),':','LineWidth',5);
        end
        ylab = [var2plot]; ylabel(ylab);
        hold on
        %
        if strcmp(var2plot,'DIC')
            ylim([1970 2140])% ylim([1970 2040])
            set(gca,'Ydir','reverse')
        elseif strcmp(var2plot,'ndic')
            ylim([1980 2120])%ylim([1980 2090])
            set(gca,'Ydir','reverse')
        elseif strcmp(var2plot,'temp')
            ylim([4 17])
        elseif strcmp(var2plot,'sigma0')
            ylim([23 26])
            set(gca,'Ydir','reverse')
        elseif strcmp(var2plot,'O2')
            ylim([160 320])
        elseif strcmp(var2plot,'pH_2mlr')
            ylim([7.8 8.15])
        elseif strcmp(var2plot,'arag_2mlr') || strcmp(var2plot,'aragonite')
            ylim([1 2.6])
        elseif strcmp(var2plot,'salt')
            ylim([31.4 33.4])
            set(gca,'Ydir','reverse')
        elseif strcmp(var2plot,'NO3')
            ylim([0 25])
            set(gca,'Ydir','reverse')
        end
        t = ['P',num2str(major_st(st)),', ',var2plot,' ', datafrom,' seas. cycle']; title(t)
        legend('10 m','20 m','30 m','40 m','50 m','60 m', '70 m', '80 m','90 m','100 m')
        hAx=gca;                                        % get axes handle
        ixmajor=find(day(tfit_seas)==1);                % indices for first of month in time vector
        hAx.XTick=tfit_seas(ixmajor);                   % set ticks appropriately
        hAx.XAxis.TickLabelRotation=20;                 % and make some room to display
        hAx.XAxis.MinorTickValues=tfit_seas(ixmajor);
        hAx.XMinorGrid='on';
        datetick('x', 'mmm');
        xlim([726834 727380])
        grid on
        
        
        
        %% === Calculate statistics for each parameter ===
        
        % 1. Calculate the means of the parameters
        b(l,st) = beta1(1)+beta1(2)*tfit(1);      % Mean intercept
        A(l,st) = sqrt(beta1(3)^2 + beta1(4)^2);  % Mean amplitude
        s(l,st) = beta1(2)*365.25;                % Mean annual slope
        
        % The variance of each parameter is the diagonal of the betacov1 matrix.
        % 2. Get the standard deviations
        % Intercept
        b_std(l,st) = sqrt(betacov1(1,1));                        % Std dev intercept
        A_std(l,st) = sqrt(betacov1(3,3) + betacov1(4,4));  % Std dev amplitude
        % How to calculate the standard deviation of the amplitude?
        % I think is the squared root of the sum of the variances of the cos and
        % sin component.
        s_std(l,st) = sqrt(betacov1(2,2))*365.25;                 % Std dev annual slope
        
        % 3. Generate 1000 sets of (a0,a1,a2,a3) or (a0,a1) based on the
        % multivariate sampling distribution of these parameters.
        n = 1000;
        R = mvnrnd(beta1,betacov1,n);
        
        if plot_timeseries == true
            % Now generate the fitted model for each of the 1000 distributions
            figure
            for r = 1:n
                % Whole fit
                yfit_tmp(r,:) = R(r,1) + R(r,2)*tfit + R(r,3)*cos(2*pi*tfit/(365.25)) + R(r,4)*sin(2*pi*tfit/(365.25));
                hold on
                plot(tfit',yfit_tmp(r,:),'color',[0.8,0.8,0.8])
                
                % Only slope
                slope_tmp(r,:) =  R(r,1) + R(r,2)*tfit;
                hold on
                plot(tfit',slope_tmp(r,:),'color',[0.8,0.8,0.8])
            end
            
            % plot the data and the fit
            hold on
            plot(tfit,yfit,'r','linewidth',3);          % Add the mean function
            hold on
            plot(tfit,mean_slope,'r','linewidth',3);    % Add the mean slope
            hold on
            plot(tdata,ydata,'bo')                      % Plot the original data
            datetick('x', 'mmmyy');
            
            if save_timeseries_plots == true
                % Need to find directory to save time series
            end
        end
        
        %% === Calculate the residuals ===
        
        res_std(l,st) = nanstd(ydata - yfit_sub);
        seas_std(l,st) = nanstd(yfit);
        total_std(l,st) = nanstd(ydata);
        
    end
    %% === Plot the amplitude vertical profile ===
    
    f1 = subplot(2,5,st+5);
    hold on
    f1 = fill([A(:,st)-A_std(:,st);flipud(A(:,st)+A_std(:,st))],[layers';flipud(layers')],[.7,.7,.7],'linestyle','none','LineWidth',0.2);
    hold on
    h1 = plot(A(:,st),layers,'k.-','LineWidth',2); % Plot the amplitude
    set(gca,'YDir','reverse')
    % %
    % %    hold on
    % %     plot(A*365*30,layers,'b.-'); % Plot the change in 30 years
    hold on
    h2 = plot(res_std(:,st),layers,'r.-','LineWidth',2); % Plot the standard deviation of residuals
    hold on
    h3 = plot(total_std(:,st),layers,'b.-','LineWidth',2); % Plot the total standard deviation
    hold on
    h4 = plot(seas_std(:,st),layers,'g.-','LineWidth',2);
    %     set(gca,'Ydir','reverse')
    if strcmp(var2plot,'O2')
        xlim([0 40]); xlabel('O2, µmol/kg')
    elseif strcmp(var2plot,'DIC')
        xlim([0 50]); xlabel('DIC, µmol/kg')
    elseif strcmp(var2plot,'ndic')
        xlim([0 40]); xlabel('nDIC, µmol/kg')
    elseif strcmp(var2plot,'temp')
        xlim([0 5]); xlabel('Temp, degC')
    elseif strcmp(var2plot,'sigma0')
        xlim([0 1.2]); xlabel('Sigma0, kgm-3')
    elseif strcmp(var2plot,'pH_2mlr')
        xlim([0 0.03]); xlabel('pH 2mlr')
    elseif strcmp(var2plot,'arag_2mlr') || strcmp(var2plot,'aragonite')
        xlim([0 0.4]); xlabel('arag 2mlr')
    elseif strcmp(var2plot,'salt')
        xlim([0 0.5]); xlabel('salt')
    elseif strcmp(var2plot,'NO3')
        xlim([0 5]); xlabel('NO3 µmol/kg')
    end
    ylim([10 100])
    set(gca,'FontSize',16)
    ylabel('Depth (m)')
    legend([h1 h2 h3 h4],'Seasonal amplitude','std dev residuals','total std dev', 'seas std')
    grid on
    %
end


