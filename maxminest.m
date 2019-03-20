function varargout = maxminest( record, dur_ratio )

% This function, written by Joseph A. Main, is a modified version of a 
% function previously written by Fahim Sadek, based on the procedure 
% described in the following paper:
% 
% Sadek, F. and Simiu, E. (2002). "Peak non-gaussian wind effects for
% database-assisted low-rise building design." Journal of Engineering
% Mechanics, 128(5), 530-539.
%
% The function computes expected maximum and minimum values of time series
% by estimating probability distributions for the peaks.
% 
% The function can be called in the following alternative forms:
%     [max_est min_est] = maxminest(record)
%     [max_est min_est] = maxminest(record, dur_ratio)
%     [max_est min_est max_std min_std] = maxminest(record, dur_ratio)
%
% INPUT ARGUMENTS:
% Each row of "record" is a time series.
% The optional input argument "dur_ratio" allows peaks to be estimated for
% a duration that differs from the duration of the record itself:
%    dur_ratio = [duration for peak estimation]/[duration of record]
%    (If unspecified, a value of 1 is used.)
%
% OUTPUT ARGUMENTS:
% "max_est" gives the expected maximum values of each row of "record"
% "min_est" gives the expected minimum values of each row of "record"
% "max_std" gives the standard deviations of the maximum value for each row of "record"
% "min_std" gives the standard deviations of the minimum value for each row of "record"
%
% "max_est", "min_est", "max_std", and "min_std" are r by 1 vectors, where 
% r is the number of rows in "record".
%
% updated 9 Nov 2006

if nargin==0 || isempty( record )
    error('Time series input expected');
elseif nargin<2 || isempty( dur_ratio )
    dur_ratio = 1; % no duration ratio supplied --> use default of 1
end

if ~ismember(nargout,[2 4])
    error(['Either 2 or 4 output arguments expected (' num2str(nargout) ' requested)']);
end    

plot_on = 0; % turns plotting on (1) [for diagnostics] or off (0)

n_cdf_pk = 1000;     % number of points in CDF for peaks
cdf_pk_min = 0.0005; % minimum value of CDF for peaks
cdf_pk_max = 0.9995; % maximum value of CDF for peaks
cdf_pk = linspace(cdf_pk_min,cdf_pk_max, n_cdf_pk); % linearly spaced CDF values for peaks

rsize = size(record);

max_est = zeros(rsize(1),1);
min_est = zeros(rsize(1),1);
if nargout==4
    max_std = zeros(rsize(1),1);
    min_std = zeros(rsize(1),1);
end

for i=1:rsize(1)
    x = record(i,:);
    n = length(x);

    mean_x = mean( x );
    std_x = std(x, 1); % biased estimate
    skew_x = sum((x - mean_x).^3)/ (n*std_x^3); % biased estimate (Matlab default for 'skewness')
    
    X = x*sign(skew_x);  % Change sign (if necessary) to have positive skewness
    sort_X = sort( X ); % sort in ascending order:
    
    mean_X = mean_x*sign(skew_x);
    std_X = std_x;
   
    CDF_X = (1:n)/(n+1); % Empirical Cumulative Distribution Function
    
    % resample CDF more coarsely for more efficient parameter estimation:
    n_coarse = min([n 1000]);
    CDF_coarse = linspace(1/(n_coarse+1),n_coarse/(n_coarse+1),n_coarse);
    X_coarse = interp1(CDF_X, sort_X, CDF_coarse);
    
    % Estimate shape parameter of gamma distribution from coarse CDF:
    mean_X_coarse = mean(X_coarse);
    std_X_coarse = std(X_coarse);

    gamma_min = 1;
    gamma_max = 125;
    n_gamma = 19; % number of trial values of gamma
    n_start = 8; % start checking with this index in gamma_list
    gamma_list = logspace(log10(gamma_min),log10(gamma_max),n_gamma);
    gam_PPCC_list = zeros(size(gamma_list));
    count = 0;
    % first try decreasing gamma:
    for j = n_start:-1:1
        count = count+1;
        % Obtain the Gamma Distribution Parameters for current gamma:
        s_gam_j = stdgaminv(CDF_coarse(:), gamma_list(j)); % standard variate
        mean_s_gam_j = mean(s_gam_j);
        % linear regression:
        beta_coarse_list(j) = (sum(s_gam_j(:).*X_coarse(:))-n_coarse*mean_s_gam_j*mean_X_coarse)/(sum(s_gam_j.^2)-n_coarse*mean_s_gam_j^2);
        mu_coarse_list(j) = mean_X_coarse - beta_coarse_list(j)*mean_s_gam_j;
        % Probability Plot Correlation Coefficient:
        gam_PPCC_list(j) = beta_coarse_list(j)*std(s_gam_j)/std_X_coarse;
        X_coarse_fit_j = mu_coarse_list(j) + beta_coarse_list(j)*s_gam_j;
        if plot_on
            figure(1)
            if count==1, clf; end
            subplot(7,3,count)
            plot(s_gam_j,X_coarse,'.',s_gam_j,X_coarse_fit_j,'-');
            title(['gamma: ' num2str(round(10*gamma_list(j))/10) '   PPCC: ' num2str(gam_PPCC_list(j))]);
        end
        if gam_PPCC_list(j)==max(gam_PPCC_list)
            gam = gamma_list(j);
            gam_PPCC_max = gam_PPCC_list(j);
        else
            break; % stop searching once the PPCC starts to decrease
        end
    end
    if gam_PPCC_list(n_start-1)<gam_PPCC_list(n_start)
        % if the PPCC decreased with decreasing gamma, try increasing gamma: 
        for j = n_start+1:n_gamma
            count = count+1;
            % Obtain the Gamma Distribution Parameters for current gamma:
            s_gam_j = stdgaminv(CDF_coarse(:), gamma_list(j)); % standard variate
            mean_s_gam_j = mean(s_gam_j);
            % linear regression:
            beta_coarse_list(j) = (sum(s_gam_j(:).*X_coarse(:))-n_coarse*mean_s_gam_j*mean_X_coarse)/(sum(s_gam_j.^2)-n_coarse*mean_s_gam_j^2);
            mu_coarse_list(j) = mean_X_coarse - beta_coarse_list(j)*mean_s_gam_j;
            % Probability Plot Correlation Coefficient:
            gam_PPCC_list(j) = beta_coarse_list(j)*std(s_gam_j)/std_X_coarse;
            X_coarse_fit_j = mu_coarse_list(j) + beta_coarse_list(j)*s_gam_j;
            if plot_on
                figure(1)
                subplot(7,3,count)
                plot(s_gam_j,X_coarse,'.',s_gam_j,X_coarse_fit_j,'-');
                title(['gamma: ' num2str(round(10*gamma_list(j))/10) '   PPCC: ' num2str(gam_PPCC_list(j))]);
            end
            if gam_PPCC_list(j)==max(gam_PPCC_list)
                gam = gamma_list(j);
                gam_PPCC_max = gam_PPCC_list(j);
            else
                break; % stop searching once the PPCC starts to decrease
            end
        end
    end
    
    % Obtain the Gamma Distribution Parameters for the best-fit gamma using all the data:
    s_gam = stdgaminv(CDF_X(:), gam); % standard variate
    mean_s_gam = mean(s_gam);
    
    beta = (sum(s_gam(:).*sort_X(:))-n*mean_s_gam*mean_X)/(sum(s_gam.^2)-n*mean_s_gam^2);
    mu = mean_X - beta*mean_s_gam;
    gam_PPCC = beta*std(s_gam)/std_X;

    X_fit = mu + beta*s_gam;
    if plot_on
        figure(1)
        subplot(7,1,6:7)
        plot(gamma_list(find(gam_PPCC_list)),gam_PPCC_list(find(gam_PPCC_list)),'.',gam,gam_PPCC_max,'o')
        xlabel('Shape parameter, gamma');
        ylabel('PPCC');
        set(gcf,'Name','Determination of shape parameter for gamma distribution (press any key to continue...)',...
            'NumberTitle','off');
        pause;
        figure(2)
        plot(s_gam,sort_X,'.',s_gam,X_fit,'-');
        xlabel('Standard variate for gamma distribution');
        ylabel('Value of time series');
        title(['gamma: ' num2str(gam) '    PPCC: ' num2str(gam_PPCC)]);
        legend({'All data','Best-fit gamma distribution'},'Location','NorthWest')
        set(gcf,'Name','Best-fit gamma distribution (press any key to continue...)',...
            'NumberTitle','off');
        pause;
    end

    % Obtain the Normal Distribution Parameters for lower portion of CDF

    CDF_split = .25; % fit normal distribution to data below this cumulative probability
    X_split = interp1(CDF_X, sort_X, CDF_split);
    ind_low = find(sort_X<X_split);
    X_low = sort_X(ind_low);
    n_low = length(X_low);
    CDF_low = CDF_X(ind_low);

    s_norm_low = stdnorminv(CDF_low); % standard normal variate
    mean_s_norm_low = mean(s_norm_low);
    mean_X_low = mean(X_low);
    % linear regression:
    sigma_low = (sum(s_norm_low(:).*X_low(:))-n_low*mean_s_norm_low*mean_X_low)/(sum(s_norm_low.^2)-n_low*mean_s_norm_low^2);
    mu_low = mean_X_low - sigma_low*mean_s_norm_low;
    X_low_fit = mu_low + sigma_low*s_norm_low;
    % Probability Plot Correlation Coefficient:
    norm_PPCC = sigma_low*std(s_norm_low)/std(X_low);

    if plot_on
        figure(3)
        plot(s_norm_low,X_low,'.',s_norm_low,X_low_fit,'-');
        xlabel('Standard normal variate');
        ylabel('Value of time series');
        title({'Normal distribution fit to lower tail',...
            ['Probability Plot Correlation Coefficient: ' num2str(norm_PPCC)]});
        legend({'Lower tail data','Best-fit Normal distribution'},'Location','NorthWest')
        set(gcf,'Name','Fit of normal distribution to lower tail (press any key to continue...)',...
            'NumberTitle','off');
        pause;
        figure(4);
        plot(sort_X,CDF_X,'k.',X_fit,CDF_X,'r-',X_low_fit,CDF_low,'g-');
        xlabel('Value of time series');
        ylabel('Cumulative Probability');
        legend({'Empirical CDF: all data','Gamma distribution fit to all data','Normal distribution fit to lower tail'},...
            'Location','SouthOutside');
        set(gcf,'Name','Cumulative Distribution Function: empirical and fitted (press any key to continue...)',...
            'NumberTitle','off');
        pause;
    end
   
    % Compute the mean zero upcrossing rate of process y(t) with standardized
    % normal probability distribution.
    X_u = mean(sort_X( find(abs(CDF_X - 0.5) == min(abs(CDF_X - 0.5))) )); % upcrossing level
    Nupcross = length(find( X(2:end)>=X_u & X(1:end-1)<X_u )); % number of upcrossings
    if Nupcross<100
        warning(['The number of median upcrossings is low (' num2str(Nupcross) ...
             '). The record may be too short for accurate peak estimation.']);
    end
    y_pk = sqrt(2.0*log(-dur_ratio*Nupcross ./ log(cdf_pk))); % maximum peaks corresponding to specified cumulative probabilities
    CDF_y = stdnormcdf(y_pk);

    % Perform the mapping procedure to compute the CDF of largest peak
    % for X(t) from y(t)
    X_max = stdgaminv(CDF_y,gam)*beta + mu;
    X_min = stdnorminv(1-CDF_y)*sigma_low + mu_low;
    
    % Compute the Mean of the Peaks for process X(t)
    pdf_pk = -y_pk .* cdf_pk .* log(cdf_pk);
    
    if sign(skew_x)>0
        max_est(i) = trapz(y_pk,pdf_pk.*X_max);
        min_est(i) = trapz(y_pk,pdf_pk.*X_min);
        if nargout==4
            max_std(i) = trapz(y_pk,(X_max-max_est(i)).^2.*pdf_pk);
            min_std(i) = trapz(y_pk,(X_min-min_est(i)).^2.*pdf_pk);
        end
    else
        max_est(i) = -trapz(y_pk,pdf_pk.*X_min);
        min_est(i) = -trapz(y_pk,pdf_pk.*X_max);
        if nargout==4
            max_std(i) = trapz(y_pk,(-X_min-max_est(i)).^2.*pdf_pk);
            min_std(i) = trapz(y_pk,(-X_max-min_est(i)).^2.*pdf_pk);
        end
    end
end

varargout(1) = {max_est};
varargout(2) = {min_est};
if nargout==4
   varargout(3) = {max_std};
   varargout(4) = {min_std};
end