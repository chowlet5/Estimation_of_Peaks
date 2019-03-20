function [max_qnt, min_qnt] = maxminqnt(record, dur_ratio, CDF_qnt)

% This function, written by Joseph A. Main, is a modified version of a 
% function previously written by Fahim Sadek, based on the procedure 
% described in the following paper:
% 
% Sadek, F. and Simiu, E. (2002). "Peak non-gaussian wind effects for
% database-assisted low-rise building design." Journal of Engineering
% Mechanics, 128(5), 530-539.
%
% The function estimates probability distributions for the peaks of a
% time series and returns maximum and minimum values corresponding to the
% specified probabilities of non-exceedence.
%
% INPUT ARGUMENTS:
% Each row of "record" is a time series.
% "dur_ratio" = (duration for expected peaks)/(duration of record).
%     If unspecified, a value of 1 is used.
% "CDF_qnt" is a vector giving the probabilities of non-exceedence for which
%     quantiles are desired. For the minima, quantiles corresponding to a
%     cumulative distribution of 1-CDF_qnt are returned in the output argument
%     "min_qnt". For example, if a value of 0.75 is specified for "CDF_qnt", 
%     then the resulting value of "max_qnt" has a 25% probability of being 
%     exceeded in the positive direction, while the resulting value of 
%     "min_qnt" has a 25% probability of being exceeded in the negative direction.
%
% OUTPUT ARGUMENTS:
% "max_qnt" gives the maximum values of each row of "record"
%     corresponding to the non-exceedance probabilities in "CDF_qnt".
% "min_qnt" gives the minimum values of each row of "record"
%     corresponding to the non-exceedance probabilities in "CDF_qnt".
%
% Both "max_qnt" and "min_qnt" are r by q matrices, where r 
%     is the number of rows in "record", and q is the number of 
%     probabilities in "CDF_qnt".
%
% updated 7 Nov 2006

if nargin~=3
    error(['Three input arguments expected (' nargin ' provided)']);
end

if isempty( record )
    error('First input argument (time series) must not be empty');
elseif isempty( CDF_qnt )
    error('Third input argument (probabilities of non-exceedance) must not be empty');
elseif isempty( dur_ratio )
    dur_ratio = 1;
end

plot_on = 1; % turns plotting on (1) [for diagnostics] or off (0)

if max(CDF_qnt)>=1 || min(CDF_qnt)<=0
    error('The range of the second input argument must be between 0 and 1.');
end

rsize = size(record);

CDF_qnt = reshape(CDF_qnt,1,length(CDF_qnt)); % make sure CDF_qnt is a row vector

max_qnt = zeros(rsize(1),length(CDF_qnt));
min_qnt = zeros(rsize(1),length(CDF_qnt));

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
    y_pk = sqrt(2.0*log(-dur_ratio*Nupcross ./ log(CDF_qnt))); % maximum peaks corresponding to specified cumulative probabilities
    CDF_y = stdnormcdf(y_pk);

    % Perform the mapping procedure to compute the CDF of largest peak
    % for X(t) from y(t)
    X_max = stdgaminv(CDF_y,gam)*beta + mu;
    X_min = stdnorminv(1-CDF_y)*sigma_low + mu_low;
    
    if sign(skew_x)>0
        max_qnt(i,:) = X_max;
        min_qnt(i,:) = X_min;
    else
        max_qnt(i,:) = -X_min;
        min_qnt(i,:) = -X_max;
    end

end