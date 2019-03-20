function x = stdgaminv( p, gam )

% This function, written by Joseph A. Main, evaluates inverse values of the 
% cumulative distribution function (CDF) for the standard gamma distribution. 
% The CDF for the standard gamma distribution is equivalent to the incomplete
% gamma function "gammainc" in MATLAB:  p = gammainc(x,gam)
%
% INPUT ARGUMENTS:
% "p" is a vector of cumulative probabilities (probabilities of non-
%     exceedance, with values between 0 and 1 ) for which values of the 
%     random variable are requested
% "gam" is the shape parameter of the standard gamma distribution, and 
%     it must be between 0.1 and 150
%
% OUTPUT ARGUMENT:
% "x" is a vector of values of the random variable corresponding to the 
%     cumulative probabilities specified in "p"
%
%  updated 9 Nov 2006

% set error tolerances for the output x:
abs_tol = 10^-3; % absolute error tolerance
rel_tol = 10^-3; % relative error tolerance

% Range in "gam" is limited by polynomial fits to max and min x values:
if gam<0.1 || gam>150
    error('The shape parameter gamma must be between 0.1 and 150');
end

% Estimate maximum and minimum values of x required for interpolation:
% initial estimate for maximum required x value:
% (4th order polynomial fit to value with exceedance probability of 10^-6)
x_max = 10.^polyval([-0.009486738 0.03376901 0.1151316 0.2358172 1.139717],log10(gam));
% increase x_max iteratively if a higher cumulative probability is required:
max_iter = 200;   % maximum number of iterations
iter = 0;
while gammainc(x_max,gam)<max(p)
    iter = iter+1;
    if iter>max_iter
        error(['Maximum specified probability is too high: ' num2str(max(p))]);
    end
    x_max = 1.5*x_max;
end
% initial estimate for minimum required x value:
% (7th order polynomial fit to value with cumulative probability of 10^-6)
x_min = 10.^polyval([-0.0854665 0.866249 -3.25511 5.14328 -0.90924 -8.09135 12.3393 -5.89628],log10(gam));
% reduce x_min iteratively if a lower cumulative probability is required:
iter = 0;
while gammainc(x_min,gam)>min(p)
    iter = iter+1;
    if iter>max_iter
        error(['Minimum specified probability is too low: ' num2str(min(p))]);
    end
    x_min = 0.1*x_min;
end

% assemble a list of values and probabilities for interpolation:
n_check = 1000;
x_check = linspace(x_min,x_max,n_check);
p_check = gammainc(x_check,gam);

[p_check, ind_u, junk] = unique(p_check); % ensure no repititions
x_check = x_check(ind_u);

% interpolate to estimate values corresponding to specified probabilities:
x_est = interp1(p_check, x_check, p);

% iterate the interpolation until error tolerances are satisfied:
max_iter = 15;   % maximum number of iterations
iter = 0;
x_step = ones(size(x_est)); % initialize vector of step sizes
while any(abs(x_step)>abs_tol) && any(abs(x_step./x_est)>rel_tol)
    iter = iter+1;
    if iter>max_iter
        warning(['Max interations reached. Max abs error: ' num2str(max(abs(x_step))) '; Max rel error: ' num2str(max(abs(x_step./x_est))) ]);
        break;
    end

    % compute probabilities associated with estimated values:
    p_est = gammainc(x_est,gam);
    
    % augment the list of values and probabilities for interpolation
    p_check = [p_check(:); p_est(:)];
    x_check = [x_check(:); x_est(:)];
    [p_check, ind_u, junk] = unique(p_check); % ensure no repetitions
    x_check = x_check(ind_u);

    x_interp = interp1(p_check, x_check, p); % interpolate values
    x_step = x_interp-x_est;                 % evaluate change in interpolated values
    x_est = x_interp;                        % update estimated values
end
x = reshape(x_est,size(p));