function [Y,f,CI,G] = sinFitFT(t,y,varargin)
% Function to calculate a Fourier transform by fitting sines and cosines of
% Fourier components. Each Fourier component has the form
%   a*cos(2*pi*fj*x) + b*sin(2*pi*fj*x)
% where fj is the Fourier frequency. Coefficients are computed for the
% Fourier frequencies from fj=0 up to the Nyquist frequency
%
% [Y,f,CI,G] = sinFitFT(t,y)
%
% Input
%   t: time vector, need not be equally spaced, function expects a column
%       vector
%   y: 1-, 2-, or 3-dimensional time series, time dimension must be first
%       dimension
%
% Output
%   Y: matrix of Fourier coefficients of the form a+ib, has the same
%       dimensions as y except the first dimension is the number of Fourier
%       components calculated
%   f: column vector of Fourier frequencies corresponding to the first
%       dimension of Y
%   CI: confidence intervals of Fourier coefficients, last 2 dimensions are
%       lower and upper confidence bounds, respectively, with real (imag)
%       components corresponding to confidence in a (b)
%   G: goodness of fit statistics, last 2 dimensions are adjusted r-squared
%       values and rmse values, respectively
%
% 13 Dec 2023
% KJW
%
% Notes:
% add ability to specify input frequencies
% add ability to specify fitting method

nan_limit = 0.5;

N = size(y,1);
T = diff(t([1 end]));

% Fourier frequencies up to Nyquist frequency
f = (1/T)*(0:floor(N/2))';
nf = length(f);


% Define fitting routine
foptions = fitoptions('method','nonlinearleastsquares','startpoint',[1 0]);
ftype = fittype('a*cos(2*pi*fj*x) + b*sin(2*pi*fj*x)','problem','fj','options',foptions);

% preallocate
ndim2 = size(y,2);
ndim3 = size(y,3);
Y = nan([nf ndim2 ndim3]);
CI = nan([nf ndim2 ndim3 2]);
G = nan([nf ndim2 ndim3 2]);

% main loop
for i = 1:ndim2
    for j = 1:ndim3
        tij = t;
        yij = y(:,i,j);
        % remove nans
        nan_idx = isnan(yij);
        if sum(nan_idx)/N>nan_limit
            continue
        end
        
        tij(nan_idx) = [];
        yij(nan_idx) = [];

        % fitting loop
        for k = 1:length(f)
            if f(k)==0
                Y(1,i,j) = mean(yij) + 1i*0; % take mean if f==0
            else
                [fijk,gijk] = fit(tij,yij,ftype,'problem',f(k)); % fit sin + cos
                cijk = confint(fijk);
                % coefficients
                Y(k,i,j) = fijk.a + 1i*fijk.b;
                % confidence intervals
                CI(k,i,j,1) = cijk(1,1) + 1i*cijk(1,2); % lower
                CI(k,i,j,2) = cijk(2,1) + 1i*cijk(2,2); % upper
                % goodness stats
                G(k,i,j,1) = gijk.adjrsquare;
                G(k,i,j,2) = gijk.rmse;
            end
        end

    end
end

% get rid of empty dimensions
CI = squeeze(CI);
G = squeeze(G);