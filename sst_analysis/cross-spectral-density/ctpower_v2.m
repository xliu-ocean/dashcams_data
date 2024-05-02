function [S,f,ci,boostCoef]=ctpower_v2(y,dt,t,varargin)
%
% [S,f,ci,boostCoef]=ctpower(y,dt,[options])
%
% computes power spectral density by the Cooley-Tookey method (S =
% abs(fft(y))^2/df).  If y is a 2xN matrix, it is treated as two signals and
% the cross-spectrum is calculated.
%
% OUTPUT: 'S' is the power (cross-)spectral density estimate at frequencies
% 'f'.  'ci' is the "confidence interval multiplier":
%
%     ci(1)*S(f)  <=  true S  <=  ci(2)*S(f),
%      at confidence level ciThresh (default 95%, see below)
%
% if windowing is used, 'boostCoef' is the variance of the windowed
% timeseries.  This can be used to boost the variance of the output
% spectrum (multiply by boostCoef).
%
%
% OPTIONS:
%
% zeroPad=n: freqeuncy resolution in the output will be a factor of n
% increased (if n>1 then it'll be more resolved).  Default zeroPad=1 (--> f
% is the fourier frequencies).
%
% sided={1,2}: 1 (default) outputs the 1-sided S.  2 outputs 2-sided
%
% dof=n: uses non-overlapping band averages to obtain n degrees of freedom.
% Default n=2 (that is, no band averaging)
%
% window={'rectangle','bartlett'}: applies different windows 
%
% ciThresh=x: outputs x-percent (two-sided) confidence interval.  Default
% x=95
%
% whiten=n: applies pre-whitening (and post-coloring) n times.  Default
% n=0.  The whitening will be applied after windowing, although this
% combination of processing would be wierd.



%--------------------------------------------------
% parse optional arguments
%--------------------------------------------------
if(mod(length(varargin),2) == 1)
  error(['problem parsing optional arguments (use format "name,value,' ...
         'name,value..."'])
end
for i = 1:2:length(varargin)
  switch varargin{i};
   case 'zeropad'
    zeropad=varargin{i+1};
   case 'sided'
    sided=varargin{i+1};
   case 'dof'
    dof=varargin{i+1};
   case 'window'
    win=varargin{i+1};
   case 'ciThresh'
    ciThresh=varargin{i+1};
   case 'whiten' 
    whiten=varargin{i+1};
   otherwise
    error(['invalid option "' varargin{i} '"'])
  end
end
% defaults...
if(~exist('zeroPad')) zeroPad=1; end
if(~exist('sided')) sided=1; end
if(~exist('dof')) dof=2; end
if(~exist('win')) win='rectangle'; end
if(~exist('ciThresh')) ciThresh=95; end
if(~exist('whiten')) whiten=0; end
% checks for invalid input
if(mod(dof,2) | ~mod(dof/2,2))
  error('dof must be 2*m where m is odd');
end
if(~(strcmp(win,'bartlett') | strcmp(win,'rectangle')))
  error(['invalid window type: "' win '".  Choices are rectangle or ' ...
         'bartlett']);
end
if(sided~=1 & sided~=2)
  error(['invalid input for sided: must be 1 or 2']);
end
if(dof>2 & sided==2)
  error('band-averaging not implemented for 2-sided S estimates')
end

% treat PSD as a special case of cross-spectral density
if(isvector(y))
  y = [y y];
  vecFlag=1;
end
if(whiten>0 & strcmp(win,'bartlett'))
  warning(['whitening and windowing accomplish the same task.  '...
           'They aren''t normally used together.'])
end

%--------------------------------------------------
% begin data processing
%--------------------------------------------------

N=size(y,1);
Np=N*zeroPad;

% apply window
varOrig = mean(y(:,1).^2);
if(strcmp(win,'rectangle'))
  ;
elseif(strcmp(win,'bartlett'))
  y(:,1)=y(:,1)-mean(y(:,1));
  y(:,2)=y(:,2)-mean(y(:,2));
  y(:,1)=myBartlett(y(:,1));
  y(:,2)=myBartlett(y(:,2));
end
varWind = mean(y(:,1).^2);
boostCoef = varOrig/varWind;

% pre-whitening
for i=1:whiten
  disp('whiten')
  y=diff(y);
  y = [[y(1,1);y(:,1)] [y(1,2);y(:,2)]];
end

% zero padding
warning off  % complains about concatenation with empty array
y=[y;zeros(Np-N,2)];
warning on

% do the fft
% use sinFitFT
% t0 = ((1:N)*dt)';
X=sinFitFT(t,y(:,1));
Y=sinFitFT(t,y(:,2));
% X=fft(y(:,1))/N;
% Y=fft(y(:,2))/N;
S = N*dt*conj(X).*Y;

% special case of vector input means we should expect real output.  But
% round-off error often makes it slightly complex.
if(exist('vecFlag'))
  S=real(S);
end

% set up the indeces to pull out of the fft, and the frequencies that
% they correspond to.  'chopInd' is the first one to wrap around to
% "negative" frequencies, when indexing the data from zero.
chopInd=ceil(Np/2);
if(sided==2)
  ind=[chopInd:(Np-1) 0:(chopInd-1)]+1;
  S = S(ind);
  f=[[[chopInd:(Np-1)]-Np] [0:(chopInd-1)]]/(Np*dt);
elseif(sided==1)
  ind=[1:(chopInd-1)]+1;  % note: omit mean-value
  S = 2*S(ind);
  f=[1:(chopInd-1)]/(Np*dt);
end

% post-coloring
GG=4*sin(pi*f*dt).^2;
for i=1:whiten
  disp('color')
  S=S./GG';
end

% band-averaging to increase dof
[f,S] = bandAverage(f,S,dof);

% confidence interval "multipliers"
alpha=1-ciThresh/100;
ci=dof./chi2inv([1-alpha/2 alpha/2],dof);
