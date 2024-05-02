function [f,gsq,phi,gsq_crit,del_phi]=crossSpec_v2(y1,y2,dt,t,dof,window)
%
% [f,gsq,phi,gsq_crit,del_phi]=crossSpec(y1,y2,dt,[dof,window])
%
% Cross-spectral analysis for coherence and phase between two timeseries
%
% INPUTS:
%
% y1,y2: two timeseries (dims 1xN), sampled at same rate
% dt: rate of sampling
% dof (optional): apply band-averaging degrees of freedom.  Note dof must be
%                 equal to 2*n, where n is an odd number.  Default: no
%                 band-averaging, dof==2.
% window (optional): if window==1, apply Bartlett windowing.  Default: no
%                    windowing
%
% OUTPUTS:
%
% f: frequencies for output vectors
% gsq: coherence (gamma^2)
% phi: phase
% gsq_crit: 95% significance level for gsq.  That is, if gsq exceeds
%           gsq_crit, then the coherence is significant at the 95%
%           confidence level
% del_phi: 95% confidence interval for phi.  That is, the true phi lies
%          within phi +/- del_phi, at the 95% confidence level.
%

% set default options
if(exist('dof')~=1)
  dof=2;
end
if(exist('window')~=1)
  window=0;
end

% confirm that inputs have proper dimensions
if(isempty(find(size(y1,1))) | isempty(find(size(y2,1))))
  error('inputs must be vectors, not matrices')
end
y1=y1(:);
y2=y2(:);

% compute spectra and cross-spectrum
if(window==1)
  [S1,f,ci]=ctpower_v2(y1,dt,t,'window','bartlett','dof',dof);
  [S2,f,ci]=ctpower_v2(y2,dt,t,'window','bartlett','dof',dof);
  [S12,f,ci]=ctpower_v2([y1 y2],dt,t,'window','bartlett','dof',dof);
else
  [S1,f,ci]=ctpower_v2(y1,dt,t,'dof',dof);
  [S2,f,ci]=ctpower_v2(y2,dt,t,'dof',dof);
  [S12,f,ci]=ctpower_v2([y1 y2],dt,t,'dof',dof);
end

% \gamma_{xy}^2
gsq = abs(S12).^2./(S1.*S2);
% \phi_{xy}
phi = angle(S12);
% 95% significance level for gsq_xy
alpha = .95;
gsq_crit = finv(alpha,2,dof-2)/(.5*(dof-2) + finv(alpha,2,dof-2));
% CI for \phi_{xy}
del_phi = asin(sqrt((2/(dof-2))*((1-gsq)./gsq)*finv(alpha,2,dof-2)));
