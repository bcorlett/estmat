function [varargout] = DoF(timeseries);
%% [N_eff(,N_eff_BH)] = DoF(timeseries);
% Determines degrees of freedom for mean within timeseries measurement
% Follows method described by Ken Brink; contains alternative option by 
%
% INPUT:
%   timeseries  - data in question (requires regular spacing in time, though allows NaNs)
%
% OUTPUT:
%   N_eff       - degrees of freedom
%

% B. Corlett, 2019

[m2,n2]=size(timeseries);
if m2<n2 %Correct orientation of timeseries
    timeseries=timeseries';
end
[m2,datasets] = size(timeseries);

% Assume Regularly Spaced Data

for i=1:datasets
    x = timeseries(:,i); % isolate timeseries
    x(~isnan(x),1) = detrend(x(~isnan(x),1)); % detrend + remove mean
    
    Sxx = NaN(m2+(m2-1),1);
    struc = NaN(m2+(m2-1),1);
    xbar = nanmean(x);
    
    % calculate variance
    N = length(find(~isnan(x)));
    Sx2 = nansum( (x - xbar).^2 ) / (N - 1);
    
    % calculate autocovariance with lags from -inf to +inf
    lags = [-m2+1:m2-1];
    for j = lags
        if j < 0
            Sxx(j+m2,1) = nansum( (x(1:(m2+j),1)-xbar).*(x((1-j):m2,1)-xbar) ) / (length(~isnan(x(1:(m2+j)))) - 1);
            struc(j+m2,1) = nansum( (x(1:(m2+j),1)-x((1-j):m2,1)).^2 ) / (length(~isnan(x(1:(m2+j)))));
        else
            Sxx(j+m2,1) = nansum( (x(1:(m2-j),1)-xbar).*(x((1+j):m2,1)-xbar) ) / (length(~isnan(x(1:(m2-j)))) - 1);
            struc(j+m2,1) = nansum( (x(1:(m2-j),1)-x((1+j):m2,1)).^2 ) / (length(~isnan(x(1:(m2-j)))));
        end
    end
    Sxx(isinf(Sxx))=NaN;
    rxx = Sxx./Sx2;
    
    % calculate from Emery & Thompson, 2nd ed. - used in K. Brink's course
    N_eff(i,1) = N / (2 * nansum(rxx(lags>=0)));
    
    if nargout>1
        % Method of Bayley and Hammersley (1946)
        Sxx = NaN(m2-1,1);
        rho_j = NaN(m2-1,1);
        for j=0:m2-1
            n = length(find(~isnan( x(1:(m2-j),1) )));
            Sxx(j+1,1) = nansum( (x(1:(m2-j),1)-xbar) .* (x((1+j):m2,1)-xbar) )/(n-1);
            rho_j(j+1,1) = (n) * (Sxx(j+1,1) / Sx2);
        end
        rho_j(isinf(rho_j)) = NaN;

        N_eff_BH(i,1) = ((1/N) + ((2/(N^2))*nansum(rho_j)))^(-1);
    end
end
      
nout = max(nargout,1);
for k = 1:nout
    if k == 1
        varargout{k} = N_eff;
    elseif k == 2
        varargout{k} = N_eff_BH;
    else
        error('too many outputs; see help file.');
    end
end

end

function unew=detrend(u,norder);
%function unew=detrend(u,[norder]);
%  if u has size u(n,m), operates on n vectors of length m
%   norder is an optional param for 2nd order or 3rd order polynomial 
%   fits to detrend with
if nargin==1; norder=1;end;
[n,m]=size(u);
flip=0;
if m==1;flip=1;u=u.'; [n,m]=size(u); end;
for i=1:n;
utrend=polyfit2(1:m,u(i,:),norder,1:m);
unew(i,:)=u(i,:)-utrend;
end;

if flip; unew=unew'; end;
end