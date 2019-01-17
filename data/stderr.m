function stderror = stderr(x);
% calculates standard error for timeseries based on DoF code

stderror = nanstd(x)/sqrt(DoF([1:length(x)],x)-1);

end

function [N] = DoF(dn,timeseries);
%% N = DoF(dn,timeseries);
% Determines degrees of freedom within timeseries measurement
%
% INPUT:
%   dn          - datenum
%   timeseries  - data in question
%
% OUTPUT:
%   N           - degrees of freedom
%

% Corlett, 2017

[~,n]=size(dn);
if n~=1 %Correct orientation of dn
    dn=dn';
end;[m,~]=size(dn);
[m2,~]=size(timeseries);
if m2~=m %Correct orientation of timeseries
    timeseries=timeseries';
end
[m2,datasets] = size(timeseries);

% Regular Spaced Data

delT = nanmedian(diff(dn));

for i=1:datasets

    timeseries(~isnan(timeseries(:,i)),i) = detrend(timeseries(~isnan(timeseries(:,i)),i));
    
    %Recalculate m based on number of data points 
    % to retaining both accuracy and position information
    m = length(find(~isnan(timeseries(:,i))));

    Sxx = NaN(m2-1,1);
    timeseries_mn = nanmean(timeseries(:,i));

    for j=0:m2-2

        dif2 = (timeseries(1:m2-j,i)-timeseries_mn).*(timeseries(1+j:m2,i)-timeseries_mn);
        if length(find(~isnan(dif2)))>1
            Sxx(j+1,1) = nansum(dif2)/(length(find(~isnan(dif2)))-1);
        else
            Sxx(j+1,1) = NaN;
        end    

    end

    dif = (timeseries(:,i) - timeseries_mn).^2;
    Sx2=nansum(dif)/(m-1);

    rxx = Sxx./Sx2;

    tx = nansum(rxx*delT);
        
    N(1,i) = abs((length(find(~isnan(timeseries(:,i))))*delT)/(2*tx));
    
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