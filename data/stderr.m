function [varagout] = stderr(x);
%% [stderror(,nonparametric_stderror)] = stderr(x);
% calculates standard error about mean for timeseries based on Degree-of-Freedom code

% Corlett, 2017

stderror = nanstd(x)/sqrt(DoF(x)-1);

nout = max(nargout,1);
varargout{1} = stderror;
if nout==2
    [PDF,xPDF]=make_mn_pdf(x);
    nonparametric_stderr = mean(abs([interp1(PDF,xPDF,0.5+0.341)-interp1(PDF,xPDF,0.5),interp1(PDF,xPDF,0.5-0.341)-interp1(PDF,xPDF,0.5)]));
    varargout{2} = nonparametric_stderr;
else 
    error('too many outputs; see help file');
end

end

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

% 
% function [N] = DoF(dn,timeseries);
% %% N = DoF(dn,timeseries);
% % Determines degrees of freedom within timeseries measurement
% %
% % INPUT:
% %   dn          - datenum
% %   timeseries  - data in question
% %
% % OUTPUT:
% %   N           - degrees of freedom
% %
% 
% % Corlett, 2017
% 
% [~,n]=size(dn);
% if n~=1 %Correct orientation of dn
%     dn=dn';
% end;[m,~]=size(dn);
% [m2,~]=size(timeseries);
% if m2~=m %Correct orientation of timeseries
%     timeseries=timeseries';
% end
% [m2,datasets] = size(timeseries);
% 
% % Regular Spaced Data
% 
% delT = nanmedian(diff(dn));
% 
% for i=1:datasets
% 
%     timeseries(~isnan(timeseries(:,i)),i) = detrend(timeseries(~isnan(timeseries(:,i)),i));
%     
%     %Recalculate m based on number of data points 
%     % to retaining both accuracy and position information
%     m = length(find(~isnan(timeseries(:,i))));
% 
%     Sxx = NaN(m2-1,1);
%     timeseries_mn = nanmean(timeseries(:,i));
% 
%     for j=0:m2-2
% 
%         dif2 = (timeseries(1:m2-j,i)-timeseries_mn).*(timeseries(1+j:m2,i)-timeseries_mn);
%         if length(find(~isnan(dif2)))>1
%             Sxx(j+1,1) = nansum(dif2)/(length(find(~isnan(dif2)))-1);
%         else
%             Sxx(j+1,1) = NaN;
%         end    
% 
%     end
% 
%     dif = (timeseries(:,i) - timeseries_mn).^2;
%     Sx2=nansum(dif)/(m-1);
% 
%     rxx = Sxx./Sx2;
% 
%     tx = nansum(rxx*delT);
%         
%     N(1,i) = abs((length(find(~isnan(timeseries(:,i))))*delT)/(2*tx));
%     
% end
%         
% end

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

function y = nanstd(varargin)
%NANSTD Standard deviation, ignoring NaNs.
%   Y = NANSTD(X) returns the sample standard deviation of the values in X,
%   treating NaNs as missing values.  For a vector input, Y is the standard
%   deviation of the non-NaN elements of X.  For a matrix input, Y is a row
%   vector containing the standard deviation of the non-NaN elements in
%   each column of X. For N-D arrays, NANSTD operates along the first
%   non-singleton dimension of X.
%
%   NANSTD normalizes Y by (N-1), where N is the sample size.  This is the
%   square root of an unbiased estimator of the variance of the population
%   from which X is drawn, as long as X consists of independent, identically
%   distributed samples and data are missing at random.
%
%   Y = NANSTD(X,1) normalizes by N and produces the square root of the
%   second moment of the sample about its mean.  NANSTD(X,0) is the same as
%   NANSTD(X).
%
%   Y = NANSTD(X,FLAG,DIM) takes the standard deviation along dimension
%   DIM of X.
%
%   See also STD, NANVAR, NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2006 The MathWorks, Inc.


% Call nanvar(x,flag,dim) with as many inputs as needed
y = sqrt(nanvar(varargin{:}));
end

function [PDF,xPDF]=make_mn_pdf(x);
%% [PDF,xPDF]=make_mn_pdf(x)
% Make PDF of mean from discrete data using bootstrapping
%
% INPUT:
%   x - dataset
%
% OUTPUT:
%
%   PDF - percentage
%  xPDF - observed value
%

% Corlett, 2017

[~,MN,~]=bootstrap(x(~isnan(x)),length(find(~isnan(x))),1000);
[n,xPDF] = hist(MN(:),1000);
PDF=cumsum(n/length(MN(:)));
[PDF,ia,~]=unique(PDF);
xPDF=xPDF(ia);

end

function [Series,Average,Std_Dev]=bootstrap(input,points,iterations);
%% [Series,Average,Std_Dev]=bootstrap(input,points,iterations);
%   Bootstrap method for creating a large composite of "samples"
%
%   INPUT:
%   input = timeseries to perform method on
%   points = number of points to be replaced
%   iterations = number of times to recreate timeseries
%
%   OUTPUT:
%   Series = string of multiple timeseries
%   Average = mean over bootstrapped dataset
%   Std_Dev = standard deviation of bootstrapped dataset
%
% Corlett, 2013.

%Ensure data is a column vector
[a,~]=size(input);
if a==1;
    input=input';
end

%Remove NaN's from dataset
input(isnan(input))=[];

%Cap maximum number of replacement points
if points>length(input)
    fprintf('The input for "points" has been capped at %d\n', length(input));
    points=length(input);
end

%Replace and append to output
for M=1:iterations
    
%     i=randsample(length(input),points);
        i=1+floor(rand(points,1)*length(input));
%     j=randsample(length(input),points);
        j=1+floor(rand(points,1)*length(input));

    NEW=input;
    NEW(i,1)=input(j,1);
    
    Series(:,M)=NEW;
    clear NEW i j;
end

    Average(:,1)=mean(Series,1);
    Std_Dev(:,1)=std(Series,0,1);
end

