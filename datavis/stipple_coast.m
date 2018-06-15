% [Px,Py] = stipple_coast(lat_overall, lon_overall [, latrng, lonrng, N, lim ] )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stippled coastline using weighted Lloyd's algorithm.
% The function receives input of the coast outline and visually iterates through stippling.
% The function produces coordinates of the stipple points.
%
% Inputs:
%   * lat_overall, lon_overall: regional coastline vectors (outline of coast)
%   * latrng, lonrng: (optional) window of interest for stippling (default: range of input vectors)
%   * N: (optional) number of dots in the stipple image (default: 20000)
%   * lim: (optional) buffer distance inland of coastline (default: 0.02 deg)
%
% Outputs:
%   * Px, Py: coordinates of stipple points
%
% Use:
% 	figure;
%   [Px,Py]=stipple_coast(lat,lon,[40.6, 40.75],[-74.19, -74.05], 20000, 0.02);
%   hold on;plot(lon,lat,'k');
%
% Requirements:
% 	* Matlab mapping toolbox (to determine land polygons)
% 	* OPTIONAL - Matlab parallel computing toolbox (to speed up iterations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. Corlett, 2018-05-15

% edited from https://github.com/mfkasim91/stippling-lloyds.git (updated 2016)
% originally based on A. Secord (2002). Weighted Voronoi Stippling. NPAR,
% Annecy, France. https://www.cs.ubc.ca/labs/imager/tr/2002/secord2002b/secord.2002b.pdf

function [Px,Py] = stipple_coast(lat_overall, lon_overall, latrng, lonrng, N, lim)
    if (nargin < 3)
        % set stippling bounds
        lonrng = [nanmin(lon_overall(:)), nanmax(lon_overall(:))];
        latrng = [nanmin(lat_overall(:)), nanmax(lat_overall(:))];
        
        % set point limit + range
        N = 20000;
        lim = 0.02;
    elseif (nargin < 5)
         % set point limit + range
        N = 20000;
        lim = 0.02;
    end
    
    % deploy initial random samples using a simple rejection method
    [Px0,Py0,map] = initial_coast(lat_overall, lon_overall, N, lim, ...
        latrng, lonrng);
        
    % upsample data for higher resolution
    [lon_overall,lat_overall] = smoothsample(lon_overall,lat_overall);
%     lon_overall = upsample(lon_overall,10);

    % find data within window of interest
    good = (lon_overall<max(lonrng) & lon_overall>min(lonrng) & lat_overall<max(latrng) & lat_overall>min(latrng));
    lat = lat_overall(good);
    lon = lon_overall(good);
    
    % use lloyds' algorithm
    options = struct();
    options.maxIter = 35;
    options.verbose = 1;
    [Px, Py] = wtd_lloyds(Px0, Py0, map, lim, lat, lon, latrng, lonrng, options);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating N random sample with probability distribution valMap using rejection algorithm.
% Input:
%   * N: number of sample
%   * valMap: discrete probability function
% Output:
%   * Px, Py: normalised coordinate of the generated sample (Nx1 each)

function [Px,Py,map] = initial_coast(lat_overall, lon_overall, N, lim, latrng, lonrng)

    if (nargin < 3)
        N = 25000;
        lim = 0.01; % max distance from line
        lonrng = [-74.2, -74.05];
        latrng = [40.55, 40.80];
    end
    
    map = intersect(polyshape([min(lonrng),max(latrng);max(lonrng),max(latrng);max(lonrng),min(latrng);min(lonrng),min(latrng);min(lonrng),max(latrng)]),polyshape([lon_overall,lat_overall]));

    % upsample data for higher resolution
    [lon_overall,lat_overall] = smoothsample(lon_overall,lat_overall);
%     lon_overall = smoothsample(lon_overall,10);
    
    % find data within window of interest
    good = (lon_overall<max(lonrng) & lon_overall>min(lonrng) & lat_overall<max(latrng) & lat_overall>min(latrng));
    
    lat = lat_overall(good);
    lon = lon_overall(good);
    
    Px = nan([N,1]);
    Py = nan([N,1]);
    reverseStr='';
    while (any(isnan(Px)))
        nanIdx = find(isnan(Px)); % unfilled position index
        nanN = length(nanIdx); % number of unfilled position
        
        % generate coordinates
        x = rand([nanN,1])*((nanmax(lon)-nanmin(lon))-1e-3)+nanmin(lon)+5e-4;
        y = rand([nanN,1])*((nanmax(lat)-nanmin(lat))-1e-3)+nanmin(lat)+5e-4;
        
        % check the range
        inside = (x >= nanmin(lon)) & (x <= nanmax(lon)) & (y >= nanmin(lat)) & (y <= nanmax(lat));
        ix = x .* inside + (1-inside);
        iy = y .* inside + (1-inside);
        % disp(any(Nx < ix));
        
        % determine the acceptance of the coordinates
        % calculate minimum distance to coast
        prob = zeros(size(ix));
        dist = ones(size(ix))*lim;
        for i = 1:length(inside)
            dist(i,1) = nanmin(sqrt( (lon - ix(i)).^2 + (lat - iy(i)).^2)); %calculate minimum distance
        end
        %calculate PDF on the fly
        prob(dist<lim) = (-log((dist(dist<lim)./lim) + 0.01)+log(1.01))./(-log(0.01)+log(1.01)); %calculate logarithmic
        prob(dist>=lim) = 0; % truncate errors from curve
%         prob(dist<lim) = ((dist(dist<lim)-lim).^4) ./ (lim.^4); %calculate quartic 
%         prob(dist>=lim) = 0; % truncate errors from curve
%         prob = (-1/lim).*dist + 1; %calculate linear 
%         prob(prob<0) = 0; %truncate errors from curve
        
        %remove values between land
        [in,~] = inpolygon(ix,iy,map.Vertices(:,1),map.Vertices(:,2));
%         in = inpoly([ix,iy],poly.Vertices);
        bad = zeros([nanN,1]);
        bad(in) = 1;
        prob = prob.*bad;
        
        r = rand([nanN,1]);
        accept = (r < prob) & inside; %find accepted points
        
        % fill the accepted coordinates
        filledIdx = nanIdx(accept);
        Px(filledIdx) = x(accept);
        Py(filledIdx) = y(accept);
        
        % Display progress
        percentDone = 100 * (double(length(find(isnan(Px)))) / double(N));
        msg = sprintf('Finding scatter points - missing: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
    
end

%% upsample via linear interpolation

function [xx]=upsample(x,n)
% increases sampling interval by 'n'
% using linear interpolation

[~,q]=size(x);
if q~=1 %Correct orientation of x
    x=x';
end

% linear equation
m = x(2:end) - x(1:end-1);
b = x(1:end-1);

% interpolation
x_new = [1:n]/(n+1);
mid = (m*x_new) + repmat(b,1,n);
xx = [x(1:end-1,1),mid]';
xx = [xx(:);x(end)];

end

function [xx,yy]=smoothsample(x,y)
% increases some sampling intervals to median interval
% using linear interpolation

[~,q]=size(x);
if q~=1 %Correct orientation of x
    x=x';
end
[~,q]=size(y);
if q~=1 %Correct orientation of x
    y=y';
end

% shift x coordinates
x_shift = nanmax(x);
x = x + x_shift;

% linear equation
mx = x(2:end) - x(1:end-1);
bx = x(1:end-1);

% linear equation
my = y(2:end) - y(1:end-1);
by = y(1:end-1);

% median difference
mid_dif = nanmedian(abs(diff(x)));
n = floor(abs(mx./mid_dif));

% save memory using sparse matrices to store variable length interpolations
Qx = sparse(length(mx),nanmax(n)+1);
Qx(:,1) = x(1:end-1,1);
Qy = sparse(length(my),nanmax(n)+1);
Qy(:,1) = y(1:end-1,1);

% interpolation
bad = find(n>1);
if ~isempty(bad)
    for i = 1:length(bad)
        x_new = [1:n(bad(i))]/(n(bad(i))+1); y_new = x_new;
        Qx(bad(i),2:n(bad(i))+1) = (mx(bad(i))*x_new) + repmat(bx(bad(i)),1,n(bad(i)));
        Qy(bad(i),2:n(bad(i))+1) = (my(bad(i))*y_new) + repmat(by(bad(i)),1,n(bad(i)));
    end
end

% reshape
Qx = Qx'; Qy = Qy';
xx = [nonzeros(Qx);x(end)]-x_shift;
yy = [nonzeros(Qy);y(end)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying Lloyd's algorithm for weighted area in parallel.
% Input:
%   * Px0, Py0: list of initial points coordinates (each numPoints x 1)
%   * valMap: value of density per pixel (Ny x Nx)
%   * options:
%       * maxIter: maximum number of iterations
% Output:
%   * Px, Py: list of final points coordinates after applying the Lloyd's algorithm (each numPoints x 1)
%   * Ap: area for each voronoi cell formed by the points (numPoints x 1)

function [Px, Py] = wtd_lloyds(Px0, Py0, map, lim, lat, lon, latrng, lonrng, options)
    

    crs = [min(lonrng),max(latrng);...
        max(lonrng),max(latrng);...
        max(lonrng),min(latrng);...
        min(lonrng),min(latrng);...
        min(lonrng),max(latrng)];
         
    if (nargin < 4)
        options = struct();
    end
    if ~isfield(options, 'maxIter'); options.maxIter = 20; end;
    if ~isfield(options, 'verbose'); options.verbose = 0; end;
    
    maxIter = options.maxIter;
    
    Px = Px0;
    Py = Py0;
    
    if (options.verbose)
        figure; pause(0.1);
    end
    
    for (i = [1:maxIter])
        disp(i);

        % remove data out of range
        bad = find(isnan(Px) | isnan(Py) | isinf(Py) | isinf(Py));
        Px(bad) = []; Py(bad) = [];
    
        [V, C] = power_bounded(Px, Py, crs);
        
%         reverseStr='';
        parfor (j = [1:length(C)])
            % Display progress
%             percentDone = 100 * (double(j) / double(length(C)));
%             msg = sprintf('Finding scatter points - missing: %3.1f', percentDone); %Don't forget this semicolon
%             fprintf([reverseStr, msg]);
%             reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
            % get the centroid
            [xp, yp] = poly2cw(V(C{j},1), V(C{j},2));
            
            % remove values between land
            [in] = inpolygon(xp,yp,map.Vertices(:,1),map.Vertices(:,2));
            bad = [1:length(xp)]';
            bad(in) = [];
            
            % shift water bounds onto coastline
            if ~isempty(bad)
                Q = intersect(polyshape([xp,yp]),map);
                xp = Q.Vertices(:,1); yp = Q.Vertices(:,2);
                dist = zeros(size(xp));
                for K = 1:length(xp)
                    dist(K,1) =  nanmin(sqrt( (lon - xp(K)).^2 + (lat - yp(K)).^2)); %calculate minimum distance for weighting
                end
            else
                dist=zeros(size(xp));
                for K = 1:length(xp)
                    dist(K,1) = nanmin(sqrt( (lon - xp(K)).^2 + (lat - yp(K)).^2)); %calculate minimum distance for weighting
    %                 dist(:,K) = (sqrt( (lon - xp(K)).^2 + (lat - yp(K)).^2)); %calculate minimum distance for weighting
                end
            end
                        
            % create distance vectors
%             distGOOD = nanmean(dist(:,good),2);
%             dist = nanmin(dist,[],1)';
            
            %calculate PDF on the fly
            prob = zeros(size(dist));
            prob(dist<lim) = (-log((dist(dist<lim)./(1.5*lim)) + 0.005)+log(1.005))./(-log(0.005)+log(1.005)); %calculate logarithmic
            prob(dist>=lim & dist<=2.5*lim) = ((-log(2015/3000)+log(1.005))./(-log(0.005)+log(1.005))) * ((-1/(1.5*lim)).*(dist(dist>=lim & dist<=2.5*lim)-lim) + 1); % truncate errors from curve - allow spread onto land
%             prob(dist>=lim & dist<3*lim) = ((-log(2015/3000)+log(1.005))./(-log(0.005)+log(1.005))) * ((-1/(2*lim)).*(dist(dist>=lim & dist<=3*lim)-lim) + 1); % truncate errors from curve - allow spread onto land
% %             prob(dist<lim) = (-log((dist(dist<lim)./lim) + 0.005)+log(1.005))./(-log(0.005)+log(1.005)); %calculate logarithmic
% % %             prob(dist<lim) = (-log((dist(dist<lim)./lim) + 0.02)+log(1.02))./(-log(0.02)+log(1.02)); %calculate logarithmic
%             prob(dist<lim & prob<5e-2) = 5e-2; % truncate errors from curve - allow spread onto land
%             prob(dist>=lim & dist<=3*lim) = 5e-2 * ((-1/(2*lim)).*(dist(dist>=lim & dist<=3*lim)-lim) + 1); % truncate errors from curve - allow spread onto land

%             prob(dist>=2*lim) = 1e-3; % truncate errors from curve - allow spread onto land
%             prob(dist<lim) = ((dist(dist<lim)-lim).^8) ./ (lim.^8); %calculate quartic 
%             prob(dist>=lim) = 0; % truncate errors from curve
    %             prob = (-1/lim).*dist + 1; %calculate linear 
    %             prob(prob<0) = 0; %truncate errors from curve


% %             [in,on] = inpoly([xp,yp],map.Vertices);
%             if length(in)<length(xp)
%                 int = intersect(polyshape([map.Vertices(:,1),map.Vertices(:,2)]),polyshape([xp,yp]));
%                 xp = int.Vertices(:,1); yp = int.Vertices(:,2);
%                 dist = zeros([length(xp),1]);
%                 for K = 1:length(xp)
%                     dist(K,1) = nanmin(sqrt( (lon - xp(K)).^2 + (lat - yp(K)).^2)); %calculate minimum distance for weighting
%                 end
%                 prob = (-1/lim).*dist + 1; %calculate linear 
%                 prob(prob<0) = 0; %truncate errors from curve
%             end
%             bad = zeros([length(xp),1]);
%             bad(in) = 1;
%             prob = prob.*bad;
            
            xc = nansum(xp.*prob)./nansum(prob);
            yc = nansum(yp.*prob)./nansum(prob);
            
            % update positions of center of mass
            Px(j) = xc;
            Py(j) = yc;
            
        end
%         fprintf('\n');
        
        if (options.verbose)
            % close all;
            hold off;
            plot(Px, Py, 'k.', 'markers', 1);
            pause(0.1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POWER_BOUNDED computes the power cells about the points (x,y) inside
% the bounding box (must be a rectangle or a square) crs.  If crs is not supplied, an
% axis-aligned box containing (x,y) is used.
% It is optimised to work fast on large number of sites (e.g. 10000 sites or more)
% Input:
%   * x, y: coordinate of the Voronoi point (numPoints x 1)
%   * wts: weights of each point (numPoints x 1)
%   * crs: vortices of the bounding box in cw order (numVert x 2)
% Output:
%   * V: x,y-coordinate of vertices of the power cells
%   * C: indices of the Voronoi cells from V
% See Matlab's voronoin for more information about the output
% Made by: Aaron Becker, atbecker@uh.edu, and Muhammad Kasim, muhammad.kasim@wolfson.ox.ac.uk

function [V,C] = power_bounded(x, y, crs)
    bnd=[min(x) max(x) min(y) max(y)]; %data bounds
    if nargin < 3
        crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]);
    end

    rgx = max(crs(:,1))-min(crs(:,1));
    rgy = max(crs(:,2))-min(crs(:,2));
    rg = max(rgx,rgy);
    midx = (max(crs(:,1))+min(crs(:,1)))/2;
    midy = (max(crs(:,2))+min(crs(:,2)))/2;

    % add 4 additional edges
    xA = [x; midx + [0;0;-0.75*rgx;+0.75*rgx]];
    yA = [y; midy + [-0.75*rgy;+0.75*rgy;0;0]];
    
    % remove duplicate points
    [~,ia,~] = unique([xA,yA],'rows');
    ia = sort(ia);
    
    % calculate voronoi polygons
    [vi,ci] = voronoin([xA(ia),yA(ia)]);
    
    % remove the last 4 cells
    C = ci(1:end-4);
    V = vi;
    % use Polybool to crop the cells
    %Polybool for restriction of polygons to domain.
    
    maxX = max(crs(:,1)); minX = min(crs(:,1));
    maxY = max(crs(:,2)); minY = min(crs(:,2));
    for ij=1:length(C)
        % thanks to http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
        Cij = C{ij};
        if (isempty(Cij)) continue; end;
        
        % first convert the contour coordinate to clockwise order:
        pts = V(Cij,:);
        if ispolycw(pts(:,1),pts(:,2))==0
%             C{ij} = flipud(Cij(ij));
            X2 = flipud(pts(:,1));
            Y2 = flipud(pts(:,2));
        else
            X2 = pts(:,1);
            Y2 = pts(:,2);
        end
        
        % second check for infinite points
        if (isinf(X2))
            X2(isinf(X2)) = midx + 2*rgx;
        end
        if (isinf(Y2))
            Y2(isinf(Y2)) = midy + 2*rgy;
        end
        
        bad = find(isnan(X2) | isnan(Y2));
        X2(bad) = 0; Y2(bad) = 0;
        
%         K = convhull(pts);
%         K = K(end-1:-1:1);
%         C{ij} = Cij(K);
%         X2 = pts(K,1);
%         Y2 = pts(K,2);
        
        % if all points are inside the bounding box, then skip it
        if (all((X2 <= maxX) & (X2 >= minX) & (Y2 <= maxY) & (Y2 >= minY))) continue; end;
        
        % if no points are not nans or inf, skip it
        if isempty(find(~isinf(X2) & ~isinf(Y2) & ~isnan(X2) & ~isnan(Y2),1)); continue; end;

%         [xb, yb] = clip_polygons(crs(:,1),crs(:,2),X2,Y2);
        int_shp = intersect(polyshape(crs),polyshape([X2(~isinf(X2) & ~isinf(Y2)),Y2(~isinf(X2) & ~isinf(Y2))]));
        xb = int_shp.Vertices(:,1); yb = int_shp.Vertices(:,2);
        
        ix=nan(1,length(xb));
        for il=1:length(xb)
            if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                ix1=find(V(:,1)==xb(il));
                ix2=find(V(:,2)==yb(il));
                for ib=1:length(ix1)
                    if any(ix1(ib)==ix2)
                        ix(il)=ix1(ib);
                    end
                end
                if isnan(ix(il))==1
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    ix(il)=lv+1;
                end
            else
                lv=length(V);
                V(lv+1,1)=xb(il);
                V(lv+1,2)=yb(il);
                ix(il)=lv+1;
            end
        end
        C{ij} = ix;
    end
end
