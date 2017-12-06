%Generate 10,000 random numbers using the linear congruential generator 
format long;
modulus = 4294967296; aa =22695477 ; const =1; seed = 66754;
snext = seed;
rnds(1) = seed;
for i = 1:10000, rnds(i+1) = mod(rnds(i)*aa+const,modulus); end
rnds = rnds/modulus;
%1)	Compute the mean value of the random numbers that you generated,
%and compare with the mean value of a uniform distribution between 0 and 1.
uni = rand([1000 1]);
mean_rnds = mean(rnds);
mean_uni = mean(uni);
%2)	Compute the mean value of the random numbers that you generated,
%and compare with the mean value of a uniform distribution between 0 and 1.
std_rnds = std(rnds);
std_uni = std(uni);
%3)	Compute the first 500 points of the autocorrelation function,
%and plot the results. Do the numbers appear to be uncorrelated?
auto_rnds = zeros(500,1);
for i = 1:500
  auto_rnds(i) = (rnds(i)-mean_rnds)* (rnds(i+1)-mean_rnds); 
end
figure;
plot(auto_rnds);
title('Autocorrleation for the generated random numbers');
saveas(gcf,'autocorr.png')
%	Using the rejection sampling method seen in class, and the above random number generator 
%defined above, generate 5,000 samples
%following an exponential distribution ?e^(-?x) with  ?=1.5. Plot the histogram of the samples.
lmda = 1.5 ;
xmin = 0; %the lower limit of our domain
xmax = 5; %the upper limit of our domain
ymax = lmda;
rnds_from_exp = zeros(5000,1);
count = 0;
innnercount = 1;
while (count <5000)

    x = xmin+ (xmax-xmin)*rnds(innnercount); %Random number between xmin & xmax http://ocw.uci.edu/upload/files/mae10_w2011_lecture13.pdf
    y =  ymax*rnds(innnercount+1); %Random number with upper limit ymax
    p_exp_x = lmda*exp(-1*x*lmda);%P(x)
    if y < p_exp_x
        count = count +1;
        rnds_from_exp(count) = x;
    end
    innnercount = mod(innnercount +2,10000);
end
figure;
histogram(rnds_from_exp);
title('Histogram of the generated random numbers from exponential distribution ');
saveas(gcf,'histexp.png')
%Same question for generating a bivariate Gaussian distribution, 
%defined as in: http://www.statisticshowto.com/bivariate-normal-distribution/, 
%with ?_1=0.5, ?_2=1, ?_1=0.75, ?_2=1, and ?=0.7    
xmin = 0; %the lower limit of our domain
xmax = 1; %the upper limit of our domain
ymax = 1; %max value of the distirbution
u1=0.5; u2=1; s1=0.75; s2=1;  ru=0.7;

rnds_from_multivariate = zeros(5000,2);
count = 0;
innnercount = 1;
while (count <5000)

    x1 = xmin+ (xmax-xmin)*rnds(innnercount); %Random number between xmin & xmax http://ocw.uci.edu/upload/files/mae10_w2011_lecture13.pdf
    x2 = xmin+ (xmax-xmin)*rnds(innnercount+1);
    y =  ymax*rnds(innnercount+2); %Random number with upper limit ymax
    x11 = x1-u1;
    x22 = x2-u2;
    z = (x11*x11)/(s1*s1) +(x22*x22)/(s2*s2) - (2*ru*x11*x22)/(s1*s2);
	p_bivar_x=(1/(2*pi*s1*s2*sqrt(1-ru*ru)))* exp( (-1*z)/(2*(1-ru*ru)));%P(x)
    if y < p_bivar_x
        count = count +1;
        rnds_from_multivariate(count,1) = x1;
        rnds_from_multivariate(count,2) = x2;

    end
    innnercount = mod(innnercount +3,10000);
    if innnercount <1
        innnercount = 1;
    end
 
end
figure;
hist(rnds_from_multivariate); 
title('Histogram of the generated random numbers from bivariate Gaussian distribution ');
saveas(gcf,'histbi.png')
%II-Uncertainty propagation:
 
%Use rejections sampling 
a0 = 0.03;b0=0.09;delt_a = 0.02;delt_b=0.02;
delt_a_b = delt_a*delt_b;
inv_delt_a_b = 1/delt_a_b;
a_up = a0+delt_a/2;
a_do = a0-delt_a/2;
b_up = b0+delt_b/2;
b_do = b0-delt_b/2;

count = 0;
innnercount = 1;
ymax = inv_delt_a_b;
rnds_from_a_b = zeros(1000,2);
while (count <1000)

    a = a_up+ (a_up-a_do)*rnds(innnercount); %Random number between xmin & xmax http://ocw.uci.edu/upload/files/mae10_w2011_lecture13.pdf
    b = b_up+ (b_up-b_do)*rnds(innnercount+1);
    y =  ymax*rnds(innnercount+2); %Random number with upper limit ymax
	p_a_b = inv_delt_a_b ;%P(a,b)
    if y < p_a_b
        count = count +1;
        rnds_from_a_b(count,1) = a;
        rnds_from_a_b(count,2) = b;

    end
    innnercount = mod(innnercount +3,10000);
    if innnercount <1
        innnercount = 1;
    end
 
end
hist(rnds_from_a_b);
title('Histogram of the generated random numbers from bivariate initial density distribution ');
saveas(gcf,'histinitdistr.png')
%%%%%%%%%%%%%%%%%%Godunov Scheme
kc1 = 0.02;
kmax = 0.2;
v = 30;
w = 5;
kc2 = ((kc1*v)/(-1*w))+ kmax;
default_k = 0;

%Godunov scheme - Solution 
cell_width = 100; %Width of each cell in x, space step
delta_x = cell_width;
cell_x = 10; %Number of cells
delta_t = 2; %Minimal time step that satisifies CFL condition
time_scale = 10; %Width of time
cell_t = time_scale / delta_t; %Number of cells in time
fac_t_x = delta_t/delta_x;
k = zeros(cell_t,cell_x,1000); %Grid, 1000 = r smaples
%k_inv = zeros(cell_t,cell_x);

%Apply initial condition 
for x = 1:cell_x
    k(1,x,:) = density_mapper(x*delta_x,rnds_from_a_b);
end

for t = 2:cell_t
    
    for x = 1:cell_x
      
        k_t_nv_one_x = k(t-1,x,:);
        s_k_t_nv_one_x = supply(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        d_k_t_nv_one_x = demand(k_t_nv_one_x,kc1,kc2,v,w,kmax);
        
        
        if x-1 < 1
            d_k_t_nv_one_x_nv_one = upstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_nv_one = k(t-1,x-1,:);
            d_k_t_nv_one_x_nv_one = demand(k_t_nv_one_x_nv_one,kc1,kc2,v,w,kmax);

        end
        if x+1 > cell_x
            s_t_nv_one_x_pl_one = downstream_boundary(t*delta_t);

        else
            k_t_nv_one_x_pl_one = k(t-1,x+1,:);
            s_t_nv_one_x_pl_one = supply(k_t_nv_one_x_pl_one,kc1,kc2,v,w,kmax);

        end
        
        for i = 1:1000
        k(t,x,i) = k_t_nv_one_x(i) +fac_t_x * (min(d_k_t_nv_one_x_nv_one(i),s_k_t_nv_one_x(i)) - min(d_k_t_nv_one_x(i),s_t_nv_one_x_pl_one(i)));
        end    
    end  
end



for j = 1:5
figure;

for i = 1:10
    subplot_tight(5,2,i, [0.09 0.09])
    %subplot(5,2,i);
    histogram(k(j,i,:));
    title(sprintf('cell%d',i));


end
currentFigure = gcf;
title(currentFigure.Children(end), ['Godunov Scheme at t =' num2str(j)]);

saveas(gcf,sprintf('G%d.png',j));

end



%FUNCTIONS


function dpd = downstream_boundary(t)
    if (t>=0)&&(t<=10)
        for i= 1:1000
            dpd(i) = 0.2;
        end
  
    end
end
function upd = upstream_boundary(t)
    if (t>=0)&&(t<=10)
        for i= 1:1000
            upd(i) = 0.1;
        end
    
    end
end
function k_0_x = density_mapper(x,rv)
%This function is a helper function to map different constant densities
%correspending to different space values 
%x is a space value, k_0_x = intial  density value at point x
    if (x>=0)&&(x<500)
        k_0_x = rv(:,1);%a distirbutin
    elseif (x>=500)%&&(x<1000)
        k_0_x = rv(:,2); %b distriubiton 
    
    end

end

function tfd = trapezoidal_fundamental_diagram(k,kc1,kc2,v,w,kmax)
%This function is the  Trapezoidal fundamental diagram, defined by ?(k)=k v for k<=kc1 ,
%?(k)=kc1 v for k>kc1 and k<=kc2  and ?(k)=-w(k-kmax) for k>kc2 
    if k <= kc1
       tfd = k*v;
    elseif (k>kc1)&&(k<=kc2)
        tfd = kc1*v;
    elseif k > kc2
        tfd = -1*w*(k-kmax);
    else
        tfd = inf;
    end
end

function dmn = demand(k,kc1,kc2,v,w,kmax)
%This function is the  demand function, defined by 
%D(k)= (?(k),if k?k_c1, ?(k_c1 )=vk_c1,if k>k_c1  and k?k_c2, ?(k_c2 )=vk_c1,if k>k_c2 
indx =1;
    if k(indx) <= kc1
       dmn(indx) = trapezoidal_fundamental_diagram(k(indx),kc1,kc2,v,w,kmax);
    elseif (k(indx)>kc1)&&(k(indx)<=kc2)
        dmn(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        dmn(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    end
for i = 2:1000
    indx = i;
    if k(indx) <= kc1
       dmn(indx) = trapezoidal_fundamental_diagram(k(indx),kc1,kc2,v,w,kmax);
    elseif (k(indx)>kc1)&&(k(indx)<=kc2)
        dmn(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        dmn(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    end    
end

end

function spy = supply(k,kc1,kc2,v,w,kmax)
%This function is the  supply function, defined by 
%S(k)= ?(k_c1 )=vk_c1,if k?k_c1,?(k_c2 )=vk_c1,if k>k_c1  and k?k_c2 , ?(k_ ),if k>k_c2 
indx =1;

    if k(indx) <= kc1
       spy(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    elseif (k(indx)>kc1)&&(k(indx)<=kc2)
        spy(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        spy(indx) = trapezoidal_fundamental_diagram(k(indx),kc1,kc2,v,w,kmax);
    end
for i = 2:1000
   indx =i;
    if k(indx) <= kc1
       spy(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    elseif (k(indx)>kc1)&&(k(indx)<=kc2)
        spy(indx) = trapezoidal_fundamental_diagram(kc1,kc1,kc2,v,w,kmax);
    else%if k > kc2
        spy(indx) = trapezoidal_fundamental_diagram(k(indx),kc1,kc2,v,w,kmax);
    end
    
end
    
end

function vargout=subplot_tight(m, n, p, margins, varargin)
%% subplot_tight
% A subplot function substitude with margins user tunabble parameter.
%
%% Syntax
%  h=subplot_tight(m, n, p);
%  h=subplot_tight(m, n, p, margins);
%  h=subplot_tight(m, n, p, margins, subplotArgs...);
%
%% Description
% Our goal is to grant the user the ability to define the margins between neighbouring
%  subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the
%  margins between subplots can reach 40% of figure area, which is pretty lavish. While at
%  the begining the function was implememnted as wrapper function for Matlab function
%  subplot, it was modified due to axes del;etion resulting from what Matlab subplot
%  detected as overlapping. Therefore, the current implmenetation makes no use of Matlab
%  subplot function, using axes instead. This can be problematic, as axis and subplot
%  parameters are quie different. Set isWrapper to "True" to return to wrapper mode, which
%  fully supports subplot format.
%
%% Input arguments (defaults exist):
%   margins- two elements vector [vertical,horizontal] defining the margins between
%        neighbouring axes. Default value is 0.04
%
%% Output arguments
%   same as subplot- none, or axes handle according to function call.
%
%% Issues & Comments
%  - Note that if additional elements are used in order to be passed to subplot, margins
%     parameter must be defined. For default margins value use empty element- [].
%  - 
%
%% Example
% close all;
% img=imread('peppers.png');
% figSubplotH=figure('Name', 'subplot');
% figSubplotTightH=figure('Name', 'subplot_tight');
% nElems=17;
% subplotRows=ceil(sqrt(nElems)-1);
% subplotRows=max(1, subplotRows);
% subplotCols=ceil(nElems/subplotRows);
% for iElem=1:nElems
%    figure(figSubplotH);
%    subplot(subplotRows, subplotCols, iElem);
%    imshow(img);
%    figure(figSubplotTightH);
%    subplot_tight(subplotRows, subplotCols, iElem, [0.0001]);
%    imshow(img);
% end
%
%% See also
%  - subplot
%
%% Revision history
% First version: Nikolay S. 2011-03-29.
% Last update:   Nikolay S. 2012-05-24.
%
% *List of Changes:*
% 2012-05-24
%  Non wrapping mode (based on axes command) added, to deal with an issue of disappearing
%     subplots occuring with massive axes.

%% Default params
isWrapper=false;
if (nargin<4) || isempty(margins)
    margins=[0.04,0.04]; % default margins value- 4% of figure
end
if length(margins)==1
    margins(2)=margins;
end

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  


height=(1-(m+1)*margins(1))/m; % single subplot height
width=(1-(n+1)*margins(2))/n;  % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width

merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
pos=[merged_left, merged_bottom, merged_width, merged_height];


if isWrapper
   h=subplot(m, n, p, varargin{:}, 'Units', 'Normalized', 'Position', pos);
else
   h=axes('Position', pos, varargin{:});
end

if nargout==1
   vargout=h;
end
end
