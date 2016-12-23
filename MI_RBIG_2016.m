% Multi-information estimation using RBIG
%
% The user can choose two orthogonal transforms:
%
%      'PCA' = PCA
%      'RND' = Random Rotations
%
% USE:
%
% [MI] = MI_RBIG_2016(dat,N_lay,transform,porc,precision)
%
% INPUTS:
% dat = data ( #dimensions , #samples );aim
% N_lay = number of layers (default N_lay = 1000);
% porc = extra domain percentage (default porc = 10)
% precision = number of points for the marginal PDFs estimation (default precision = 1000)
% transformation = linear transformation applied ('RND','PCA' default transformation = 'PCA')
%
% OUTPUTS
% MI = Multi-information
% MIs = Multi-information reduction at each layer.
% datT = Gaussianized data.
%
% e.g.
%
% dat = rand(5)*(rand(5,1000).^2);
% N_lay = 50;
% porc = 1;
% precision = 1000;
% transformation = 'PCA';
% [MI, MIs] = MI_RBIG_2016(dat,N_lay,transformation,porc,precision);
%
%
% Citation:
% Iterative Gaussianization: from ICA to Random Rotations. 
% V. Laparra, G. Camps & J. Malo 
% IEEE Transactions on Neural Networks, 22:4, 537 - 549, (2011)
%


function [MI MIs datT] = MI_RBIG_2016(dat,N_lay,transformation,porc,precision)

if nargin < 5, precision = 1000; end
if nargin < 4, porc = 10; end
if nargin < 3, transformation = 'PCA'; end
if nargin < 2, N_lay = 1000; end

DIM = size(dat,1);
Nsamples = size(dat,2);

for rep=1:1000
    [p R] = hist(randn(1,round(Nsamples)),round(sqrt(Nsamples)));
    delta = R(3)-R(2);
    hx = entropy_mm(p)+log2(delta);
    
    [p R] = hist(randn(1,round(Nsamples)),round(sqrt(Nsamples)));
    delta = R(3)-R(2);
    hy = entropy_mm(p)+log2(delta);
    ee(rep) = hy - hx;
end
tol_m = mean(ee);
tol_d = std(ee);
%pause(10)

for n = 1:N_lay
    % marginal gaussianization
    for dim = 1:DIM
        [dat(dim,:) T]= marginal_gaussianization(dat(dim,:),porc,precision);
    end

    dat_aux = dat;

    % rotation
    if transformation == 'RND'
        V = rand(DIM);
        V = V * inv(sqrtm(V'*V)); % orthogonalization
        V = V / (abs(det(V))^(1/size(V,1))); % normalization
        dat = V'*dat;

    elseif transformation == 'PCA'
        C = dat*dat'/size(dat,2);
        [V D] = eig(C);
        dat = V'*dat;
    end

    % multi-information reduction
    Res(n).I = information_reduction_LT(dat,dat_aux,tol_m,tol_d);
    if (n>100),
        auxi = cat(1,Res.I);
        if sum(abs(auxi(end-90:end))) == 0
            break
        end
    end
end
datT = dat;

MIs = cat(1,Res.I);
MI = sum(MIs);



function I=information_reduction_LT(X,Y,tol_m,tol_d)

[DIM,Nsamples] = size(X);

for n=1:DIM
    [p R]=hist(X(n,:),sqrt(Nsamples));
    delta = R(3)-R(2);
    hx(n)=entropy_mm(p)+log2(delta);

    [p R]=hist(Y(n,:),sqrt(Nsamples));
    delta = R(3)-R(2);
    hy(n)=entropy_mm(p)+log2(delta);
end

%disp(hy)
%disp(hx)
I = sum(hy) - sum(hx);
II = sqrt(sum((hy - hx).^2));
p = 0.25;
%disp(sqrt(DIM*((p*tol_d.^2))))
%pause(1)
if abs(II)<sqrt(DIM*((p*tol_d.^2)))
    I=0;
end
    

function H = entropy_mm(p)
% mle estimator with miller-maddow correction
c = 0.5 * (sum(p>0)-1)/sum(p);  % miller maddow correction
p = p/sum(p);               % empirical estimate of the distribution
idx = p~=0;
H = -sum(p(idx).*log2(p(idx))) + c;     % plug-in estimator of the entropy with correction

function [x_gauss T] = marginal_gaussianization(x,porc,precision) 

[x_unif T] = marginal_uniformization(x,porc,precision);
x_gauss = norminv(x_unif);

function [x_lin T] = marginal_uniformization(x,porc,precision)

if nargin == 2, precision = 1000; end

p_aux = (porc/100)*abs(max(x)-min(x));
R_aux = linspace(min(x),max(x),2*sqrt(length(x))+1);

R = mean([R_aux(1:end-1); R_aux(2:end)]);

[p,R] = hist(x,R);   

delta_R = R(3)-R(2); % posiblemente haiga que dividirlo por 2
T.R_ant = [R(1)-delta_R R R(end)+delta_R];
T.p_ant = [0 p./(sum(p)*(R(4)-R(3))) 0];

C = [(cumsum(p))];
N = max(C);

C = (1-1/N)*C/N;

incr_R = (R(2)-R(1))/2;

R = [min(x)-p_aux min(x)  R(1:end)+incr_R max(x)+p_aux+incr_R];
C = [0 1/N C 1]; 

Range_2 = linspace(R(1),R(end),precision);
C_2 = made_monotonic(interp1(R,C,Range_2));
C_2 = C_2/max(C_2);
x_lin = interp1(Range_2,C_2,x); 

T.C = C_2;
T.R = Range_2;


function fn = made_monotonic(f)

fn=f;
for nn=2:length(fn)
    if(fn(nn)<=fn(nn-1))
        if abs(fn(nn-1))>1e-14
            fn(nn)=fn(nn-1)+1.0e-14;
        elseif fn(nn-1)==0
            fn(nn)= 1e-80;
        else
            fn(nn)=fn(nn-1)+10^(log10(abs(fn(nn-1))));
        end
    end
end
