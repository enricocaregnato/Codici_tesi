
function demo_sphtri

%--------------------------------------------------------------------------
% Object:
% Demo on hyperinterpolation on a fixed domain over 
% spherical triangle on the unit sphere S2.
% It checks function reconstruction, assuming that 
% noise may be present.
% If no noise is required, set "a=0; sigma=0;" below.
% It performs: hyperinterpolation, filtered hyp.
%--------------------------------------------------------------------------
% Usage:
% >> demo_sphtri
%--------------------------------------------------------------------------
% Dates:
% Written on November 26, 2022: E. Caregnato and A. Sommariva.
%--------------------------------------------------------------------------

clear; clf;

domain_example=1;

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
nV=1:10;

%--------------------------------------------------------------------------
% Noise and lasso parameter.
%--------------------------------------------------------------------------

noise=1;          % 0: no noise, 1: noise
a=2.5;              % defining impulse noise (in experiment 2)
sigma=0; %0.02;   % defining gaussian noise (in experiment 2)
%lambda=10^(-2);   % defining lasso parameter
pos=0;            % extraction type.
domain_structure.domain='spherical-triangle';

% ........................ Main code below ................................


% ....... Special settings .......

% Approximation type parameter "pts_type".
%     case 1, pts_type='Hyperinterpolation full set';
%     case 2, pts_type='Hyperinterpolation compressed set';
%
% Note: the case "2" should be used for mild values of "n", say at most 15.
pts_type=1;

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=1;

% Number of tests for reconstructing polygons.
ntests=10;

% testing functions:
% 1. polynomial of degree n,
% 2. polynomial of degree floor(n/2)-1
% 3. gaussian like exponential
funct_example=3;


% ....... Apply settings to define domain, pointsets, and functions .......

% Domain
vertices=define_domain(domain_example);

% Test points
% [XYZWR,dbox]=define_cub_rule(domain_structure,30);
P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);
XYZWR = cub_sphtri(30,P1',P2',P3',pos);

XR=XYZWR(:,1:end-1); WR=XYZWR(:,end);

% ........ Numerical approximation, varying the degree in "nV" ............
fprintf('\n \t ');
AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

for k=1:length(nV)
    n=nV(k);

    dimpoly=(n+1)^2;

    % ... extract hyperinterpolation set (notice that ade=2*n) ...
    if pts_type == 2 % compressed set
        % fprintf('\n \t * Compressed set')
        XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
        [pts,weights,momerr,dbox] =...
            dCATCH(2*n,XYZW(:,1:3),XYZW(:,4));
    else % full set
        % fprintf('\n \t * Full set')
        XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
        pts=XYZW(:,1:end-1); weights=XYZW(:,end);
    end

    % .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..

    poly_coeffs=[];

    for j=1:ntests

        switch funct_example

            case 1 % test exactness hyperinterpolation
                nexp=n;
                c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
                g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

            case 2 % test exactness filt. hyperinterpolation
                nexp=floor(n/2);
                c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
                g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

            case 3 % exponential type
                c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
                g=@(x,y,z) exp(-(c0*x.^2+c1*y.^2+c2*z.^2+c3));

        end

        % ... evaluate function to approximate ...
        g_pts=feval(g,pts(:,1),pts(:,2),pts(:,3));

        % ... Add noise (if present) ...

        if noise
            % add impulse noise
            pert_impulse=0;
            if a > 0
                pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
                while norm(pert_impulse) == 0
                    pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
                end
            end

            % add gaussian noise
            pert_gauss=0;
            if sigma > 0
                var=sigma^2;
                pert_gauss=sqrt(var)*randn(size(g_pts));
                while norm(pert_gauss) == 0
                    pert_gauss=sqrt(var)*randn(size(g_pts));
                end
            end

            % add gaussian + impulse noise
            pert=pert_impulse+pert_gauss;

            % perturbed values
            g_pts_pert=g_pts+pert;
        else
            g_pts_pert=g_pts;
        end

        % ... determine polynomial hyperinterpolant ...
        [coeff0,R,jvec,dbox,degs] = dHYPERFIT2(n,pts,weights,g_pts_pert,...
            [],[],domain_structure,dimpoly);

        if iscell(jvec), degs=degs(jvec{1}); else, degs=degs(jvec); end

        % test hyperinterpolant with or withour filters.
        for ktest=1:2

            switch ktest
                case 1
                    hypermode='hyperinterpolation';
                    parms.lambda=[]; parms.mu=[];
                    coeff=coeff0;
                case 2
                    hypermode='filtered';
                    parms.lambda=[]; parms.mu=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
            end

            % ... evaluate hyperinterpolant at initial pointset ...
            p_XR=dPOLYVAL2(n,coeff,XR,R,jvec,dbox,domain_structure,dimpoly);

            % ... estimating hyperinterpolant error ...
            g_XR=feval(g,XR(:,1),XR(:,2),XR(:,3));

            tolAE=10^(-12)*ones(size(p_XR));
            AEinfV(ktest,j)=norm(g_XR-p_XR,inf); % absolute error (inf norm)
            AE2V(ktest,j)=sqrt(WR'*((g_XR-p_XR).^2)); % absolute error (2 norm)
            beta0V(ktest,j)=sum(abs(coeff) > 0);

        end

    end

    % averages of the errors (vectors 2 x 1)
    AEinfM=mean(AEinfV,2);
    AE2M=mean(AE2V,2);
    beta0M=mean(beta0V,2);

    fprintf('\n       ........ table at degree: %2.0f ........ \n \n ',n);
    HypType=categorical({'hyperint.'; 'filtered'});
    T = table(HypType,AEinfM,AE2M,beta0M); disp(T)

    AEinfMV=[AEinfMV AEinfM]; AE2MV=[AE2MV AE2M]; beta0MV=[beta0MV beta0M];

    % ... parameters for statistics ...
    card_pts(k)=length(weights);

    card_polyspace(k)=(n+1)^2;

end


% ............................ Plot nodes .................................
domain_parms=vertices;
plot_s2(domain_example,domain_parms,pts,[],'',[]);


% ............................ General Statistics .........................


fprintf('\n \t .........................................................................');
fprintf('\n \t ** Domain: %2.0f'); disp(domain_example);
if pts_type == 1
    fprintf('\n \t ** Hyperinterpolation on full pointset \n \t');
else
    fprintf('\n \t ** Hyperinterpolation on compressed pointset \n \t');
end

fprintf('\n \t .........................................................................');
fprintf('\n \t |  n  |  hypset  |  dimpoly  |  ');
fprintf('\n \t .........................................................................');
for k=1:length(nV)
    fprintf('\n \t | %3.0f |   %5.0f  |   %5.0f   |',...
        nV(k),card_pts(k),card_polyspace(k));
end
fprintf('\n \t .........................................................................');
fprintf('\n \t * SPECIFIC PARAMETERS');
fprintf('\n \t .........................................................................');
fprintf('\n \t impulse noise parameter,      a: %1.15g',a);
fprintf('\n \t gaussian noise parameter, sigma: %1.15g',sigma);
fprintf('\n \t averages on # tests            : %3.0f',ntests);
fprintf('\n \t .........................................................................');
fprintf('\n \n');



fprintf('\n \n');








function vertices=define_domain(example)


switch example

    case 0
        % large domain
        vertices=[ 1 0 0; 
            0 1 0;
            0 0 1];

    case 1
        a=0.5; b=0.1; % MEDIUM-LARGE SIZE (Australia)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];

end

