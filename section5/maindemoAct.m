clear all
tic

% Open a csv data file and print the header line
fileID = fopen('maindemoAct.csv', 'w');
fprintf(fileID, 'bmax,mu,betaC1,betaC2,');
fprintf(fileID, 'base_c,diff_c,prVert,');
fprintf(fileID, 'alphaf,alpham,gammaf,gammam,p_ESS?,h_ESS?,NormGrad,Iter,');
fprintf(fileID, 'Sfbar,Smbar,Ifbar,Imbar\n');

% set parameter values outside of parameter loop
mu=0.1;
bmax = 3.5;
c0 = 1;
v=0.25;
dc=-0.05;

% establish tradeoff faced by pathogen
betaC1 = 3.1;
betaC2 = 4.0;
betaXf =@(a) betaC1 * a./(a + betaC2);
betaXm =@(a) betaC1 * a./(a + betaC2);

% establish tradeoff faced by host
% pX is 1 when gamma = 0
pf =@(gamma) exp( -(c0+dc) * gamma.^2 );
pm =@(gamma) exp( -(c0-dc) * gamma.^2 );

% begin with a poor guess at ESS values
alphaf = rand(); 
alpham = rand(); 
gammaf = rand();
gammam = rand();

% update guesses iteratively to find convergence-stable 4-tuple 
tol = 1e-07;
maxiter = 1e05;
iter=0;

grad = ones(4,1);
while norm(grad) > tol
    iter = iter + 1;
    if iter > maxiter
        break
    end

    % Resident Populations
    %%%%%%%%%%%%%%%%%%%%%%%

    % establish pathogen transmissibilities
    betaff=betaXf(alphaf);
    betamf=betaXf(alphaf);

    betafm=betaXm(alpham);
    betamm=betaXm(alpham);

    % establish host fecundity b, modified!
    % b = [ [ bSfSm, SfIm ]; [bIfSm , bImIf] ]
    b = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gammam)] ];

    % determine resident equilibrium values
    [y,dydt]=ResEquilAct(b(1,1), b(1,2), b(2,1), b(2,2) ,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm); % modified
    Sf = y(1);
    Sm = y(2);
    If = y(3);
    Im = y(4);

    % Pathogen 
    %%%%%%%%%%%

    % determine mutant pathogen traits
    delta = 1e-01;
    afup = alphaf + 0.5 * delta;
    afdn = alphaf - 0.5 * delta;
    amup = alpham + 0.5 * delta;
    amdn = alpham - 0.5 * delta;

    while or(afdn<0, amdn<0) 
        delta = 0.1 * delta;
        afup = alphaf + 0.5 * delta;
        afdn = alphaf - 0.5 * delta;
        amup = alpham + 0.5 * delta;
        amdn = alpham - 0.5 * delta;
    end

    betaffup=betaXf(afup);
    betaffdn=betaXf(afdn); 

    betamfup=betaXf(afup);
    betamfdn=betaXf(afdn);

    betafmup=betaXm(amup);
    betafmdn=betaXm(amdn);

    betammup=betaXm(amup);
    betammdn=betaXm(amdn);

    % determine selection gradient acting on alphaf 
    % need different b functions but only infected females needed
    % b = [ [ bSfSm, SfIm ]; [bIfSm , bImIf] ]
    Wpfup = WpAct(b(2,1), b(2,2),mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im); % modified!
    Wpfdn = WpAct(b(2,1), b(2,2),mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im); % modified!
    grad(1) = Wpfup - Wpfdn;

    % determine selection gradient acting on alpham
    Wpmup = WpAct(b(2,1), b(2,2),mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im); % modified!
    Wpmdn = WpAct(b(2,1), b(2,2),mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im); % modified!
    grad(2) = Wpmup - Wpmdn;

    % Host
    %%%%%%%

    % determine host mutant traits
    delta = 1e-01;
    gfup = gammaf + 0.5 * delta;
    gfdn = gammaf - 0.5 * delta;
    gmup = gammam + 0.5 * delta;
    gmdn = gammam - 0.5 * delta;

    while or(gfdn<0, gmdn<0) 
        delta = 0.1 * delta;
        gfup = gammaf + 0.5 * delta;
        gfdn = gammaf - 0.5 * delta;
        gmup = gammam + 0.5 * delta;
        gmdn = gammam - 0.5 * delta;
    end

    % need matrices for mutant hosts
    % b = [ [ bSfSm, SfIm ]; [bIfSm , bImIf] ]
    bfup = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gfup), bmax * pf(gfup) * pm(gammam)] ];
    bfdn = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gfdn), bmax * pf(gfdn) * pm(gammam)] ];

    bmup = [ [bmax, bmax * pm(gmup)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gmup)] ];
    bmdn = [ [bmax, bmax * pm(gmdn)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gmdn)] ];

    % determine selection gradient acting on gammaf
    Whfup = WhAct(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
    Whfdn = WhAct(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
    grad(3) = Whfup - Whfdn;

    % determine selection gradient acting on gammam
    Whmup = WhAct(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
    Whmdn = WhAct(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
    grad(4) = Whmup - Whmdn;

    % Update Pathogen and Host
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    step = 1;
    alphaf = alphaf + step * grad(1);
    alpham = alpham + step * grad(2);
    gammaf = gammaf + step * grad(3);
    gammam = gammam + step * grad(4);
end

% check that pathogen is at ESS

% establish transmissibilities (redundant here but can cut paste)
betaff=betaXf(alphaf);
betamf=betaXf(alphaf);

betafm=betaXm(alpham);
betamm=betaXm(alpham);

% determine resident equilibrium values (redundant here but can cut paste)
[y,dydt]=ResEquilAct(b(1,1), b(1,2), b(2,1), b(2,2) ,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm); % modified
Sf = y(1);
Sm = y(2);
If = y(3);
Im = y(4);

% determine mutant traits
delta = 1e-01;
afup = alphaf + 0.5 * delta;
afdn = alphaf - 0.5 * delta;
amup = alpham + 0.5 * delta;
amdn = alpham - 0.5 * delta;

while or(afdn<0, amdn<0) 
    delta = 0.1 * delta;
    afup = alphaf + 0.5 * delta;
    afdn = alphaf - 0.5 * delta;
    amup = alpham + 0.5 * delta;
    amdn = alpham - 0.5 * delta;
end

betaffup=betaXf(afup);
betaffdn=betaXf(afdn); 

betamfup=betaXf(afup);
betamfdn=betaXf(afdn);

betafmup=betaXm(amup);
betafmdn=betaXm(amdn);

betammup=betaXm(amup);
betammdn=betaXm(amdn);

Wp0 =  WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); %modified
% Wp0 should be 1, which can be verified by un-commenting next line
% display(Wp0);

Wpfup = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im);%modified
Wpfdn = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im);%modified
Wpmup = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im);%modified
Wpmdn = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im);%modified

Wpupup = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,afup,amup,betaffup,betafmup,betamfup,betammup,Sf,Sm,If,Im);%modified
Wpdndn = WpAct(b(2,1),b(2,2),mu,v,gammaf,gammam,afdn,amdn,betaffdn,betafmdn,betamfdn,betammdn,Sf,Sm,If,Im);%modified

H=zeros(2);
H(1,1) = (Wpfup - 2 * Wp0 + Wpfdn) / ((0.5 * delta)^2);
H(2,2) = (Wpmup - 2 * Wp0 + Wpmdn) / ((0.5 * delta)^2);
H(1,2) = (Wpupup - Wpfup - Wpmup + 2 * Wp0 - Wpfdn - Wpmdn + Wpdndn) / ( 2 * 0.5 * delta * 0.5 *delta );
H(2,1) = (Wpupup - Wpfup - Wpmup + 2 * Wp0 - Wpfdn - Wpmdn + Wpdndn) / ( 2 * 0.5 * delta * 0.5 *delta ); 
evals = eig(H);
if and(evals(1)<0, evals(2)<0)
    pESS = 'True';
else
    pESS = 'False';
end

% check Host ESS

% establish fecundity(redundant here but can cut paste) modified
% b = [ [ bSfSm, SfIm ]; [bIfSm , bImIf] ]
b = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gammam)] ];

% determine resident equilibrium values (redundant here but can cut paste)
[y,dydt]=ResEquilAct(b(1,1), b(1,2), b(2,1), b(2,2) ,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm); % modified
Sf = y(1);
Sm = y(2);
If = y(3);
Im = y(4);

% determine mutant traits
delta = 1e-01;
gfup = gammaf + 0.5 * delta;
gfdn = gammaf - 0.5 * delta;
gmup = gammam + 0.5 * delta;
gmdn = gammam - 0.5 * delta;

while or(gfdn<0, gmdn<0) 
    delta = 0.1 * delta;
    gfup = gammaf + 0.5 * delta;
    gfdn = gammaf - 0.5 * delta;
    gmup = gammam + 0.5 * delta;
    gmdn = gammam - 0.5 * delta;
end

% b = [ [ bSfSm, SfIm ]; [bIfSm , bImIf] ]
bfup = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gfup), bmax * pf(gfup) * pm(gammam)] ];
bfdn = [ [bmax, bmax * pm(gammam)]; [bmax * pf(gfdn), bmax * pf(gfdn) * pm(gammam)] ];

bmup = [ [bmax, bmax * pm(gmup)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gmup)] ];
bmdn = [ [bmax, bmax * pm(gmdn)]; [bmax * pf(gammaf), bmax * pf(gammaf) * pm(gmdn)] ];

Wh0 = WhAct(b,b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
% Wp0 should be 1, which can be verified by un-commenting next line
% display(Wh0);

Whfup = WhAct(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
Whfdn = WhAct(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);% modified
Whmup = WhAct(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im); % modified
Whmdn = WhAct(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);% modified

Whupup= WhAct(bfup,bmup,mu,v,gfup,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);% modified
Whdndn= WhAct(bfdn,bmdn,mu,v,gfdn,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);% modified

H=zeros(2);
H(1,1) = (Whfup - 2 * Wh0 + Whfdn) / ((0.5 * delta)^2);
H(2,2) = (Whmup - 2 * Wh0 + Whmdn) / ((0.5 * delta)^2);
H(1,2) = (Whupup - Whfup - Whmup + 2 * Wh0 - Whfdn - Whmdn + Whdndn) / ( 2 * 0.5 * delta * 0.5 *delta );
H(2,1) = (Whupup - Whfup - Whmup + 2 * Wh0 - Whfdn - Whmdn + Whdndn) / ( 2 * 0.5 * delta * 0.5 *delta ); 
evals = eig(H);
if and(evals(1)<0, evals(2)<0)
    hESS = 'True';
else
    hESS = 'False';
end

% write the data to the datafile
fprintf(fileID, '%5.3f,%5.3f,%5.3f,%5.3f,', bmax, mu, betaC1, betaC2);
fprintf(fileID, '%8.6f,%8.6f,%8.6f,', c0, dc, v);
fprintf(fileID, '%8.6f,%8.6f,%8.6f,%8.6f,%5s,%5s,%g,%7d,', ...
    alphaf,alpham,gammaf,gammam,pESS,hESS,norm(grad),iter);
fprintf(fileID, '%8.6f,%8.6f,%8.6f,%8.6f\n', Sf,Sm,If,Im);

toc
% close the datafile
fclose(fileID);