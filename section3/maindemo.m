clear all
% maindemo.m

tic

% Open a csv data file and print the header line
fileID = fopen('maindemo.csv', 'w');
fprintf(fileID, 'bmax,mu,betaC1,betaC2,');
fprintf(fileID, 'base_c,diff_c,prVert,');
fprintf(fileID, 'alphaf,alpham,gammaf,gammam,p_ESS?,h_ESS?,NormGrad,Iter,');
fprintf(fileID, 'Sfbar,Smbar,Ifbar,Imbar\n');

% set parameter values outside of parameter loop
mu=0.1;
bmax = 2;
c0 = 0.1;
v=0.25;
dc = -0.05; % twice the difference between cf and cm

        
% establish tradeoff faced by pathogen
betaC1 = 2.0;
betaC2 = 4.0;
betaXf =@(a) betaC1 * a./(a + betaC2);
betaXm =@(a) betaC1 * a./(a + betaC2);

% establish tradeoff faced by host
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

    % establish host fecundity b
    b = bmax * pf(gammaf) * pm(gammam);

    % determine resident equilibrium values
    [y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
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
    Wpfup = Wp(b,mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im);
    Wpfdn = Wp(b,mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im);
    grad(1) = Wpfup - Wpfdn;

    % determine selection gradient acting on alpham
    Wpmup = Wp(b,mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im);
    Wpmdn = Wp(b,mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im);
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

    bfup = bmax * pf(gfup) * pm(gammam);
    bfdn = bmax * pf(gfdn) * pm(gammam);

    bmup = bmax * pf(gammaf) * pm(gmup);
    bmdn = bmax * pf(gammaf) * pm(gmdn);

    % determine selection gradient acting on gammaf
    Whfup = Wh(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    Whfdn = Wh(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    grad(3) = Whfup - Whfdn;

    % determine selection gradient acting on gammam
    Whmup = Wh(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    Whmdn = Wh(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
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
[y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
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

Wp0 =  Wp(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
% Wp0 should be 1, which can be verified by un-commenting next line
% display(Wp0);

Wpfup = Wp(b,mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im);
Wpfdn = Wp(b,mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im);
Wpmup = Wp(b,mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im);
Wpmdn = Wp(b,mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im);

Wpupup = Wp(b,mu,v,gammaf,gammam,afup,amup,betaffup,betafmup,betamfup,betammup,Sf,Sm,If,Im);
Wpdndn = Wp(b,mu,v,gammaf,gammam,afdn,amdn,betaffdn,betafmdn,betamfdn,betammdn,Sf,Sm,If,Im);

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

% establish fecundity(redundant here but can cut paste)
b = bmax * pf(gammaf) * pm(gammam);

% determine resident equilibrium values (redundant here but can cut paste)
[y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
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

bfup = bmax * pf(gfup) * pm(gammam);
bfdn = bmax * pf(gfdn) * pm(gammam);

bmup = bmax * pf(gammaf) * pm(gmup);
bmdn = bmax * pf(gammaf) * pm(gmdn);

Wh0 = Wh(b,b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
% Wp0 should be 1, which can be verified by un-commenting next line
% display(Wh0);

Whfup = Wh(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whfdn = Wh(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whmup = Wh(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whmdn = Wh(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);

Whupup= Wh(bfup,bmup,mu,v,gfup,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whdndn= Wh(bfdn,bmdn,mu,v,gfdn,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);

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