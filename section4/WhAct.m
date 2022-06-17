function y = WhAct(bf, bm, mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im)
%Wp calculates the fitness of the host at equilibrium Sf, Sm, If, Im
%   gammas and bf and bm (formerly just b) are _mutant_ values 
F=zeros(4);
Vinv=zeros(4);
N=Sf+Sm+If+Im;

% now bf and bm are arrays
bSfSm1 = bf(1,1);
bSfIm1 = bf(1,2); 
bIfSm1 = bf(2,1);
bIfIm1 = bf(2,2);

bSfSm2 = bm(1,1);
bSfIm2 = bm(1,2);
bIfSm2 = bm(2,1);
bIfIm2 = bm(2,2);

% populate F column by column

% for this col bSfSm and bSfIm with mutant arg in 1st pos (mutant f trait)
F(1,1)= (bSfSm1 * Sm + bSfIm1 * Im) / N;
F(2,1)= (bSfSm1 * Sm + bSfIm1 * Im) / N;

% for this col bSfSm and bIfSm with mutant arg in 2nd pos (mutant m trait)
F(1,2)= (Sf * bSfSm2 + (1-v)*If * bIfSm2) / N;
F(2,2)= (Sf * bSfSm2 + (1-v)*If * bIfSm2) / N;
F(3,2)= v *If * bIfSm2 / N;
F(4,2)= v *If * bIfSm2 / N;

% for this col bIfSm and bIfIm with mutant arg in 1st pos (mutant f trait)
F(1,3)= (1-v) * ( bIfSm1 * Sm + bIfIm1 * Im) / N;
F(2,3)= (1-v) * ( bIfSm1 * Sm + bIfIm1 * Im) / N;
F(3,3)= v * ( bIfSm1 * Sm + bIfIm1 * Im) / N;
F(4,3)= v * ( bIfSm1 * Sm + bIfIm1 * Im) / N;

% for this col bSfIm and bIfIm with mutant arg in 2nd pos (mutant m trait)
F(1,4)= (Sf * bSfIm2 + (1-v) * If * bIfIm2) / N;
F(2,4)= (Sf * bSfIm2 + (1-v) * If * bIfIm2) / N;
F(3,4)= v * If * bIfIm2 / N;
F(4,4)= v * If * bIfIm2 / N;

F = 0.5 * F;

% populate Vinv column by column
tauf=( mu*N*gammaf + (mu*N+alphaf) * (mu*N + betaff*If + betafm*Im) )^(-1);
taum=( mu*N*gammam + (mu*N+alpham) * (mu*N + betamf*If + betamm*Im) )^(-1);

Vinv(1,1)= ( mu * N + alphaf + gammaf ) * tauf;
Vinv(3,1)= ( betaff * If + betafm * Im ) * tauf;

Vinv(2,2)= ( mu * N + alpham + gammam ) * taum;
Vinv(4,2)= ( betamf * If + betamm * Im ) * taum;

Vinv(1,3)= gammaf * tauf;
Vinv(3,3)= ( mu * N + betaff * If + betafm * Im ) * tauf;

Vinv(2,4)= gammam * taum;
Vinv(4,4)= ( mu * N + betamf * If + betamm * Im ) * taum;

y = max(abs(eig(F*Vinv)));
end

