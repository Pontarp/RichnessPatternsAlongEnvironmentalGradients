%This function calculates fitness for each of the individuals in the system
%One fitness value per individual will be outputed in a vector. 

function fitness = fitfunc(zp,Kp,Up,zz,N,r,sig_a,sig_K)

%Calculate Neff, the effective no of competitior
alphaij=exp(-(zp-zz).^2/2/sig_a^2);

Neff=sum(alphaij.*N); 

%Beräkna finess
fitness= 1 + r*(1-Neff/Kp);
fitness=fitness*Up; 
