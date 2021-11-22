function [E_prop, E_ft, E_ft_prop, kx]=propagator(E0,k0,delta_z,pad,res)

E0=[zeros(pad,1); E0(:,1); zeros(pad,1)];

N=length(E0);
npos=floor((N-1)/2);
nneg=N-1-npos;

dx=res;         

E_ft=fft(E0);                   % spectrum at z=0 distance
dk=2*pi/(N*dx); 
kx=([0:npos -nneg:-1])*dk;      % k vector

q = sqrt(k0.^2 - kx.^2).';        

P = exp(1i*q*delta_z);
P(abs(kx) > real(k0)) = 0;
E_ft_prop = E_ft.*P;
E_prop = ifft(E_ft_prop);

end