function Result=MSEfunction_G(theta,g,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2,Ns,D)

I=eye(Ns);
II=eye(Nr);
Theta=diag(exp(j*[theta(1:M)]));
H=Hre_irs*Theta*Hirs_u+Hre_u;
G=[g(1:Nr)+1i*g(Nr+1:2*Nr);g(2*Nr+1:3*Nr)+1i*g(3*Nr+1:4*Nr)];
Result=real(trace(G*(H*P*I*H'+sigma2*II)*G'-G*H*P*I*D'-D*P*I*H'*G'+D*P*I*D'));

end