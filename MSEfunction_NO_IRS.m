function Result=MSEfunction_NO_IRS(G,Nr,P,Hre_u,sigma2,Ns,D)

I=eye(Ns);
II=eye(Nr);
H=Hre_u;
Result=real(trace(G*(P*H*H'+sigma2*II)*G'-G*P*H*D'-D*H'*P*G'+2*P*I));

end