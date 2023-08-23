function Result=MSEfunction_NO_IRS_NO_PNC(G,Nr,P,Hre_u,sigma2,Ns)

I=eye(Ns);
II=eye(Nr);
H=Hre_u;
Result=real(trace(G*(H*P*I*H'+sigma2*II)*G'-G*P*H-H'*P*G'+P*I));

end