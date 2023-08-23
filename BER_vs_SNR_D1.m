close all;
clear all;
clc;
%Parameter
P=1; %power of the signal
Nr=2; %number of antennas at the relay node
Ns=2; %number of the single antenna users
Nu=1; %number of antennas in single user 
M=32; %number of the reflected elements in IRS
I=eye(Nr);
II=eye(Ns);
D=[1 1; 1 -1];
Ngau=1000; %number of gaussian randomization vectors
SNR_dB=-25:5:20; %dB
SNR=10.^(SNR_dB./10);
sigma2=P./SNR; %sigma squarec
options = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxFunctionEvaluations',1000);
iter=2;
Length=100000; %length of stream of bits
lN=eps; %lowest number is added to avoid NaN
error=1e-4;
error_q=1e-3;
max_iteration=100;
gamma2=1;
% BER=zeros(iter,length(SNR));
MSE_iter_total=zeros(iter,max_iteration);
MSE_iter_total1=zeros(iter,max_iteration);
MSE_iter_total3=zeros(iter,max_iteration);
for ii=1:iter
%channels with zero mean and unit variance
    
%channels from the users 1 and 2 to relay via IRS   
% Hre_u=sqrt(1/2)*(randn(Nr,Ns)+1i*randn(Nr,Ns))/sqrt(Nr);
% Hirs_u=sqrt(1/2)*(randn(M,Ns)+1i*randn(M,Ns))/sqrt(M);
% Hre_irs=sqrt(1/2)*(randn(Nr,M)+1i*randn(Nr,M))/sqrt(Nr*M);  

Hre_u=sqrt(1/2)*(randn(Nr,Ns)+1i*randn(Nr,Ns));
Hirs_u=sqrt(1/2)*(randn(M,Ns)+1i*randn(M,Ns));
Hre_irs=sqrt(1/2)*(randn(Nr,M)+1i*randn(Nr,M));  

H3=Hre_u;
Phi1_her=Hre_irs(1,:)*diag(Hirs_u(:,1));
Phi2_her=Hre_irs(1,:)*diag(Hirs_u(:,2));
Phi3_her=Hre_irs(2,:)*diag(Hirs_u(:,1));
Phi4_her=Hre_irs(2,:)*diag(Hirs_u(:,2));


%channel from the relay to user 3
% Hr_D1=sqrt(1/2)*(randn(Nu,Nr)+1i*randn(Nu,Nr))/sqrt(Nu*Nr);
Hr_D1=sqrt(1/2)*(randn(Nu,Nr)+1i*randn(Nu,Nr));


%channels from user 1 to user 3 via IRS
% Hs1d1=sqrt(1/2)*(randn(Nu,Nu)+1i*randn(Nu,Nu))/sqrt(Nu*Nu);
Hs1d1=sqrt(1/2)*(randn(Nu,Nu)+1i*randn(Nu,Nu));

theta2=(2*pi)*rand(1,M);
%% optimization
for jj=1:length(SNR)
  
%Plotting BER


%Info bits
bits=randi([0,1],1,2*Length);

% Generating BPSK symbols
BPSK_symbols=2*bits-1;

X=reshape(BPSK_symbols,Ns,Length); %BPSK symbols to send 

% calculating true xor at the transmission
x1_xor_x2_transmitted_actual=(-1)*(X(1,:).*X(2,:));
x1_xor_x2_transmitted_actual_bits=(x1_xor_x2_transmitted_actual==1)*1;

x1_actual=(X(1,:)==1)*1; 
x2_actual=(X(2,:)==1)*1; 

%%%% Cases
%Case butterfly network with PNC and optimal phases
%Case1 butterfly network with PNC and zero phases
%Case2 butterfly network with PNC and random phases
%Case3 butterfly network with PNC and quantized phases
%Case4 butterfly network with PNC and without IRS
%Case5 butterfly network with NNC and without IRS


%%defining MSE for each case 

%%%%%%%Case optimizing G and theta SDR Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart=tic;
MSE_previous=0;
MSEdecreasepercent=Inf;
ph=pi*ones(1,M);
% v_phases=exp(1i*ph');
iteration=0;
% while((MSEdecreasepercent>=error)&&(iteration<=max_iteration))
while(iteration<=max_iteration)

theta=ph;
Theta=diag(exp(1i*theta));
H_total=Hre_irs*Theta*Hirs_u+Hre_u;
% H_hat=H_total*inv(D);
% G3=inv(sigma2(jj)*II+H_hat2'*H_hat2)*H_hat2';
G=P*D*H_total'*inv(H_total*H_total'+sigma2(jj)*I);


A=[Phi3_her'*Phi1_her+Phi2_her'*Phi4_her+Phi1_her'*Phi3_her+Phi4_her'*Phi2_her Phi3_her'*H3(1,1)+Phi2_her'*H3(2,2)+Phi1_her'*H3(2,1)+Phi4_her'*H3(1,2); H3(2,1)'*Phi1_her+H3(1,2)'*Phi4_her+H3(1,1)'*Phi3_her+H3(2,2)'*Phi2_her 0];
Part1=(norm(G(1,1))^2+norm(G(2,1))^2)*[Phi1_her'*Phi1_her Phi1_her'*H3(1,1); H3(1,1)'*Phi1_her 0];
Part2=(norm(G(1,1))^2+norm(G(2,1))^2)*[Phi2_her'*Phi2_her Phi2_her'*H3(1,2); H3(1,2)'*Phi2_her 0];
Part3=(G(1,2)'*G(1,1)+G(2,2)'*G(2,1))*A;
Part4=(norm(G(1,2))^2+norm(G(2,2))^2)*[Phi3_her'*Phi3_her Phi3_her'*H3(2,1); H3(2,1)'*Phi3_her 0];
Part5=(norm(G(1,2))^2+norm(G(2,2))^2)*[Phi4_her'*Phi4_her Phi4_her'*H3(2,2); H3(2,2)'*Phi4_her 0];
Part6=[zeros(M,M) (Phi1_her'+Phi2_her')*G(1,1)'; (Phi1_her+Phi2_her)*G(1,1)  0];
Part7=[zeros(M,M) (Phi1_her'-Phi2_her')*G(2,1)'; (Phi1_her-Phi2_her)*G(2,1)  0];
Part8=[zeros(M,M) (Phi3_her'+Phi4_her')*G(1,2)'; (Phi3_her+Phi4_her)*G(1,2)  0];
Part9=[zeros(M,M) (Phi3_her'-Phi4_her')*G(2,2)'; (Phi3_her-Phi4_her)*G(2,2)  0];
M_matrix=Part1+Part2+Part3+Part4+Part5-Part6-Part7-Part8-Part9;


% V_final=[v_phases; 1]*[v_phases' 1];
% V_hat=zeros(M+1,M+1);
% inner_iteration=1;
% while(((norm(V_final-V_hat,'fro')^2)>=error)&&(inner_iteration<=max_iteration))
% 
% V_hat=V_final;
% [Eigen_vectors,D_eigenvalues] = eigs(V_hat);
% [maxValue,index_max]=max(diag(D_eigenvalues));  %# The maximum eigenvalue and its index
% u1=Eigen_vectors(:,index_max);           %# The associated eigenvector in V    
%     
    
cvx_setup
cvx_begin 
    variable V(M+1,M+1) complex;
    minimize(real(trace(M_matrix*V)));
    subject to
    %V>=0;
    V==hermitian_semidefinite(M+1);
    eye(M+1)*diag(V)==ones(M+1,1);
cvx_end

V_final=V;
% 
% inner_iteration=inner_iteration+1;
% 
% end
if(rank(V_final)>1)
    [U,E]=eig(V_final);
mu=zeros(1,M+1);
Varint=Inf;
Icov=eye(M+1);
R=mvnrnd(mu,Icov,Ngau);
R=conj(R');
for l=1:Ngau
   v_iter=U*sqrt(E)*R(:,l);
   Var=v_iter'*M_matrix*v_iter; 
    if(Var<Varint)
        Varint=Var;
        Lmin=l;
    end
end
r=U*sqrt(E)*R(:,Lmin);
rnormalized=r/r(M+1);
v_next=rnormalized(1:M);
else
   v_next=V_final(M+1,1:M)';
end
ph=angle(v_next);
MSE_iter_total(ii,iteration+1)=MSEfunction(ph,G,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);

MSEdecreasepercent=norm(MSE_iter_total(ii,iteration+1)-MSE_previous);
MSE_previous=MSE_iter_total;
% 
% MSE1_total2=Objective_function_given_G_CVX(v_phases,G2,Phi1_her,Phi2_her,Phi3_her,Phi4_her,P,H3,sigma2(jj),M,Nr);
% MSEdecreasepercent2=abs(((MSE1_total2-MSE_previous2)/MSE_previous2)*100);
% MSE_previous2=MSE1_total2;
iteration=iteration+1;

end
tEnd=toc(tStart);
theta=ph;

Theta=diag(exp(j*theta));
H=Hre_irs*Theta*Hirs_u+Hre_u;
H_hat=H*inv(D);
% G=P*D*H_hat'*inv(H_hat*H_hat'+sigma2(jj)*I);
Check_identity(ii,jj)=trace(G*H_hat)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%  Case 1 .  Theoretical Optimization of the 
%%%%%%%%%%%%%%%%%%%%%%%%%%%  phases given G Proposed Algorithm
tStart1=tic;
MSE_previous1=0;
MSEdecreasepercent1=Inf;
ph1=pi*ones(1,M);
v_phases=exp(1i*ph1');
iteration1=0;
while((MSEdecreasepercent1>=error)&&(iteration1<=max_iteration))

theta1=ph1;
Theta1=diag(exp(1i*theta1));
H_total1=Hre_irs*Theta1*Hirs_u+Hre_u;
% H_hat1=H_total1*inv(D);
% G3=inv(sigma2(jj)*II+H_hat2'*H_hat2)*H_hat2';
G1=P*D*H_total1'*inv(H_total1*H_total1'+sigma2(jj)*I);

A=[Phi3_her'*Phi1_her+Phi2_her'*Phi4_her+Phi1_her'*Phi3_her+Phi4_her'*Phi2_her Phi3_her'*H3(1,1)+Phi2_her'*H3(2,2)+Phi1_her'*H3(2,1)+Phi4_her'*H3(1,2); H3(2,1)'*Phi1_her+H3(1,2)'*Phi4_her+H3(1,1)'*Phi3_her+H3(2,2)'*Phi2_her 0];
Part1=(norm(G1(1,1))^2+norm(G1(2,1))^2)*[Phi1_her'*Phi1_her Phi1_her'*H3(1,1); H3(1,1)'*Phi1_her 0];
Part2=(norm(G1(1,1))^2+norm(G1(2,1))^2)*[Phi2_her'*Phi2_her Phi2_her'*H3(1,2); H3(1,2)'*Phi2_her 0];
Part3=(G1(1,2)'*G1(1,1)+G1(2,2)'*G1(2,1))*A;
Part4=(norm(G1(1,2))^2+norm(G1(2,2))^2)*[Phi3_her'*Phi3_her Phi3_her'*H3(2,1); H3(2,1)'*Phi3_her 0];
Part5=(norm(G1(1,2))^2+norm(G1(2,2))^2)*[Phi4_her'*Phi4_her Phi4_her'*H3(2,2); H3(2,2)'*Phi4_her 0];
Part6=[zeros(M,M) (Phi1_her'+Phi2_her')*G1(1,1)'; (Phi1_her+Phi2_her)*G1(1,1)  0];
Part7=[zeros(M,M) (Phi1_her'-Phi2_her')*G1(2,1)'; (Phi1_her-Phi2_her)*G1(2,1)  0];
Part8=[zeros(M,M) (Phi3_her'+Phi4_her')*G1(1,2)'; (Phi3_her+Phi4_her)*G1(1,2)  0];
Part9=[zeros(M,M) (Phi3_her'-Phi4_her')*G1(2,2)'; (Phi3_her-Phi4_her)*G1(2,2)  0];
M_matrix=Part1+Part2+Part3+Part4+Part5-Part6-Part7-Part8-Part9;


V_final=[v_phases; 1]*[v_phases' 1];
V_hat=zeros(M+1,M+1);
inner_iteration=1;
while(((norm(V_final-V_hat,'fro')^2)>=error)&&(inner_iteration<=max_iteration))

V_hat=V_final;
[Eigen_vectors,D_eigenvalues] = eigs(V_hat);
[maxValue,index_max]=max(diag(D_eigenvalues));  %# The maximum eigenvalue and its index
u1=Eigen_vectors(:,index_max);           %# The associated eigenvector in V    
    
    
cvx_setup
cvx_begin 
    variable V(M+1,M+1) complex;
    minimize(real(trace(M_matrix*V)+gamma2*(trace(V)-trace(V*(u1)*(u1)'))));
    subject to
    %V>=0;
    V==hermitian_semidefinite(M+1);
    eye(M+1)*diag(V)==ones(M+1,1);
cvx_end

V_final=V;

inner_iteration=inner_iteration+1;

end

v_next=V_final(M+1,1:M)';
% v_next=v_next./abs(v_next);
v_phases=v_next;
ph1=angle(v_next);

MSE_iter_total1(ii,iteration1+1)=MSEfunction(ph1,G1,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);

MSEdecreasepercent1=norm(MSE_iter_total1(ii,iteration1+1)-MSE_previous1);
MSE_previous1=MSE_iter_total1;
% 
% MSE1_total2=Objective_function_given_G_CVX(v_phases,G2,Phi1_her,Phi2_her,Phi3_her,Phi4_her,P,H3,sigma2(jj),M,Nr);
% MSEdecreasepercent2=abs(((MSE1_total2-MSE_previous2)/MSE_previous2)*100);
% MSE_previous2=MSE1_total2;
iteration1=iteration1+1;

end

tEnd1=toc(tStart1);
theta1=ph1;
Theta1=diag(exp(1i*theta1));
H_total1=Hre_irs*Theta1*Hirs_u+Hre_u;
H_hat1=H_total1*inv(D);

Check_identity1(ii,jj)=trace(G1*H_hat1)/2;

theta_q=(((1*pi/2)<theta1).*(theta1<(3*pi/2)))*pi;

Theta_q=diag(exp(j*theta_q));
H_q=Hre_irs*Theta_q*Hirs_u+Hre_u;
H_q_hat=H_q*inv(D);

G_q=P*D*H_q'*inv(H_q*H_q'+sigma2(jj)*I);

Check_identity_q(ii,jj)=trace(G_q*H_q_hat)/2;






% %%%%%%%Case 1
% theta1=zeros(1,M); %zeros phases
% Theta1=diag(exp(j*theta1));
% H1=Hre_irs*Theta1*Hirs_u+Hre_u;
% H_hat1=H1*inv(D);
% 
% G1=P*D*H1'*inv(H1*H1'+sigma2(jj)*I);

%%%%%%%Case 2 Random phases 
Theta2=diag(exp(j*theta2));
H2=Hre_irs*Theta2*Hirs_u+Hre_u;
H_hat2=H2*inv(D);

G2=P*D*H2'*inv(H2*H2'+sigma2(jj)*I);

Check_identity2(ii,jj)=trace(G2*H_hat2)/2;




%%%%%%%Case optimizing G and theta iteratively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart3=tic;
iteration3=0;
MSE_previous3=Inf;
MSEdecreasepercent3=Inf;
ph3=(pi/2)*ones(1,M);
while((MSEdecreasepercent3>error)&&(iteration3<=max_iteration))

theta3=ph3;
Theta_iter=diag(exp(j*theta3));
H_iter=Hre_irs*Theta_iter*Hirs_u+Hre_u;
% H_hat3=H_iter*inv(D);
% G=inv(H_hat'*H_hat+sigma2(jj)*II)*H_hat';
G3=P*D*H_iter'*inv(H_iter*H_iter'+sigma2(jj)*I);
% t_iter_sol=[real(G_iter(1,1:Nr)) imag(G_iter(1,1:Nr)) real(G_iter(2,1:Nr)) imag(G_iter(2,1:Nr))];
% MSE_iter_givenTheta=@(t) MSEfunction_G(theta,t,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
% t_0=ones(1,2*Ns*Nr);
% A_t=[];
% b_t=[];
% t_iter_sol=fmincon(MSE_iter_givenTheta,t_0,A_t,b_t,[],[],[],[],[],options);
% G=[t_iter_sol(1:Nr)+1i*t_iter_sol(Nr+1:2*Nr); t_iter_sol(2*Nr+1:3*Nr)+1i*t_iter_sol(3*Nr+1:4*Nr)];


MSE_givenG=@(ph) MSEfunction(ph3(1:M),G3,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);

lb=zeros(1,M);
ub=2*pi*ones(1,M);
ph0=(lb+ub)/2;
ph3=fmincon(MSE_givenG,ph0,[],[],[],[],lb,ub,[],options);
MSE_iter_total3(ii,iteration3+1)=MSEfunction(ph3,G3,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);

MSEdecreasepercent3=norm(MSE_iter_total3(ii,iteration3+1)-MSE_previous3);
MSE_previous3=MSE_iter_total3;
iteration3=iteration3+1;
end
% 
tEnd3=toc(tStart3);
theta3=ph3;

Theta3=diag(exp(j*theta3));
H_total3=Hre_irs*Theta3*Hirs_u+Hre_u;
H_hat3=H_total3*inv(D);
% G=P*D*H_hat'*inv(H_hat*H_hat'+sigma2(jj)*I);
Check_identity3(ii,jj)=trace(G3*H_hat3)/2;


%%%%%%%Case 4
H4=Hre_u;
H_hat4=H4*inv(D);

% G4=inv(H_hat4'*H_hat4+sigma2(jj)*II)*H_hat4';
G4=P*D*H4'*inv(H4*H4'+sigma2(jj)*I);

Check_identity4(ii,jj)=trace(G4*H_hat4)/2;
%%%%%%%Case 5

G5=P*Hre_u'*inv(Hre_u*Hre_u'+sigma2(jj)*I);

Check_identity5(ii,jj)=trace(G5*Hre_u)/2;
% +sigma2(jj)*II

%%%MSE plotting
Final_MSE(ii,jj)=MSEfunction(theta,G,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
Final_MSE_q(ii,jj)=MSEfunction(theta_q,G_q,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
Final_MSE1(ii,jj)=MSEfunction(theta1,G1,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
Final_MSE2(ii,jj)=MSEfunction(theta2,G2,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
Final_MSE3(ii,jj)=MSEfunction(theta3,G3,Nr,M,P,Hre_irs,Hirs_u,Hre_u,sigma2(jj),Ns,D);
Final_MSE4(ii,jj)=MSEfunction_NO_IRS(G4,Nr,P,Hre_u,sigma2(jj),Ns,D);
Final_MSE5(ii,jj)=MSEfunction_NO_IRS_NO_PNC(G5,Nr,P,Hre_u,sigma2(jj),Ns);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noise between sources and the relay
Noise=sqrt(sigma2(jj))*(randn(Nr,Length)+1i*randn(Nr,Length));
%Noise between S1 and D1
Noise_s1d1=sqrt(sigma2(jj))*(randn(1,Length)+1i*randn(1,Length));

%Noise between relay and d1
Noise_D1=sqrt(sigma2(jj))*(randn(1,Length)+1i*randn(1,Length));


X_hat=D*X; %input is common for all cases




%received signals for each case

R=H_hat*X_hat+Noise;
RSNR(jj)=10*log10((norm(H_hat*X_hat)^2)/(sigma2(jj)*Length));
y=G*R;

R_q=H_q_hat*X_hat+Noise;
RSNR_q(jj)=10*log10((norm(H_q_hat*X_hat)^2)/(sigma2(jj)*Length));
y_q=G_q*R_q;


R1=H_hat1*X_hat+Noise;
RSNR1(jj)=10*log10((norm(H_hat1*X_hat)^2)/(sigma2(jj)*Length));
y1=G1*R1;


R2=H_hat2*X_hat+Noise;
RSNR2(jj)=10*log10((norm(H_hat2*X_hat)^2)/(sigma2(jj)*Length));
y2=G2*R2;

R3=H_hat3*X_hat+Noise;
RSNR3(jj)=10*log10((norm(H_hat3*X_hat)^2)/(sigma2(jj)*Length));
y3=G3*R3;

R4=H_hat4*X_hat+Noise;
RSNR4(jj)=10*log10((norm(H_hat4*X_hat)^2)/(sigma2(jj)*Length));
y4=G4*R4;

R5=Hre_u*X+Noise;
RSNR5(jj)=10*log10((norm(Hre_u*X)^2)/(sigma2(jj)*Length));
y5=G5*R5;



Rs1d1=Hs1d1*X(1,:)+Noise_s1d1;
RSNR5_s1d1(jj)=10*log10((norm(Hs1d1*X(1,:))^2)/(sigma2(jj)*Length));
% Gs1d1=inv(Hs1d1'*Hs1d1+sigma2(jj))*Hs1d1';
Gs1d1=inv(Hs1d1);
ys1d1=Gs1d1*Rs1d1;

sigma2_s1d1=Gs1d1'*Gs1d1*sigma2(jj);

% G's and sigmas for each case

G_G=G*G';

sigma12_check(ii,jj)=G_G(1,1)*sigma2(jj);
sigma22_check(ii,jj)=G_G(2,2)*sigma2(jj);
sigma12=sigma12_check(ii,jj);
sigma22=sigma22_check(ii,jj);

G_G_q=G_q*G_q';

sigma12_q=G_G_q(1,1)*sigma2(jj);
sigma22_q=G_G_q(2,2)*sigma2(jj);

G_G1=G1*(G1)';

sigma12_1=G_G1(1,1)*sigma2(jj);
sigma22_1=G_G1(2,2)*sigma2(jj);


G_G2=G2*(G2)';

sigma12_2=G_G2(1,1)*sigma2(jj);
sigma22_2=G_G2(2,2)*sigma2(jj);

G_G3=G3*(G3)';

sigma12_3=G_G3(1,1)*sigma2(jj);
sigma22_3=G_G3(2,2)*sigma2(jj);


G_G4=G4*G4';

sigma12_4=G_G4(1,1)*sigma2(jj);
sigma22_4=G_G4(2,2)*sigma2(jj);


G_G5=G5*G5';

sigma12_5=G_G5(1,1)*sigma2(jj);
sigma22_5=G_G5(2,2)*sigma2(jj);

% Final_MSE(ii,jj)=norm(y-X_hat)/length(y);
% Final_MSE_q(ii,jj)=norm(y_q-X_hat)/length(y);
% Final_MSE1(ii,jj)=norm(y1-X_hat)/length(y);
% Final_MSE2(ii,jj)=norm(y2-X_hat)/length(y);
% Final_MSE4(ii,jj)=norm(y4-X_hat)/length(y);
% Final_MSE5(ii,jj)=norm(y5-X)/length(y);

% L=log((exp((2/sigma12)-(2/sigma22))+lN).*((cosh((2.*real(y(2,:)))./sigma22)+lN)./(cosh((2.*real(y(1,:)))./sigma12)+lN)+lN));
L=2*((1/sigma12)-(1/sigma22))+log(cosh((2.*real(y(2,:)))./sigma22)+lN)-log(cosh((2.*real(y(1,:)))./sigma12)+lN);
L_q=2*((1/sigma12_q)-(1/sigma22_q))+log(cosh((2.*real(y_q(2,:)))./sigma22_q)+lN)-log(cosh((2.*real(y_q(1,:)))./sigma12_q)+lN);
L1=2*((1/sigma12_1)-(1/sigma22_1))+log(cosh((2.*real(y1(2,:)))./sigma22_1)+lN)-log(cosh((2.*real(y1(1,:)))./sigma12_1)+lN);
L2=2*((1/sigma12_2)-(1/sigma22_2))+log(cosh((2.*real(y2(2,:)))./sigma22_2)+lN)-log(cosh((2.*real(y2(1,:)))./sigma12_2)+lN);
L3=2*((1/sigma12_3)-(1/sigma22_3))+log(cosh((2.*real(y3(2,:)))./sigma22_3)+lN)-log(cosh((2.*real(y3(1,:)))./sigma12_3)+lN);
L4=2*((1/sigma12_4)-(1/sigma22_4))+log(cosh((2.*real(y4(2,:)))./sigma22_4)+lN)-log(cosh((2.*real(y4(1,:)))./sigma12_4)+lN);
Ls1d1=(2.*real(ys1d1))./sigma2_s1d1;
% Ls1d1=exp(((2.*real(Rs1d1*Hs1d1))./(sigma2(jj)))+lN);



%%%%%% data transmitted from the relay
x1_xor_x2_transmitted=(L>=0)*1+(L<0)*(-1); 

x1_xor_x2_transmitted_q=(L_q>=0)*1+(L_q<0)*(-1); 

x1_xor_x2_transmitted1=(L1>=0)*1+(L1<0)*(-1); %+isnan(L1)*t1;

x1_xor_x2_transmitted2=(L2>=0)*1+(L2<0)*(-1); %+isnan(L2)*t2;

x1_xor_x2_transmitted3=(L3>=0)*1+(L3<0)*(-1); %+isnan(L2)*t2;

x1_xor_x2_transmitted4=(L4>=0)*1+(L4<0)*(-1); %+isnan(L4)*t4;

%%%%%%%%%%%%%% estimation of x1 and x2 separately at the relay
L5_x1=(2.*real(y5(1,:)))./sigma12_5;
L5_x2=(2.*real(y5(2,:)))./sigma22_5;
x1_NNC_rel=(L5_x1>=0)*1+(L5_x1<0)*(-1);
x2_NNC_rel=(L5_x2>=0)*1+(L5_x2<0)*(-1);
x1_xor_x2_transmitted5=(-1).*(x1_NNC_rel.*x2_NNC_rel);

%%% Calculating BER at the relay
x1_xor_x2_transmitted_bits=(x1_xor_x2_transmitted==1)*1;
x1_xor_x2_transmitted_bits_q=(x1_xor_x2_transmitted_q==1)*1;
x1_xor_x2_transmitted_bits1=(x1_xor_x2_transmitted1==1)*1;
x1_xor_x2_transmitted_bits2=(x1_xor_x2_transmitted2==1)*1;
x1_xor_x2_transmitted_bits3=(x1_xor_x2_transmitted3==1)*1;
x1_xor_x2_transmitted_bits4=(x1_xor_x2_transmitted4==1)*1;
BER_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER_q_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits_q-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER1_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits1-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER2_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits2-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER3_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits3-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER4_SR(ii,jj)=sum((x1_xor_x2_transmitted_bits4-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
% BER_SR_th(ii,jj)=(1/2)*qfunc(1/sqrt(sigma12))+(1/2)*(qfunc(2/sqrt(sigma12))*qfunc(-3/sqrt(sigma12))+qfunc(-2/sqrt(sigma12))*qfunc(1/sqrt(sigma12)))+(1/2)*qfunc(1/sqrt(sigma22))+(1/2)*(qfunc(2/sqrt(sigma22))*qfunc(-3/sqrt(sigma22))+qfunc(-2/sqrt(sigma22))*qfunc(1/sqrt(sigma22)));
% BER_SR_th(ii,jj)=(1/2)*qfunc(1/sqrt(sigma22))+(1/2)*(qfunc(2/sqrt(sigma22))*qfunc(1/sqrt(sigma22))+qfunc(-2/sqrt(sigma22))*qfunc(-3/sqrt(sigma22)));
% BER_SR_th(ii,jj)=(7/2)*((qfunc(-2/sqrt(sigma12))+qfunc(-2/sqrt(sigma22)))*qfunc(sqrt(1+sigma22/sigma12)/sqrt(sigma22))+qfunc(2/sqrt(sigma22))*qfunc((-3+sigma22/sigma12)/(sqrt(sigma22)*sqrt(1+sigma22/sigma12)))+qfunc(2/sqrt(sigma12))*qfunc((-3+sigma12/sigma22)/(sqrt(sigma12)*sqrt(1+sigma12/sigma22))));
%BER_SR_th(ii,jj)=(1/2)*qfunc(2/sqrt(sigma22))*qfunc(((2/sigma12)+(2/sigma22))/sqrt((4/sigma12)+(4/sigma22)))+(1/2)*qfunc(-2/sqrt(sigma22))*qfunc(((2/sigma12)-(6/sigma22))/sqrt((4/sigma12)+(4/sigma22)))+(1/2)*qfunc(2/sqrt(sigma12))*qfunc(((2/sigma22)+(2/sigma12))/sqrt((4/sigma12)+(4/sigma22)))+(1/2)*qfunc(-2/sqrt(sigma12))*qfunc(((2/sigma22)-(6/sigma12))/sqrt((4/sigma12)+(4/sigma22)));
% BER_SR_th(ii,jj)=(1/2)*qfunc(2/sqrt(sigma22))*qfunc(((2/sigma12)+(2/sigma22))/sqrt((4/sigma12)+(4/sigma22)))+(1/2)*qfunc(-2/sqrt(sigma22))*qfunc(((2/sigma12)-(6/sigma22))/sqrt((4/sigma12)+(4/sigma22)))+(1/2)*qfunc(2/sqrt(sigma12))*qfunc(((2/sigma22)+(2/sigma12))/sqrt((4/sigma12)+(4/sigma22)));
BER_SR_th(ii,jj)=(7/2)*(qfunc(-2/sqrt(sigma12))+qfunc(-2/sqrt(sigma22)))*qfunc(sqrt(1+sigma22/sigma12)/sqrt(sigma22))+qfunc(2/sqrt(sigma22))*qfunc((-3+sigma22/sigma12)/(sqrt(sigma22)*sqrt(1+sigma22/sigma12)))+qfunc(2/sqrt(sigma12))*qfunc((-3+sigma12/sigma22)/(sqrt(sigma12)*sqrt(1+sigma12/sigma22)));

% BER_SR_th(ii,jj)=qfunc(2/sqrt(sigma22))*qfunc(((1/sigma12)-(3/sigma22))/sqrt((1/sigma12)+(1/sigma22)))+qfunc(-2/sqrt(sigma22))*qfunc(((1/sigma12)+(1/sigma22))/sqrt((1/sigma12)+(1/sigma22)))+qfunc(2/sqrt(sigma12))*qfunc(((1/sigma22)-(3/sigma12))/sqrt((1/sigma12)+(1/sigma22)))+qfunc(-2/sqrt(sigma12))*qfunc(((1/sigma22)+(1/sigma12))/sqrt((1/sigma12)+(1/sigma22)));
BER_q_SR_th(ii,jj)=(7/2)*(qfunc(-2/sqrt(sigma12_q))+qfunc(-2/sqrt(sigma22_q)))*qfunc(sqrt(1+sigma22_q/sigma12)/sqrt(sigma22_q))+qfunc(2/sqrt(sigma22_q))*qfunc((-3+sigma22_q/sigma12_q)/(sqrt(sigma22_q)*sqrt(1+sigma22_q/sigma12_q)))+qfunc(2/sqrt(sigma12_q))*qfunc((-3+sigma12_q/sigma22_q)/(sqrt(sigma12_q)*sqrt(1+sigma12_q/sigma22_q)));
BER1_SR_th(ii,jj)=(7/2)*(qfunc(-2/sqrt(sigma12_1))+qfunc(-2/sqrt(sigma22_1)))*qfunc(sqrt(1+sigma22_1/sigma12_1)/sqrt(sigma22_1))+qfunc(2/sqrt(sigma22_1))*qfunc((-3+sigma22_1/sigma12_1)/(sqrt(sigma22_1)*sqrt(1+sigma22_1/sigma12_1)))+qfunc(2/sqrt(sigma12_1))*qfunc((-3+sigma12_1/sigma22_1)/(sqrt(sigma12_1)*sqrt(1+sigma12_1/sigma22_1)));
BER2_SR_th(ii,jj)=(7/2)*(qfunc(-2/sqrt(sigma12_2))+qfunc(-2/sqrt(sigma22_2)))*qfunc(sqrt(1+sigma22_2/sigma12_2)/sqrt(sigma22_2))+qfunc(2/sqrt(sigma22_2))*qfunc((-3+sigma22_2/sigma12_2)/(sqrt(sigma22_2)*sqrt(1+sigma22_2/sigma12_2)))+qfunc(2/sqrt(sigma12_2))*qfunc((-3+sigma12_2/sigma22_2)/(sqrt(sigma12_2)*sqrt(1+sigma12_2/sigma22_2)));
BER3_SR_th(ii,jj)=(7/2)*(qfunc(-2/sqrt(sigma12_3))+qfunc(-2/sqrt(sigma22_3)))*qfunc(sqrt(1+sigma22_3/sigma12_3)/sqrt(sigma22_3))+qfunc(2/sqrt(sigma22_3))*qfunc((-3+sigma22_3/sigma12_3)/(sqrt(sigma22_3)*sqrt(1+sigma22_3/sigma12_3)))+qfunc(2/sqrt(sigma12_3))*qfunc((-3+sigma12_3/sigma22_3)/(sqrt(sigma12_3)*sqrt(1+sigma12_3/sigma22_3)));
BER4_SR_th(ii,jj)=(1/2)*(qfunc(-2/sqrt(sigma12_4))+qfunc(-2/sqrt(sigma22_4)))*qfunc(sqrt(1+sigma22_4/sigma12_4)/sqrt(sigma22_4))+qfunc(2/sqrt(sigma22_4))*qfunc((-3+sigma22_4/sigma12_4)/(sqrt(sigma22_4)*sqrt(1+sigma22_4/sigma12_4)))+qfunc(2/sqrt(sigma12_4))*qfunc((-3+sigma12_4/sigma22_4)/(sqrt(sigma12_4)*sqrt(1+sigma12_4/sigma22_4)));


x1_xor_x2_transmitted5_bits=(x1_xor_x2_transmitted5==1)*1;
BER5_SR(ii,jj)=sum((x1_xor_x2_transmitted5_bits-x1_xor_x2_transmitted_actual_bits)~=0)/(Length);
BER5_SR_th(ii,jj)=qfunc(1/sqrt(sigma12_5))+qfunc(1/sqrt(sigma22_5))-2*qfunc(1/sqrt(sigma12_5))*qfunc(1/sqrt(sigma22_5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%beamforming for all cases
m=1/norm(Hr_D1);
W=Hr_D1'./norm(Hr_D1);



%received signal at D1 for each case
Us_rec=Hr_D1*W*x1_xor_x2_transmitted+Noise_D1;
USNR(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted)^2)/(sigma2(jj)*Length));
Output=m*Us_rec;

Us_rec_q=Hr_D1*W*x1_xor_x2_transmitted_q+Noise_D1;
USNR_q(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted_q)^2)/(sigma2(jj)*Length));
Output_q=m*Us_rec_q;

Us_rec1=Hr_D1*W*x1_xor_x2_transmitted1+Noise_D1;
USNR1(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted1)^2)/(sigma2(jj)*Length));
Output1=m*Us_rec1;

Us_rec2=Hr_D1*W*x1_xor_x2_transmitted2+Noise_D1;
USNR2(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted2)^2)/(sigma2(jj)*Length));
Output2=m*Us_rec2;

Us_rec3=Hr_D1*W*x1_xor_x2_transmitted3+Noise_D1;
USNR3(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted3)^2)/(sigma2(jj)*Length));
Output3=m*Us_rec3;

Us_rec4=Hr_D1*W*x1_xor_x2_transmitted4+Noise_D1;
USNR4(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted4)^2)/(sigma2(jj)*Length));
Output4=m*Us_rec4;


Us_rec5=Hr_D1*W*x1_xor_x2_transmitted5+Noise_D1;
USNR5(jj)=10*log10((norm(Hr_D1*W*x1_xor_x2_transmitted5)^2)/(sigma2(jj)*Length));
Output5=m*Us_rec5;




%%%%%%%%%%%%%%%%%%%% LLRs at D1
% sigma_output=sigma2(jj)/(norm(Hr_D1)^2);
sigma2_output=(m^2)*sigma2(jj);
Lout=(2.*real(Output))./sigma2_output;
Lout_q=(2.*real(Output_q))./sigma2_output;
Lout1=(2.*real(Output1))./sigma2_output;
Lout2=(2.*real(Output2))./sigma2_output;
Lout3=(2.*real(Output3))./sigma2_output;
Lout4=(2.*real(Output4))./sigma2_output;
Lout5=(2.*real(Output5))./sigma2_output;


%received signals at D1 from relay
x1_xor_x2_received=(Lout>=0)*1+(Lout<0)*(-1); 
x1_xor_x2_received_q=(Lout_q>=0)*1+(Lout_q<0)*(-1); 
x1_xor_x2_received1=(Lout1>=0)*1+(Lout1<0)*(-1); 
x1_xor_x2_received2=(Lout2>=0)*1+(Lout2<0)*(-1); 
x1_xor_x2_received3=(Lout3>=0)*1+(Lout3<0)*(-1);
x1_xor_x2_received4=(Lout4>=0)*1+(Lout4<0)*(-1); 
x1_xor_x2_received5=(Lout5>=0)*1+(Lout5<0)*(-1); 

%received signals at D1 from S1
x1_received_s1d1=(Ls1d1>=0)*1+(Ls1d1<0)*(-1); 

%%%%% x2 received at D1
x2_received=(-1)*x1_received_s1d1.*x1_xor_x2_received;
x2_received_q=(-1)*x1_received_s1d1.*x1_xor_x2_received_q;
x2_received2=(-1)*x1_received_s1d1.*x1_xor_x2_received2;
x2_received3=(-1)*x1_received_s1d1.*x1_xor_x2_received3;
x2_received4=(-1)*x1_received_s1d1.*x1_xor_x2_received4;
x2_received1=(-1)*x1_received_s1d1.*x1_xor_x2_received1;
x2_received5=(-1)*x1_received_s1d1.*x1_xor_x2_received5;

%%%%% x2 received at D1 in bits
x2_received_bits=(x2_received==1)*1;
x2_received_bits_q=(x2_received_q==1)*1;
x2_received_bits1=(x2_received1==1)*1;
x2_received_bits2=(x2_received2==1)*1;
x2_received_bits3=(x2_received3==1)*1;
x2_received_bits4=(x2_received4==1)*1;
x2_received_bits5=(x2_received5==1)*1;

x1_received_s1d1_bits=(x1_received_s1d1==1)*1;


%BER at D1 from S1 x1
BER_s1d1(ii,jj)=sum((x1_received_s1d1_bits-x1_actual)~=0)/(Length);
BER_s1d1_th(ii,jj)=qfunc(1/sqrt(sigma2_s1d1));
% BER_s1d1_th(ii,jj)=qfunc(real(Hs1d1)/sqrt(sigma2(jj)));

%BER from relay to D1
x1_xor_x2_received_bits=(x1_xor_x2_received==1)*1;
BER_rd1(ii,jj)=sum((x1_xor_x2_received_bits-x1_xor_x2_transmitted_bits)~=0)/(Length);
BER_rd1_th(ii,jj)=qfunc(1/sqrt(sigma2_output));

%BER at D1 for x2
BER(ii,jj)=sum((x2_received_bits-x2_actual)~=0)/(Length);
BER_q(ii,jj)=sum((x2_received_bits_q-x2_actual)~=0)/(Length);
BER1(ii,jj)=sum((x2_received_bits1-x2_actual)~=0)/(Length);
BER2(ii,jj)=sum((x2_received_bits2-x2_actual)~=0)/(Length);
BER3(ii,jj)=sum((x2_received_bits3-x2_actual)~=0)/(Length);
BER4(ii,jj)=sum((x2_received_bits4-x2_actual)~=0)/(Length);
BER5(ii,jj)=sum((x2_received_bits5-x2_actual)~=0)/(Length);


BER5_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER5_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER5_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER5_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER5_SR_th(ii,jj);
BER_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER_SR_th(ii,jj);
BER_q_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER_q_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER_q_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER_q_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER_q_SR_th(ii,jj);
BER1_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER1_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER1_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER1_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER1_SR_th(ii,jj);
BER2_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER2_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER2_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER2_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER2_SR_th(ii,jj);
BER3_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER3_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER3_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER3_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER3_SR_th(ii,jj);
BER4_th(ii,jj)=BER_s1d1_th(ii,jj)*(1-BER_rd1_th(ii,jj))*(1-BER4_SR_th(ii,jj))+BER_s1d1_th(ii,jj)*BER_rd1_th(ii,jj)*BER4_SR_th(ii,jj)+(1-BER_s1d1_th(ii,jj))*BER_rd1_th(ii,jj)*(1-BER4_SR_th(ii,jj))+(1-BER_s1d1_th(ii,jj))*(1-BER_rd1_th(ii,jj))*BER4_SR_th(ii,jj);

end
ii
end
%%
%%%At D1
BER_er=mean(BER);
BER_th_er=mean(BER_th);
BER_q_er=mean(BER_q);
BER_q_th_er=mean(BER_q_th);
BER1_er=mean(BER1);
BER1_th_er=mean(BER1_th);
BER2_er=mean(BER2);
BER2_th_er=mean(BER2_th);
BER3_er=mean(BER3);
BER3_th_er=mean(BER3_th);
BER4_er=mean(BER4);
BER4_th_er=mean(BER4_th);
BER5_er=mean(BER5);
BER5_th_er=mean(BER5_th);

BER_s1d1_er=mean(BER_s1d1);
BER_s1d1_th_er=mean(BER_s1d1_th);

BER_rd1_er=mean(BER_rd1);
BER_rd1_th_er=mean(BER_rd1_th);

%%% at relay
BER_SR_er=mean(BER_SR);
BER_SR_th_er=mean(BER_SR_th);
BER_q_SR_er=mean(BER_q_SR);
BER_q_SR_th_er=mean(BER_q_SR_th);
BER1_SR_er=mean(BER1_SR);
BER1_SR_th_er=mean(BER1_SR_th);
BER2_SR_er=mean(BER2_SR);
BER2_SR_th_er=mean(BER2_SR_th);
BER3_SR_er=mean(BER3_SR);
BER3_SR_th_er=mean(BER3_SR_th);
BER4_SR_er=mean(BER4_SR);
BER4_SR_th_er=mean(BER4_SR_th);
BER5_SR_er=mean(BER5_SR);
BER5_SR_th_er=mean(BER5_SR_th);


%%% MSE ergodic
Final_MSE_er=mean(Final_MSE);
Final_MSE_q_er=mean(Final_MSE_q);
Final_MSE1_er=mean(Final_MSE1);
Final_MSE2_er=mean(Final_MSE2);
Final_MSE3_er=mean(Final_MSE3);
Final_MSE4_er=mean(Final_MSE4);
Final_MSE5_er=mean(Final_MSE5);




Check_identity_er=mean(Check_identity);
Check_identity_q_er=mean(Check_identity_q);
Check_identity1_er=mean(Check_identity1);
Check_identity2_er=mean(Check_identity2);
Check_identity3_er=mean(Check_identity3);
Check_identity4_er=mean(Check_identity4);
Check_identity5_er=mean(Check_identity5);
%%
%   cd '/Users/kafizoa/Dropbox/Amanat/Amanat_Course_Project/Paper_figures/BER_vs_SNR_relay_MMSE'
%   save('BER_vs_SNR_relay_workspace','SNR_dB','BER_SR_er','BER_SR_th_er','BER_q_SR_er','BER_q_SR_th_er','BER2_SR_er','BER2_SR_th_er','BER4_SR_er','BER4_SR_th_er','BER5_SR_er','BER5_SR_th_er');
%%
%  cd '/Users/kafizoa/Dropbox/Amanat/Amanat_Course_Project/Paper_figures/MSE_vs_SNR_relay_MMSE'
% save('MSE_vs_SNR_relay_workspace','Final_MSE_er','Final_MSE_q_er','Final_MSE2_er','Final_MSE4_er','Final_MSE5_er')
%%
%  cd '/Users/kafizoa/Dropbox/Amanat/Amanat_Course_Project/Paper_figures/BER_vs_SNR_D1_MMSE'
%  save('BER_vs_SNR_D1_workspace','SNR_dB','BER_er','BER_th_er','BER_q_er','BER_q_th_er','BER2_er','BER2_th_er','BER4_er','BER4_th_er','BER5_er','BER5_th_er');

%%
figure(1)
semilogy(SNR_dB,BER_SR_er,'-.r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("BER",'FontSize', 25);
title("BER at relay P(x1xorx2hat!=x1xorx2) vs SNR",'FontSize', 25);
set(gca,'FontSize',25)
ylim([10e-5 1])
hold on;
semilogy(SNR_dB,BER_SR_th_er,'-r','LineWidth',5);
hold on;
semilogy(SNR_dB,BER_q_SR_er,'--c+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER_q_SR_th_er,'-c+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER1_SR_er,':bs','LineWidth',5);
hold on;
semilogy(SNR_dB,BER1_SR_th_er,'-bs','LineWidth',5);
hold on;
semilogy(SNR_dB,BER2_SR_er,'--g+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER2_SR_th_er,'-g+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER3_SR_er,'--y+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER3_SR_th_er,'-y+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER4_SR_er,'--b+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER4_SR_th_er,'-b+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER5_SR_er,'--k+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER5_SR_th_er,'-k+','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER_quant,'-y+','LineWidth',5);
lgd=legend("PNC and optimal phases Gauss","PNC and optimal phases Gauss th","PNC and quantized optimal phases","PNC and quantized optimal phases th","PNC and optimal phases Prop","PNC and optimal phases Prop th","PNC and random phases","PNC and random phases th","PNC and optimal phases IPM","PNC and optimal phases IPM th","PNC and NO IRS","PNC and NO IRS th","NNC and NO IRS","NNC and NO IRS th");
lgd.FontSize = 25;
%%
figure(2)
semilogy(SNR_dB,BER_er,'--r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("BER",'FontSize', 25);
title("BER at user 3 P(x2hat!=x2) vs SNR",'FontSize', 25);
set(gca,'FontSize',25)
ylim([10e-4 1])
xlim([-10 20])
hold on;
semilogy(SNR_dB,BER_th_er,'-r*','LineWidth',5);
hold on;
semilogy(SNR_dB,BER_q_er,'--y*','LineWidth',5);
hold on;
semilogy(SNR_dB,BER_q_th_er,'-y*','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER2_er,':bs','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER1_th_er,'-bs','LineWidth',5);
% hold on;
semilogy(SNR_dB,BER2_er,'--g+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER2_th_er,'-g+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER4_er,'--c+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER4_th_er,'-c+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER5_er,'--k+','LineWidth',5);
hold on;
semilogy(SNR_dB,BER5_th_er,'-k','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER_quant,'-y+','LineWidth',5);
% lgd=legend("PNC and optimal phases","PNC and quantized phases","PNC and random phases","PNC and No IRS","NNC and NO IRS");
lgd=legend("PNC and optimal phases","PNC and optimal phases th","PNC and quantized optimal phases","PNC and quantized optimal phases th","PNC and random phases","PNC and random phases th","PNC and NO IRS","PNC and NO IRS th","NNC and NO IRS","NNC and NO IRS th");
lgd.FontSize = 25;
%%

figure(3)
semilogy(SNR_dB,Final_MSE_er,'-.r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("MSE",'FontSize', 25);
title("MSE at relay vs SNR",'FontSize', 25);
hold on;
semilogy(SNR_dB,Final_MSE_q_er,':ys','LineWidth',5);
set(gca,'FontSize',25)
% hold on;
% semilogy(SNR_dB,Final_MSE1_er,':bs','LineWidth',5);
hold on;
semilogy(SNR_dB,Final_MSE2_er,'--g+','LineWidth',5);
hold on;
semilogy(SNR_dB,Final_MSE4_er,'--c+','LineWidth',5);
hold on;
semilogy(SNR_dB,Final_MSE5_er,'--k+','LineWidth',5);
lgd=legend("PNC and optimal phases","PNC and quantized optimal phases","PNC and random phases","PNC and No IRS","NNC and NO IRS");
lgd.FontSize = 25;
%%
figure(4)
semilogy(SNR_dB,BER_s1d1_er,'-.r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("BER",'FontSize', 25);
title("BER at D1 P(x1hat!=x1) vs SNR",'FontSize', 25);
set(gca,'FontSize',25)
hold on;
semilogy(SNR_dB,BER_s1d1_th_er,'--mo','LineWidth',5);
lgd=legend("Monte-Carlo","Theoretical");
lgd.FontSize = 25;
% 
%%
figure(5)
semilogy(SNR_dB,BER_rd1_er,'-.r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("BER",'FontSize', 25);
title("BER at D1 P(x1xorx2hat!=x1xorx2) vs SNR",'FontSize', 25);
set(gca,'FontSize',25)
hold on;
semilogy(SNR_dB,BER_rd1_th_er,'--mo','LineWidth',5);
lgd=legend("Monte-Carlo","Theoretical");
lgd.FontSize = 25;
%%
figure(6)
semilogy(SNR_dB,BER_SR_er,'-.r*','LineWidth',5);
xlabel("SNR (dB)",'FontSize', 25);
ylabel("BER",'FontSize', 25);
title("BER at relay P(x1xorx2hat!=x1xorx2) vs SNR",'FontSize', 25);
set(gca,'FontSize',25)
hold on;
semilogy(SNR_dB,BER_SR_th_er,'--mo','LineWidth',5);
lgd=legend("Monte-Carlo","Theoretical");
lgd.FontSize = 25;
ylim([10e-8 1])
% %%
% figure(5)
% semilogy(SNR_dB,BER5_SR_er,'-.r*','LineWidth',5);
% xlabel("SNR (dB)",'FontSize', 25);
% ylabel("BER",'FontSize', 25);
% title("BER at relay P(x1xorx2hat!=x1xorx2) vs SNR",'FontSize', 25);
% set(gca,'FontSize',25)
% hold on;
% semilogy(SNR_dB,BER5_SR_th_er,'--mo','LineWidth',5);
% lgd=legend("Monte-Carlo","Theoretical");
% lgd.FontSize = 25;
% %%
% Gain=BER3_er_wo./BER4_er_wo; 
% 
% Gain1=BER3_er_wo./BER3_er;
% 
% Gain2=BER3_er_wo./BER4_er;
% 
% 
% figure(3)
% semilogy(SNR_dB,Gain,'-.r*','LineWidth',5);
% xlabel("SNR (dB)",'FontSize', 25);
% ylabel("Gain",'FontSize', 25);
% title("Gain in terms of BER vs SNR at user 3",'FontSize', 25);
% set(gca,'FontSize',25)
% hold on;
% semilogy(SNR_dB,Gain1,'--mo','LineWidth',5);
% hold on;
% semilogy(SNR_dB,Gain2,'-c+','LineWidth',5);
% lgd=legend("when only IRS1 is added","when only IRS2 is added","when both IRS1 and IRS2 are added");
% lgd.FontSize = 25;
% grid on;
%%
figure(4)
semilogy([2:iteration],MSE_iter_total(2:end),'-.r*','LineWidth',5);
xlabel("Iteration",'FontSize', 36);
ylabel("MSE",'FontSize', 36);
%title("BER vs SNR at user 3",'FontSize', 25);
set(gca,'FontSize',36)
hold on;
semilogy([2:iteration3],MSE_iter_total3(2:end),'--mo','LineWidth',5);
lgd=legend("SDR with Gaussian Randomization","Proposed Algorithm");
lgd.FontSize = 36;
grid on;

%%
% cd '/Users/kafizoa/Dropbox/Amanat/Amanat_Course_Project/Paper_figures/MSE_vs_SNR_relay_MMSE'
% save('BER_vs_SNR_Butterfly_methods_workspace_new','SNR_dB','iteration','iteration3','MSE_iter_total','MSE_iter_total3')
%%
% close all;
% clear all;
% clc;
% cd '/Users/kafizoa/Dropbox/Amanat/Amanat_Course_Project/Paper_figures/MSE_vs_SNR_relay_MMSE'
% A=load('MSE_vs_SNR_relay_workspace');
%%
% Data_PNC_optimal_th=[A.SNR_dB(1:end); A.BER_th_er(1,1:end)]';
% Data_PNC_optimal_monte=[A.SNR_dB(1:end); A.BER_er(1,1:end)]';
% Data_PNC_random=[A.SNR_dB(1:end); A.BER1_er(1,1:end)]';
% Data_PNC_random_th=[A.SNR_dB(1:end); A.BER1_th_er]';
% Data_PNC_quant=[A.SNR_dB(1:end); A.BER_q_er(1,1:end)]';
% Data_PNC_quant_th=[A.SNR_dB(1:end); A.BER_q_th_er]';
% Data_PNC_no_IRS=[A.SNR_dB(1:end); A.BER4_er(1,1:end)]';
% Data_PNC_no_IRS_th=[A.SNR_dB(1:end); A.BER4_th_er]';
% Data_NNC_no_IRS_monte=[A.SNR_dB(1:end); A.BER5_er(1,1:end)]';
% Data_NNC_no_IRS_th=[A.SNR_dB(1:end); A.BER5_th_er(1,1:end)]';
% % 
% Data_PNC_optimal_th_relay=[A.SNR_dB(1:end); A.BER_SR_er(1,1:end)]';
% Data_PNC_optimal_monte_relay=[A.SNR_dB(1:end); A.BER_SR_th_er(1,1:end)]';
% Data_PNC_random_relay=[A.SNR_dB(1:end); A.BER1_SR_er(1,1:end)]';
% Data_PNC_random_th_relay=[A.SNR_dB(1:end); A.BER1_SR_th_er]';
% Data_PNC_quant_relay=[A.SNR_dB(1:end); A.BER_q_SR_er(1,1:end)]';
% Data_PNC_quant_th_relay=[A.SNR_dB(1:end); A.BER_q_SR_th_er]';


% Data_PNC_optimal=[SNR_dB(1:end); A.Final_MSE_er(1,1:end)]';
% Data_PNC_random=[SNR_dB(1:end); A.Final_MSE2_er(1,1:end)]';
% Data_PNC_quant=[SNR_dB(1:end); A.Final_MSE_q_er(1,1:end)]';
% Data_PNC_no_IRS=[SNR_dB(1:end); A.Final_MSE4_er(1,1:end)]';
% Data_NNC_no_IRS=[SNR_dB(1:end); A.Final_MSE5_er(1,1:end)]';
%%
% figure(1)
% semilogy(A.SNR_dB,A.BER2_er,'-r','LineWidth',5);
% xlabel("SNR (dB)",'FontSize', 25);
% ylabel("BER",'FontSize', 25);
% title("BER at user 3 P(x2hat!=x2) vs SNR",'FontSize', 25);
% set(gca,'FontSize',25)
% hold on;
% semilogy(A.SNR_dB,BER2_er_th,'--r','LineWidth',5);
% hold on;
% semilogy(A.SNR_dB,A.BER_quant,'-b','LineWidth',5);
% hold on;
% semilogy(A.SNR_dB,BER_quant_th,'--b','LineWidth',5);
% hold on;
% semilogy(A.SNR_dB,A.BER4_er,'-c','LineWidth',5);
% hold on;
% semilogy(A.SNR_dB,BER4_er_th,'--c','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER5_th_er,'-mo','LineWidth',5);
% hold on;
% semilogy(SNR_dB,BER_quant,'-y+','LineWidth',5);
% lgd=legend("PNC and optimal phases","PNC and optimal phases: iterative","PNC and zero phases","PNC and random phases","PNC and No IRS","NNC and NO IRS","Theor: NNC and NO IRS");
% lgd.FontSize = 25;
%%
%Data_wt_quantized=[A.SNR_dB(1:end); BER_quantized(1,1:end)]';
%%
% save('Data_PNC_optimal_th.dat','Data_PNC_optimal_th','-ascii');
% save('Data_PNC_optimal_monte.dat','Data_PNC_optimal_monte','-ascii');
% save('Data_PNC_random.dat','Data_PNC_random','-ascii');
% save('Data_PNC_random_th.dat','Data_PNC_random_th','-ascii');
% save('Data_PNC_quant.dat','Data_PNC_quant','-ascii');
% save('Data_PNC_quant_th.dat','Data_PNC_quant_th','-ascii');
% save('Data_PNC_no_IRS.dat','Data_PNC_no_IRS','-ascii');
% save('Data_PNC_no_IRS_th.dat','Data_PNC_no_IRS_th','-ascii');
% save('Data_NNC_no_IRS_monte.dat','Data_NNC_no_IRS_monte','-ascii');
% save('Data_NNC_no_IRS_th.dat','Data_NNC_no_IRS_th','-ascii');
% save('Data_PNC_optimal_th_relay.dat','Data_PNC_optimal_th_relay','-ascii');
% save('Data_PNC_optimal_monte_relay.dat','Data_PNC_optimal_monte_relay','-ascii');
% save('Data_PNC_random_relay.dat','Data_PNC_random_relay','-ascii');
% save('Data_PNC_random_th_relay.dat','Data_PNC_random_th_relay','-ascii');
% save('Data_PNC_quant_relay.dat','Data_PNC_quant_relay','-ascii');
% save('Data_PNC_quant_th_relay.dat','Data_PNC_quant_th_relay','-ascii');


%
% save('Data_PNC_optimal.dat','Data_PNC_optimal','-ascii');
% save('Data_PNC_random.dat','Data_PNC_random','-ascii');
% save('Data_PNC_quant.dat','Data_PNC_quant','-ascii');
% save('Data_PNC_no_IRS.dat','Data_PNC_no_IRS','-ascii');
% save('Data_NNC_no_IRS.dat','Data_NNC_no_IRS','-ascii');
