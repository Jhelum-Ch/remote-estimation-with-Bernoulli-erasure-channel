%This code computes the constrained performance of the remote-state estimation
%system with one transmitter and one remote estimator connected by a
%Bernoulli packet drop channel
%Written by Jhelum Chakravorty
%Please see for theoretical reference the paper Jhelum Chakravorty and Aditya Mahajan, "Remote-state 
%estimation with packet drop", IFAC workshop, NecSys 2016, Tokyo, Japan

clear all
clc
format long

epsilon_vec = [0 0.3 0.7]; %probability of packet drop

S = 30; %Truncation size. 
p = 0.3;
%betaVec=0.9; %it can be a vector such as betaVec=[0.9;1.0];
betaVec = [0.9;1.0];
err=10^(-4);

alpha = [linspace(0.01,0.1, 9) linspace(0.15, 0.95, 17)]; 

K=1;



l_K_vec = [];
m_K_vec = [];

k_opt = ones(length(epsilon_vec),length(alpha)); % optimal thresholds corresponding to given alpha
color = ['b', 'm', 'k'];


for b=1:length(betaVec)
beta=betaVec(b);
for ep=1:length(epsilon_vec)
    epsilon = epsilon_vec(ep);
    for ind = 1: length(alpha)
        for k=1:S %thresholds
            for K=k:S
                %transition probability matrix of size 2*K+1 for birth-death markov chain
                P=zeros(2*K+1);
                P(1,:)=[1-2*p p zeros(1,2*K-1)];
                P(2*K+1,:)=[zeros(1,2*K-1) p 1-2*p];

                for J=2:2*K
                    P(J,J-1)=p;
                    P(J,J)=1-2*p;
                    P(J,J+1)=p;
                end
     
                h=zeros(2*K+1,1);
                
                h(-k+K+2:k+K) = ones(2*k-1,1);
                h(1:-k+K+1)= epsilon*ones(K-k+1,1);
                h(k+K+1:2*K+1) = h(1:-k+K+1);
   
                P_had=hadamard_prod(h,P);
  
                d= zeros(2*K+1,1);
                ell= zeros(2*K+1,1);
                
                for j=1:2*K+1
                    if j<k+K+1 && j>-k+K+1
                        d(j) = abs(j-K-1); % per-step distortion func
                    else d(j) = epsilon*abs(j-K-1);
                    end
                end
                
              for j=1:2*K+1
                  if j<k+K+1 && j>-k+K+1
                      ell(j) = 0; % per-step distortion func
                  else ell(j) = 1;
                  end
              end

              L_vec = (eye(2*K+1) - beta*P_had)^(-1)*d;
              M_vec = (eye(2*K+1) - beta*P_had)^(-1)*h;
              K_ell_vec = (eye(2*K+1) - beta*P_had)^(-1)*ell;

              l_add(k,K) = L_vec(K+1);
              m_add(k,K) = M_vec(K+1);
	          K_ell_add(k) = K_ell_vec(K+1); 


            end
              trun_L(ep,k) = l_add(k,K);
              trun_M(ep,k) = m_add(k,K);
			  trun_K_ell(ep,k) = K_ell_add(k) ; 	
              trun_D(ep,k) = l_add(k,K)/m_add(k,K);
              trun_N(ep,k) = trun_K_ell(ep,k)/m_add(k,K);
              
         end
        
       
        for l=1:size(trun_M,2)
            if l > 1 
                
                if trun_N(ep,l) <= alpha(ind) && trun_N(ep,l-1)> alpha(ind)
                
                k_opt(ep,ind)= l-1;  %optimal threshold
               
                %Optimal convex combination
                theta_star(ind) = (alpha(ind)-trun_N(ep,k_opt(ep,ind)+1))/(trun_N(ep,k_opt(ep,ind))-trun_N(ep,k_opt(ep,ind)+1));
                %Optimal distortion
                D_opt(ep,ind) = theta_star(ind)*trun_D(ep,k_opt(ep,ind))+ (1-theta_star(ind))*trun_D(ep,k_opt(ep,ind)+1);  %dist_trans func
                break;
                else D_opt(ep,ind) = trun_D(ep,1);
                end
            end
        end
          
           
    end
 

end


%Plots of optimal distortion, D_opt, versus alpha for different packet drop
%prob epsilon
figure
for i=1:length(epsilon_vec)
plot(alpha(1:size(D_opt,2)),D_opt(i,:),color(i));
hold on
end
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$\alpha$','Interpreter','latex');
ylabel(' $D^*_\beta(\alpha)$','Interpreter','latex');
title(['For $\beta$ =', num2str(beta)], 'Interpreter','latex');
legend('\epsilon = 0','\epsilon = 0.3','\epsilon = 0.7','Location','northeast','Interpreter','latex');

%Optional
%Write the results into CSV files
csvwrite(strcat('DT_beta_', num2str(beta), '_trunSize_', num2str(S), '.csv'), D_opt);
csvwrite(strcat('M_beta_', num2str(beta), '_trunSize_', num2str(S), '.csv'), trun_M);
csvwrite(strcat('K_ell_beta_', num2str(beta), '_trunSize_', num2str(S), '.csv'), trun_K_ell);
csvwrite(strcat('N_beta_', num2str(beta), '_trunSize_', num2str(S), '.csv'), trun_N);
end


