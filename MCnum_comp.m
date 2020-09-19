%% Advanced communication assignment
% Steward Qie 22669183 Group 6
%% Initial
clear all
close all
clc
iteration = 10;
Pcrossover = 0.005:0.01:1; 
Monte_Carlo = [1 5 10 50 100]; 
%% Simulation: Monto Carlo& Sum-product method
for i = 1:length(Monte_Carlo)
    for mc = 1:Monte_Carlo(i)
        cer(mc,:) = sumproduct(iteration,Pcrossover); 
    end
    avgcer{i} = mean(cer,1);
end
%% Plot
figure()
for i = 1:length(Monte_Carlo)
    semilogy(Pcrossover,avgcer{i},'Linewidth',1.5);
    xlabel('prob.Crossover');
    ylabel('CER');
    hold on
end
grid on;
legend('Monte_Carlo 1','Monte_Carlo 5','Monte_Carlo 10','Monte_Carlo 50','Monte_Carlo 100')