clear all
close all
clc
iteration = [1 5 10 50 100];
Pcrossover = 0.0001:0.01:1; 
Monte_Carlo = 1; 
%% Simulation: Monto Carlo& Sum-product method
for i = 1:length(iteration)
    iter = iteration(i); % iteration
    for mc = 1:Monte_Carlo
        cer(mc,:) = sumproduct(iter,Pcrossover); 
    end
    avgcer{i} = mean(cer,1);
end
%% Plot
figure()
for i = 1:length(iteration)
    semilogy(Pcrossover,avgcer{i},'Linewidth',1.5);
    xlabel('prob.Crossover');
    ylabel('CER');
    hold on
end
grid on;
legend('iteration times 1','iteration times 5','iteration times 10','iteration times 50','iteration times 100')