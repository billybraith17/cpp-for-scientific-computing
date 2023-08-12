data = importdata('output.dat') ; 
N = data(:,1) ;
inf = data(:,2) ; % give required column of index 
MSE = data(:,3) ;
% plot(N, MSE);
loglog(N, MSE, N, inf, N, N.^(-2));
% plot(x, exp(-x))
xlabel('$N$','interpreter','latex', 'FontSize', 24)
ylabel('Error','interpreter','latex', 'FontSize', 24)
legend('$L^2$ error', '$L^\infty$ error','$N^{-2}$','interpreter','latex', 'FontSize', 24)