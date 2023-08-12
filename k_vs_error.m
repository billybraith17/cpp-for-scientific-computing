data = importdata('output2.dat') ; 
k = data(:,1) ;
L2 = data(:,2) ; % give required column of index 
inf = data(:,3) ;
loglog(k, L2, k, inf);
xlabel('$k$','interpreter','latex', 'FontSize', 24)
ylabel('Error','interpreter','latex', 'FontSize', 24)
legend('$L^2$ error','$L^\infty$ error','interpreter','latex', 'FontSize', 24)