f = 0:100:20000;
a = 0.015;
rho = 1.19875;
c = 344;       

Zr = rad_imped_circular_piston(f,a,rho,c);
ka = 2*pi*f*a;

figure, 
hold on
plot(ka, real(Zr), 'k')
plot(ka, imag(Zr), 'r')

xlabel('ka')
ylabel('Rad imped, real and imag parts')


figure
hold on
plot(ka, abs(Zr))
plot(ka, ka)
plot([ka(1) ka(end)], rho*c*[1 1])


ylim([0 1.1*max(abs(Zr))])
xlabel('ka')
ylabel('Magnitude of rad imped')

legend('|Zr|', 'ka', 'rho*c')
