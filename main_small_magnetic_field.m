%% this program contain three parts
% 1 Get the classical motion and classical varibles
% 2 Get the fringes on the detached plane
% 3 Fit the fringes and extrapolate the data outsides the classical
% boundary
%% First part Classical motion and varibles
syms F B A m e t omega z0 E R 
syms x(t) y(t) z(t) 
syms p_x(t) p_y(t) p_z(t) 
velocity = sqrt(2*m*E)        ;
H(t) =p_z(t)^2/(2*m) +1/2*m*omega^2*(z(t)+p_x(t)/(m*omega))^2 +p_y(t)^2/(2*m)+F*e*z(t);%Hamitonian
syms theta phi
% Hamitonian Equations
equl(1) = diff(x(t),t)   == functionalDerivative(H(t), p_x(t));
equl(2) = diff(y(t),t)   == functionalDerivative(H(t), p_y(t));
equl(3) = diff(z(t),t)   == functionalDerivative(H(t),p_z(t));
equl(4) = diff(p_x(t),t) == -functionalDerivative(H(t), x(t));
equl(5) = diff(p_y(t),t) == -functionalDerivative(H(t), y(t));
equl(6) = diff(p_z(t),t) == -functionalDerivative(H(t),z(t));
vx(t) = diff(x(t), t);
vy(t) = diff(y(t), t);
vz(t) = diff(z(t), t);
% Initial Conditions
cond(1) = x(0) ==  R*sin(theta)*cos(phi);
cond(2) = y(0) ==  R*sin(theta)*sin(phi);
cond(3) = z(0) ==  z0+R*cos(theta);
cond(4) = vx(0) == velocity*sin(theta)*cos(phi);
cond(5) = vy(0) == velocity*sin(theta)*sin(phi);
cond(6) = vz(0) == velocity*cos(theta);
% Solution: and get the expressions of coordinate velocity and momentum
Bs = dsolve(equl(:), cond(:));
xSol(t) = Bs.x;
ySol(t) = Bs.y;
zSol(t) = Bs.z;
vx(t) = diff(xSol,t);
vy(t) = diff(ySol,t);
vz(t) = diff(zSol,t);
p_x = Bs.p_x ;
p_y = Bs.p_y ;
p_z = Bs.p_z ;
% Get the time reaching the detecting plane
tmax = solve(zSol(t) ==0, t);% two solution ,second result above zero can be used
%%
fun_s = p_x*real(vx) + p_y*vy + real(p_z)*real(vz);%-e*zSol*vx*B;
S = int( fun_s ,t, 0,real(tmax(2)));
display(S);
xSol_sin=rewrite(xSol,'sincos');
xSol_real= simplify(xSol_sin,'Criterion','preferReal');
zSol_sin=rewrite(zSol,'sincos');
zSol_real= simplify(zSol_sin,'Criterion','preferReal');
vx_real = simplify(rewrite(vx,'sincos'),'Criterion','preferReal');
p_z_real = simplify(rewrite(p_z,'sincos'),'Criterion','preferReal');
vz_real = simplify(rewrite(vz,'sincos'),'Criterion','preferReal');
Jt= det(jacobian([xSol_real ySol zSol_real],[t theta phi]));
display(Jt);
%%
E = 8.352*10^-5/27.211385     ;  %a.u. 83.52 uev
F = 2.91*10^2/(5.142206*10^11);  %a.u. 291 v/m
B = 1.3*10^-6/(2.35*10^5)    ;  %a.u. 1.1 uT
L_au =  5.2917721092*10^-11;     %(m)
R = 0                        ;  %(a0) as radius of sphere
z0 = 0.5 / L_au               ;  %a.u. 0.5 m
e= 1;
m =1;
omega = B*e/m                 ;  %period
%%
xSol_num = matlabFunction(subs(xSol_real));%(t, phi, theta)
ySol_num = matlabFunction(subs(ySol));%(t, phi, theta)
zSol_num = matlabFunction(subs(zSol));
tmax_sym = subs(tmax(2));
tmax_num = matlabFunction(tmax_sym);%(phi,theta)
xDetect = matlabFunction(subs(xSol(tmax(2))));%(phi,theta)
yDetect = matlabFunction(subs(ySol(tmax(2)))); %this number is in S.I.
action =matlabFunction(subs(S));%(phi, theta)
R =100;
Jacob = matlabFunction(subs(Jt));%%(t phi theta)
%% second part fringes on the detated plane
T_cpu_1 = cputime;
num_phi = 1000;
num_theta   = 200;
num_point = num_phi*num_theta;
initial_angle = zeros(num_point,2);%[theta phi]
Matrix_a = linspace(0,2*pi,num_phi);
Matrix_b = [];
ma = linspace(pi/2000,pi/2*930/1000,num_theta); %%set theta range pi/20000, 
% ma = linspace(pi/2*800/1000,pi/2*970/1000,num_theta);
mb = [];
for j = 1:num_theta
    Matrix_b = [Matrix_b Matrix_a];
    mb = [mb ma(j)*ones(1,num_phi)];
end
initial_angle(:,2)= Matrix_b';
initial_angle(:,1)= mb';
theta_2 = zeros(1,num_point);
phi_2   = zeros(1,num_point);
xmn     = zeros(1,num_point);
ymn     = zeros(1,num_point);
amp1    = zeros(1,num_point); %phi <1/2*pi,trajectory amptitude
amp2    = zeros(1,num_point); %phi >1/2*pi
amp_airya= zeros(1,num_point); %airy approximation
s1      = zeros(1,num_point);
s2      = zeros(1,num_point);
delta_s = zeros(1,num_point);
delta_a = zeros(1,num_point);
flag    = zeros(1,num_point);
sum_s   = zeros(1,num_point);
gamma   = zeros(1,num_point);
delta   = zeros(1,num_point);
zeta    = zeros(1,num_point);
amp_m   = zeros(1,num_point);
amp_p   = zeros(1,num_point);

options = optimoptions('fsolve','Display','none');
% options.FunctionTolerance=10^-10;
% options.OptimalityTolerance = 10^-10 ;
% options.StepTolerance = 10^-10;
for i = linspace(1,num_point,num_point)
    tx = real(tmax_num(initial_angle(i,2),initial_angle(i,1)));
    xmn(i) = xSol_num(tx,initial_angle(i,2),initial_angle(i,1));
    ymn(i) = ySol_num(tx,initial_angle(i,2),initial_angle(i,1));
    amp1(i) = abs(Jacob(0,initial_angle(i,2),initial_angle(i,1))/Jacob(tx,initial_angle(i,2),initial_angle(i,1)))^(1/2);
    s1(i) = action(initial_angle(i,2), initial_angle(i,1));
    fun1 = @(ang)[real(xDetect(ang(1), ang(2))-xmn(i)) ;
        real(yDetect(ang(1), ang(2))-ymn(i))];
    x0 = [initial_angle(i,2), pi-initial_angle(i,1)];
    [mm1,MML,flag(i),MMM] = fsolve(fun1,x0, options);%%mm1 [phi theta]
    phi_2(i) = mm1(1);
    theta_2(i) = mm1(2);
    T2 = real(tmax_num(phi_2(i), theta_2(i)));
    amp2(i)=abs(Jacob(0,phi_2(i),theta_2(i))/Jacob(T2,phi_2(i),theta_2(i)))^(1/2);
    s2(i) = action(phi_2(i), theta_2(i));
    delta_s(i) = s1(i)-s2(i);
    sum_s(i) = s2(i)+s1(i);
    zeta(i) = (3/4*delta_s(i))^(2/3);
% % %     s-wave
%     amp_m(i) = amp2(i) -amp1(i);
%     amp_p(i) = amp2(i) +amp1(i);
% %   p0-wave
%     amp_m(i) = amp2(i)*cos(theta_2(i)) - amp1(i)*cos(initial_angle(i,1));
%     amp_p(i) = amp2(i)*cos(theta_2(i)) + amp1(i)*cos(initial_angle(i,1));
% %     p1-wave
    amp_m(i) = amp2(i)*sin(theta_2(i))*exp(1i*phi_2(i)) -amp1(i)*sin(initial_angle(i,1))*exp(1i*initial_angle(i,2));
    amp_p(i) = amp2(i)*sin(theta_2(i))*exp(1i*phi_2(i)) +amp1(i)*sin(initial_angle(i,1))*exp(1i*initial_angle(i,2));
    gamma(i) = sqrt(pi)*exp(1i*sum_s(i))*amp_p(i)*zeta(i)^(1/4);
    delta(i) = sqrt(pi)*exp(1i*sum_s(i))*amp_m(i)*zeta(i)^(-1/4);
    amp_airya(i) = abs(gamma(i)*airy(-zeta(i))-1i*delta(i)*airy(1,-zeta(i)))^2;
end
%% third part fitting part the fitting area is chosen by the determinating the number of the initial_theta.

num_phi_fit =num_phi ; 
num_theta_fit =200;
number_fit_data = num_theta_fit *num_phi_fit;
xmn_test = zeros(1,number_fit_data);
ymn_test = zeros(1,number_fit_data);
initial_theta = linspace(720/1000*pi/2, 930/1000*pi/2,num_theta_fit);
flag    = zeros(1,number_fit_data);
flag2   = zeros(1,number_fit_data);
% amp_test_1 = zeros(1,number_fit_data);
% amp_test_2 = zeros(1,number_fit_data);
% s_test_1   = zeros(1,number_fit_data);
% s_test_2   = zeros(1,number_fit_data);
sum_s_test = zeros(1,number_fit_data);
zeta_test  = zeros(1,number_fit_data);
amp_test_m = zeros(1,number_fit_data);
amp_test_p = zeros(1,number_fit_data);
% gamma_test = zeros(1,number_fit_data);
% delta_test = zeros(1,number_fit_data);
amp_airya_test = zeros(1,number_fit_data);

options = optimoptions('fsolve','Display','none');
options.Algorithm = 'levenberg-marquardt';
options.StepTolerance = 1.00e-8;
for j = linspace(0,num_phi_fit-1,num_phi_fit)
    phi_fit = j*2*pi/(num_phi_fit-1);
    for i  = linspace(1,num_theta_fit,num_theta_fit)
        tx_test = real(tmax_num(0,initial_theta(i )));
        radius = xSol_num(tx_test,0,initial_theta(i )) - xmn(1);
        xmn_test(i + j*num_theta_fit) = xmn(1) + radius *cos(phi_fit);
        ymn_test(i + j*num_theta_fit) = radius*sin(phi_fit);
        fun2 = @(ang)[real(xDetect(ang(1), ang(2))-xmn_test(i + j*num_theta_fit)) ;
        real(yDetect(ang(1), ang(2))-ymn_test(i + j*num_theta_fit))];
        x0 = [phi_fit , initial_theta(i )];
        x1 = [ phi_fit, pi - initial_theta(i )];
        [mm1,MML,flag(i + j*num_theta_fit),MMM] = fsolve(fun2,x0, options);
        [mm2,MML2,flag2(i + j*num_theta_fit),MMM2] = fsolve(fun2,x1,options);
        tx_test_1 = real(tmax_num(mm1(1),mm1(2)));
        tx_test_2 = real(tmax_num(mm2(1),mm2(2)));
    %     variables below can be saved for debuging, Later can just be taken as  intermediate variables
        amp_test_1 = abs(Jacob(0,mm1(1),mm1(2))/Jacob(tx_test_1,mm1(1),mm1(2)))^(1/2);
        amp_test_2 = abs(Jacob(0,mm2(1),mm2(2))/Jacob(tx_test_2,mm2(1),mm2(2)))^(1/2);
        s_test_1   = action(mm1(1),mm1(2));
        s_test_2   = action(mm2(1), mm2(2));
        sum_s_test(i + j*num_theta_fit)  = s_test_1 +s_test_2;
        zeta_test(i + j*num_theta_fit) = (3/4*abs(s_test_1 - s_test_2))^(2/3);
    %     vaiables below is useful
    % % for s-wave firstly
    %     amp_test_m(i + j*num_theta_fit) = amp_test_2(i + j*num_theta_fit) - amp_test_1(i + j*num_theta_fit);
    %     amp_test_p(i + j*num_theta_fit) = amp_test_2(i + j*num_theta_fit) + amp_test_1(i + j*num_theta_fit);
    % %     for p0-wave
        amp_test_m(i + j*num_theta_fit) = amp_test_2*cos(mm2(2))-amp_test_1*cos(mm1(2));
        amp_test_p(i + j*num_theta_fit) = amp_test_2*cos(mm2(2))+amp_test_1*cos(mm1(2));

        gamma_test = sqrt(pi)*exp(1i*sum_s_test(i + j*num_theta_fit))*amp_test_p(i + j*num_theta_fit)*zeta_test(i + j*num_theta_fit)^(1/4);
        delta_test = sqrt(pi)*exp(1i*sum_s_test(i + j*num_theta_fit))*amp_test_m(i + j*num_theta_fit)*zeta_test(i + j*num_theta_fit)^(-1/4);
        amp_airya_test(i + j*num_theta_fit) = abs(gamma_test*airy(-zeta_test(i + j*num_theta_fit))-1i*delta_test*airy(1,-zeta_test(i + j*num_theta_fit)))^2;
    end
end
warning('off','all')
amp_fit_boundary =[];
x_fit_boundary =[];
y_fit_boundary =[];
for i = linspace(0,num_phi_fit-1,num_phi_fit)
    fit_test_p = amp_test_p(1+i*num_theta_fit:(i+1)*num_theta_fit).*zeta_test(1+i*num_theta_fit:(i+1)*num_theta_fit).^(1/4);  % 
    fit_test_m = amp_test_m(1+i*num_theta_fit:(i+1)*num_theta_fit).*zeta_test(1+i*num_theta_fit:(i+1)*num_theta_fit).^(-1/4);
    radius_test = ((xmn_test(1+i*num_theta_fit:(i+1)*num_theta_fit)-xmn(1)*ones(1,num_theta_fit)).^2*L_au^2+ymn_test(1+i*num_theta_fit:(i+1)*num_theta_fit).^2*L_au^2).^(1/2);
    [fit_test_p1, radius_test1] = prepareCurveData(fit_test_p,radius_test);
    [fit_test_m1 ,radius_test2] = prepareCurveData(fit_test_m,radius_test);
    [sum_s_test1 , radius_test3] = prepareCurveData(sum_s_test(1+i*num_theta_fit:(i+1)*num_theta_fit),radius_test);
    [zeta_test1 ,radius_test4] = prepareCurveData(zeta_test(1+i*num_theta_fit:(i+1)*num_theta_fit),radius_test);
    
    
    [f_test_p ,gof_test_p] = fit(radius_test1,fit_test_p1,'poly2') ;
%     figure();
%     plot(f_test_p,radius_test,fit_test_p)
%     hold on;
%     plot(radius_exp,f_test_p(radius_exp))
%     hold off;
    [f_test_m ,gof_test_m] = fit(radius_test2,fit_test_m1,'poly2');
%     figure();
%     hold on;
%     plot(f_test_m,radius_test,fit_test_m)
%     plot(radius_exp,f_test_m(radius_exp));
%     hold off;
    [fs_test_p , gof_test_sp] = fit(radius_test3,sum_s_test1,'poly1');
    [fs_test_m , gof_test_sm] = fit(radius_test4,zeta_test1,'poly1');
%     radius_exp = linspace(radius_test(1),3*(radius_test(end) - radius_test(1))+radius_test(end),300);
    radius_exp = linspace(radius_test(1),0.9*10^-3,500);
    x_fit = radius_exp*cos(i*2*pi/(num_phi_fit-1))+xmn(1)*L_au;
    y_fit = radius_exp*sin(i*2*pi/(num_phi_fit-1));
    x_fit_boundary = [x_fit_boundary x_fit];
    y_fit_boundary = [y_fit_boundary y_fit];
    amp_test_fit = abs(sqrt(pi)*exp(1i*fs_test_p(radius_exp)).*f_test_p(radius_exp) ...
    .*airy(-fs_test_m(radius_exp))-1i*sqrt(pi)*exp(1i*fs_test_p(radius_exp)).*f_test_m(radius_exp) ...
    .*airy(1,-fs_test_m(radius_exp))).^2;
    amp_fit_boundary = [amp_fit_boundary amp_test_fit'];
end
% x_axis = ymn*L_au;
% y_axis = ymn*L_au;
T_cpu_2 = cputime;
% display(T_cpu_2 - T_cpu_1);
%% Last part save the data and plot out the results
% save fit2d_p1_B11.mat
x_plot = [xmn*L_au x_fit_boundary];
y_plot = [ymn*L_au y_fit_boundary];
amp_plot = [amp_airya amp_fit_boundary];
figure();
scatter(x_plot,y_plot,3,amp_plot,'filled')
mapa = gray;
mapb = flipud(mapa);
colormap(mapb);