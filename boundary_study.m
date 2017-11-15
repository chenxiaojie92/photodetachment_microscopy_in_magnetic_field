%% Profile of the pattern and approximation
% When the profile of the pattern is far from a circle,The algorithm is not 
% useful when appling to arbitrary direction  and we need to study the profile 
% of the pattern in order to set the initial value for the 'fsolve' and get the 
% initial information of $\theta$ and $\phi$.
%% Classical motion

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
fun_s = p_x*real(vx) + p_y*vy + real(p_z)*real(vz);%-e*zSol*vx*B;
S = int( fun_s ,t, 0,real(tmax(2)));

xSol_sin=rewrite(xSol,'sincos');
xSol_real= simplify(xSol_sin,'Criterion','preferReal');
zSol_sin=rewrite(zSol,'sincos');
zSol_real= simplify(zSol_sin,'Criterion','preferReal');
vx_real = simplify(rewrite(vx,'sincos'),'Criterion','preferReal');
p_z_real = simplify(rewrite(p_z,'sincos'),'Criterion','preferReal');
vz_real = simplify(rewrite(vz,'sincos'),'Criterion','preferReal');
Jt= det(jacobian([xSol_real ySol zSol_real],[t theta phi]));
%%
E = 8.352*10^-5/27.211385     ;  %a.u. 83.52 uev
F = 2.91*10^2/(5.142206*10^11);  %a.u. 291 v/m
B = 60*10^-6/(2.35*10^5)    ;  %a.u. 1.1 uT
L_au =  5.2917721092*10^-11;     %(m)
R = 0                        ;  %(a0) as radius of sphere
z0 = 0.5 / L_au               ;  %a.u. 0.5 m
e= 1;
m =1;
omega = B*e/m                 ;  %period
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
%% Study the shape of the pattern
% After we find that the phi_pattern will not change much when I choose a number 
% of phi in the system.Then I make spline connecting the phi_sys and phi_pattern(they 
% can taken as the same when the pattern is 

x_boun =[];
y_boun =[];
for i = linspace(0,pi/2,101) % num_of_theta 100
    for j = linspace(0,pi,101) % num_of_phi 100
        t_boun = real(tmax_num(j , i));
        x_boun = [x_boun ;xSol_num(t_boun,j,i)*L_au i j];
        y_boun = [y_boun ;ySol_num(t_boun,j,i)*L_au i j];
    end
end
phi_pattern =[];
for i_phi = linspace(0,pi,101)
    num_boun = find(x_boun(:,3) ==i_phi);
    phi_pat = atan(y_boun(num_boun(end),1)/(x_boun(num_boun(end),1)-x_boun(num_boun(1),1)));
    phi_pat = phi_pat *(phi_pat >= 0)+(phi_pat +pi)*(phi_pat < 0);
    phi_pattern =[phi_pattern phi_pat];
end
pattern_sys = spline(phi_pattern,linspace(0,pi,101));
%% test
% plot(ppval(pattern_sys,linspace(25/40*pi,35/40*pi,20)))
%%