%% cinemática inversa para la generación de un círculo y generación de parámetros para que la
%% trayectoria de la artículación se identifique por un HGO y sea el exosistema
clc
clearvars
close all

FS=14;

 l=0.5;
 m1=0.4;
 m2=0.3;
 g=9.81;

 t1_op=deg2rad(85);  %theta 1 o x1 de operación para el atractor
 t2_op=deg2rad(-109);   %theta 2 o x3 de operación para el atractor

 Ts=0.0005;    
 epsilon=0.02; 
 TF=6;

 t1_ini=1.4941; 
 t2_ini=-1.9134; 

%% cálculo del espacio de trabajo
l1 = .5;
l2 = .5; 

theta1 = deg2rad(15):0.1:deg2rad(120); 
theta2 = deg2rad(-150):0.1:deg2rad(-10); 

[THETA1,THETA2] = meshgrid(theta1,theta2); 

X = l1 * cos(THETA1) + l2 * cos(THETA1 + THETA2); 
Y = l1 * sin(THETA1) + l2 * sin(THETA1 + THETA2); 



%% C I R C U L O
%  n=Tf/Ts;
% tiempo=0:Ts:(Ts*(n-1));
%  xy = circle([0.35 0.75], 0.15, 'n', n);
%  px=xy(1,:);
%  py=xy(2,:);



% R E C T A N G U L O
segmentos=5;
n=(TF/segmentos)/Ts;
tiempo=0:Ts:(Ts*(n*segmentos-1));

T0=[0.9133 0.4073 0 0.4949
   -0.4073 0.9133 0 0.2949
    0 0 1 0
    0 0 0 1];

T1=[0.9742 -0.2258 0 0.6
    0.2258  0.9742 0 0.6
    0 0 1 0
    0 0 0 1];

T2=[0.9239 -0.3827 0 0.75
    0.3827  0.9239 0 0.6
    0 0 1 0
    0 0 0 1];

T3=[0.5876 0.8091 0 0.75
    0.8091 0.5876 0 -0.2
    0 0 1 0
    0 0 0 1];

T4=[0.3552 0.9348 0 0.6
   -0.9348 0.3552 0 -0.2
    0 0 1 0
    0 0 0 1];

tc0 = ctraj(T0, T1, n);
tc1 = ctraj(T1, T2, n);
tc2 = ctraj(T2, T3, n);
tc3 = ctraj(T3, T4, n);
tc4 = ctraj(T4, T1, n);

for i=1:n
    px0(i)=tc0(1,4,i);
    py0(i)=tc0(2,4,i);

    px1(i)=tc1(1,4,i);
    py1(i)=tc1(2,4,i);

    px2(i)=tc2(1,4,i);
    py2(i)=tc2(2,4,i);

    px3(i)=tc3(1,4,i);
    py3(i)=tc3(2,4,i);

    px4(i)=tc4(1,4,i);
    py4(i)=tc4(2,4,i);
end

px=[px0 px1 px2 px3 px4];
py=[py0 py1 py2 py3 py4];

% px=[px1 px2 px3 px4];
% py=[py1 py2 py3 py4];

%% grafica espacio de trabajo y trayectoria
figure(1)
  hold on
  grid on 
  plot(X(:),Y(:),'r.');
  axis equal;
  xlabel('p_x (m)','fontsize',FS)
  ylabel('p_y (m)','fontsize',FS)
  plot(px,py,'b-','LineWidth',1)

  xl1=l*cos(t1_ini);
  yl1=l*sin(t1_ini);
  xl2=l*cos(t1_ini)+l*cos(t1_ini+t2_ini);
  yl2=l*sin(t1_ini)+l*sin(t1_ini+t2_ini);

  line([0,xl1],[0,yl1],'Color','black','LineWidth',1.5)
  line([xl1,xl2],[yl1,yl2],'Color','black','LineWidth',1.5)

  legend('','p_{x\midref}, p_{y\midref}','fontsize',FS)

%% solución para theta_1=alfa-beta

% px=[0.4949 0.6 0.75 0.75 0.6] 
% py=[0.2949 0.6 0.6 -0.2 -0.2]

alpha=atan2(py,px);

cos_beta=(px.^2+py.^2+l1^2-l2^2)./(2*l1*sqrt(px.^2+py.^2));
sin_beta=sqrt(1-cos_beta.^2);

beta=atan2(sin_beta,cos_beta);

theta_1a=alpha-beta;
theta_1b=alpha+beta;

%% solución para theta_2

cos_theta2=(px.^2+py.^2-l1^2-l2^2)./(2*l1*l2);
sin_theta2=sqrt(1-cos_theta2.^2);

theta_2a=atan2(sin_theta2,cos_theta2);
theta_2b=-atan2(sin_theta2,cos_theta2);

%% r e s u l t a d o s
EU=[theta_1b' theta_2b'];
ED=[theta_1a' theta_2a'];


theta1=[tiempo',theta_1b'];
theta2=[tiempo',theta_2b'];
PX=[tiempo',px'];
PY=[tiempo',py'];
trayectoria=[px',py'];

figure(2)
  hold on
  grid on   
   % axis equal;
  plot(tiempo, theta_1b,'r--')
  xlabel('tiempo (s)','fontsize',FS)
  ylabel('\Psi_1 (rad)','fontsize',FS)



figure(3)
  hold on
  grid on 
%   axis equal;
  plot(tiempo, theta_2b, 'r--')
    xlabel('tiempo (s)','fontsize',FS)
  ylabel('\Psi_2 (rad)','fontsize',FS)

% figure(4)
%   hold on
%   grid on   
% %    axis equal;
%   plot(tiempo, theta_1b,'b--')
%   plot(tiempo, theta_2b, 'r--')
%   xlabel('time (s)','fontsize',FS)
%   ylabel('\Psi (rad)','fontsize',FS)
%   legend('\psi_1', '\psi_2','fontsize',FS)


%% Generación de parámetros para el HGO en simulink
% Ts=0.05;
% 
% S_hg=[0 1;0 0];
% S_disc=(eye(2)+Ts*S_hg);
% 
% Q_hg=[1 0];
% 
% r_p=[-.05 -.1];
% pol=poly(r_p);
% a1=pol(2);
% a2=pol(3);
% 
% 
% e=0.01;
% h1=a1/e;
% 
% h2=a2/(e^2);
% H=[h1;h2];
% 
% H_dis=Ts*H;

  %% ecuaciones de estado

 syms t1 t2 dt1 dt2 x1 x2 x3 x4 u1 u2 

M=[(1/3)*m1*(l^2) + (4/3)*(m2)*(l^2) + m2*cos(t2)*l^2                 (1/3)*(m2)*(l^2) + 0.5*m2*(l^2)*cos(t2)
      (1/3)*m2*(l^2) + 0.5*m2*(l^2)*cos(t2)                              (1/3)*m2*(l^2)];
      
C=[-0.5*m2*sin(t2)*(l^2)*(dt2^2) - m2*sin(t2)*(l^2)*dt1*dt2
         0.5*m2*sin(t2)*(l^2)*(dt1^2)];
     
G=[0.5*m1*g*l*cos(t1) + 0.5*m2*g*l*cos(t1+t2) + m2*g*l*cos(t1)
          0.5*m2*g*l*cos(t1+t2)];
      
      alpha=M\([u1;u2] - C - G);
    
%%
      
dx1=dt1; %x3
dx2=alpha(1); %x4
dx3=dt2; %dx3 
dx4=alpha(2); %dx4

dx=[dx1
    dx2
    dx3
    dx4];
y=[x1;x2;0;0];
%% Cálculo de las matrices A B C

A=jacobian(dx,[t1 dt1 t2 dt2]);
B=jacobian(dx,[u1 u2]);
C=jacobian(y,[t1 dt1 t2 dt2]);

%% Linealización de las matrices A B C

t1=(t1_op); 
t2=(t2_op);
dt1=0; 
dt2=0;
u=(solve(dx2,dx4,[u1,u2]));
U1=u.u1;
U2=u.u2;
u1=eval(U1);
u2=eval(U2);
ue=[u1;u2];    % torques en el punto de operación
A=eval(A); % matriz A en tiempo continuo
B=eval(B); % matriz B en tiempo continuo

Ad=eye(4)+Ts*A;
Bd=Ts*B;

C1=[1 0];
C2=[1 0];

Cm=blkdiag(C1,C2);


Kd=place(Ad,Bd,[0.62 0.6 0.62 0.6]);    %ganancias que estabilizan al sistema
%Kd=place(Ad,Bd,[0.8 0.95 0.9 0.95]) %para el lineal funcional mas o menos


% Kd=place(Ad,Bd,[0.9 0.95 0.95 0.9])
% Kd=place(Ad,Bd,[0.97 0.98 0.97 0.98])


%% Contrucción del observador de alta ganancia en este caso son 2, uno por cada articulación
% Matrices A del HGO según yo esto se utliza cuando se conoce psi
% Ah1=[0 1
%      0 0];
% 
% Ah2=[0 1
%      0 0];
% 
% SA = blkdiag(Ah1,Ah2);
% 
% Ah1d=eye(2)+Ts*Ah1;
% Ah2d=eye(2)+Ts*Ah2;
% 
% SAd = blkdiag(Ah1d,Ah2d);

% Matrices H del HGO

pol_1=poly([-1, -1]); % obtener el polinomio de Hurtwitz que tenga las raices -1
pol_2=poly([-1, -1]); % obtener el polinomio de Hurtwitz que tenga las raices -1
% 
alpha0_1=pol_1(1,1);
alpha1_1=pol_1(1,2);
alpha2_1=pol_1(1,3);
alpha0_2=pol_2(1,1);
alpha1_2=pol_2(1,2);
alpha2_2=pol_2(1,3);

H1=[alpha1_1/epsilon alpha2_1/epsilon^2]';
H2=[alpha1_2/epsilon alpha2_2/epsilon^2]';

H=[H1;H2];

S1=[-alpha1_1/epsilon   1
    -alpha2_1/epsilon^2 0];

S2=[-alpha1_2/epsilon   1
    -alpha2_2/epsilon^2 0];

%SH=blkdiag(H1,H2); % si se trabaja en tiempo continuo se usan estas
%Si=blkdiag(S1,S2);

% S_i=[zeros(1) eye(1); zeros(1,2)]-[H1 zeros(2,1)] así las genera el dr. % Meda

%discretización de las matrices Si & Hi del HGO

S1d=eye(2)+Ts*S1;
S2d=eye(2)+Ts*S2;

H1d=Ts*H1;
H2d=Ts*H2;

Sid=blkdiag(S1d,S2d);
SHd=blkdiag(H1d,H2d);

Q1=[1 0];
Q2=[1 0];

Q=blkdiag(Q1,Q2);

%% resolver el problema de regulación
syms P_11 P_12 P_13 P_14 P_21 P_22 P_23 P_24 P_31 P_32 P_33 P_34 P_41 P_42 P_43 P_44 G_11 G_12 G_13 G_14 G_21 G_22 G_23 G_24

PI=[P_11 P_12 P_13 P_14
    P_21 P_22 P_23 P_24
    P_31 P_32 P_33 P_34
    P_41 P_42 P_43 P_44];

GAMMA=[G_11 G_12 G_13 G_14
       G_21 G_22 G_23 G_24];

EQ1= Ad*PI + Bd*GAMMA - PI*Sid - PI*SHd*Q;
EQ2 = Cm*PI-Q;

sol=solve(EQ1,EQ2,[P_11 P_12 P_13 P_14 P_21 P_22 P_23 P_24 P_31 P_32 P_33 P_34 P_41 P_42 P_43 P_44 G_11 G_12 G_13 G_14 G_21 G_22 G_23 G_24]);

Pi=eval([sol.P_11 sol.P_12 sol.P_13 sol.P_14
    sol.P_21 sol.P_22 sol.P_23 sol.P_24
    sol.P_31 sol.P_32 sol.P_33 sol.P_34
    sol.P_41 sol.P_42 sol.P_43 sol.P_44]);

Gamma=eval([sol.G_11 sol.G_12 sol.G_13 sol.G_14
       sol.G_21 sol.G_22 sol.G_23 sol.G_24]);

%% comprobación
EQ1=( Ad*Pi + Bd*Gamma - Pi*Sid - Pi*SHd*Q); 
EQ2= ( Cm*Pi-Q);

rn=height(px')
traj=[px' py' zeros(rn,1)];
traj_adj=traj(1:10:end,:);



