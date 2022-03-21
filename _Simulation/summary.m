% summary of results
 
f0=0.02; s=0.1; 

% R_NL=N L 1-Closs for s=0.1, M=3, f0=0.02
R_N=[1e2 200 0.12; 1e3 200 0.34; 1e4 200 0.8; 1e5 200 1];
R_L=[ 1e3 4e3 0.065; 1e3 1e3 0.15; 1e3 2e2 0.35; 1e3 50 0.9; 1e3 10 1];
    
% R_M= M 1-Closs for N=1000, L=1000, s=0.1, f0=0.02
R_M=[1 0.1; 3 0.15; 15 0.23; 60 0.26; 300 0.26];

% R_s= s 1-Closs for N=1000, L=1000, M=3, f0=0.02
R_s=[0.4 0.11; 0.2 0.11; 0.1 0.15; 0.05 0.15; 0.025 0.2];

 
%% 
figure(7)
loglog(R_M(:,1),R_M(:,2),'r',R_s(:,1),R_s(:,2),'b')
ylabel('Survived fraction, 1-C_{loss}')
xlabel('Number of crossovers, M, and selection coefficient, s')
title(sprintf('L=1000, s=%g, M=3, N=1000, M=3, f_0=%g',s,f0))

% fitting parameters

a=0.5;  
%%
  figure(1) 
subplot(2,2,1)
N=R_N(:,1); 
surv=R_N(:,3);
x=[min(N) max(N)];
semilogx(N, surv ,'bo', x,0.145*log(x)-0.6,'k')
ylabel('Survived fraction')
xlabel('N')
title(sprintf('L=200, s=%g, M=3, f_0=%g',s,f0))

subplot(2,2,2)
surv  = R_L(:,3);
xx=R_L(:,2).^(-1);
x=[min(xx), max(xx)];
loglog(xx , surv,'bo',x, 5.2*x.^a,'k')
title(sprintf('a=%g',a))
ylabel('Survived fraction')
xlabel('1/L ')
title(sprintf('N=1000, s=%g, M=3, f_0=%g, a=%g',s,f0,a))

subplot(2,2,3)
x =log(R_L(:,1)*f0)./R_L(:,2).^a;
y=R_L(:,3);
xx=[0 0.52];
x1=log(R_N(:,1)*f0)./R_N(:,2).^a;
y1=R_N(:,3);
plot(x, y,'bo',xx,1.9*xx,'k',x1,y1,'r+')
xlabel('log(N f_0)/L^a')
ylabel('Survived fraction')


%% effect of f0 on 1-Closs ar N*f0 >> 1
s=0.1; N=1000; M=3; L=200;
ff=[0.005 0.01 0.02 0.05 0.1];
x =log(N*ff)./L.^a;
y=[0.28 0.25 0.25;0.33 0.35 0.31;0.34 0.40 0.41; 0.51 0.52 0.51; 0.54 0.63 0.66]';
hold on; plot(x, mean(y),'mx'); hold off
title(sprintf('Variation N + L o f0 x \n s=%g, M=3, L=%g, N=%g, f_0=%g, a=%g',s,L, N,f0,a))

%% time point where C + Closs = 1
subplot(2,2,4)
N=1000;
L=[4000 1000 200];

t1loc=(1/s)*log(1/f0);
tast1=[10 18 36]; 
 x1=log(N*f0)./L.^a;
%x1=log(N*s)./L.^a;
%
L=200;
N=[100 1000];
tast2=[10 36];
 x2=log(N*f0)./L^a;
%x2=log(N*s)/L^a;
xx=[0 0.25];
plot(x1, tast1/t1loc,'bo',xx,180*xx/t1loc,'k',x2,tast2/t1loc,'r+')
xlabel('log(N f_0)/L^a')
ylabel('Time where C+C_{loss}=1, t_{1loc}')
title(sprintf('s=%g, M=3, f_0=%g, a=%g',s,f0,a))

%% fastest sites
figure(2)

N=1000;
L=[4000 1000 200 50];

t1loc=(1/s)*log(1/f0);
t1=[5 4 5 5 3.5 2.7 6 3.7 2.9 3.5; 8 7 11 5 8 8 9 8.5 6 7.5; 17.5 14 10 10 12 13.5 17 12.5 12 15; 18 24 33 18 25 17.5 20 21 15 21]/t1loc; 
t1av=mean(t1');
t1std=std(t1');
t1up=t1av+t1std; t1low=t1av-t1std; 
%x1=log(N*f0)*L.^(-a);
x1=log(N*s)*L.^(-a);
%

N=10000;
L=[1000 200 50];

t3=[5.8 4.4 8.2 9.7 10 ; 20 18 15 12.5 17.5; 31 28 26 29 28.5]/t1loc; 
t3av=mean(t3');
t3std=std(t3');
t3up=t3av+t3std; t3low=t3av-t3std; 
%x3=log(N*f0)*L.^(-a);
x3=log(N*s)*L.^(-a);
%
L=200;
N=100;
t2=[9.5 7 3.5 7.5 12 7 8 7.5 12.5 4.5]/t1loc;
t2av=mean(t2);
t2std=std(t2);
t2up=t2av+t2std; t2low=t2av-t2std; 
% x2=log(N*f0)/L^a;
x2=log(N*s)/L^a;

xx=0:0.01: 1;
plot(x3, t3av,'m*',x3, t3up,'m*',x3, t3low,'m*',...
    x1, t1av,'bo',x1, t1up,'bo',x1, t1low,'bo',x2,t2av,'r+',...
    x2,t2low,'r+',x2,t2up,'r+',xx,0.75* xx.^0.75,'k')
xlabel('log(N s)/L^a')
ylabel('Time when max(fsite)=0.5,  t_{1loc}')
title(sprintf(' Fastest sites and 67 percent error interval   s=%g, M=3, f_0=%g, a=%g',s,f0,a))

%%  Probability of fixation at N*f0 << and >> 1

s=0.1; N=1000; M=3; L=200;
ff=[ 0.0001 0.0005 0.001 0.002 0.007 0.02 0.07];
y=[0.005 0.01 0.015; 0.056 0.056 0.062; 0.097 0.108 0.11;0.19 0.16 0.16; 0.285 0.305 0.275;0.335 0.395 0.405; 0.57 0.525 0.6];
y= y' ./(ones(3,1)*ff)/N/s;  % p_fix in terms of s 
% %L=2000
% ff1=[0.001 0.001 0.001]; 
% y1=[0.037 0.037 0.041]./ff1/N/s; 
% %L=4000
% ff2=[0.001 0.001]; 
% y2=[0.025 0.024]./ff2/N/s; 
% %L=500
% ff3=[0.001 0.001 0.001]; 
% y3=[0.074 0.08 0.072]./ff3/N/s;
% %


s=0.1; N=10000; M=3; L=200;

ff1=[1e-4 0.0005 0.001 0.002 0.02];
y1=[0.165 0.1 0.11; 0.385 0.393 0.4; 0.545 0.535 0.54; 0.655 0.655 0.655;0.8 0.825 0.81];
y1= y1' ./(ones(3,1)*ff1)/N/s;

s=0.1; N=1000; M=3; L=1000;

ff2=[1e-4 0.0005 0.001 0.002 0.007 0.02];
y2=[0.008 0.011 0.012; 0.032 0.039 0.037; 0.049 0.049 0.057; 0.069 0.07 0.073; 0.112 0.097 0.104; 0.154 0.151 0.158];
y2= y2' ./(ones(3,1)*ff2)/N/s;


%% Transition to dilute  limit
figure(5)

loglog(ff,[mean(y);mean(y)+std(y);mean(y)-std(y)],'b',...
    ff2,[mean(y2);mean(y2)+std(y2);mean(y2)-std(y2)],'k',...
            ff1,[mean(y1);mean(y1)+std(y1);mean(y1)-std(y1)],'m',[1e-4 0.1],[1 1],'g')
title(sprintf('s=%g, M=%g',s ,M))
xlabel('Initial allele frequency per locus, f_0')
ylabel(' Fixation probability (1-C_{loss})/(f_0 N s)')
%axis([1e-4 2e-2 1e-1 2])
box off

%% meant and stdt for normalized coalescent N* d log Nlineage/dt at N=1000, s=0.1, M=3, T =150
figure(10)
subplot(2,1,1)
L=[2000 500 200 50 10]; N=1e3;
meant=[21 33 25; 35 28 47; 43 36 57; 32 25 21; 22 28 29]';
stdt=[18 24 17; 22 23 27; 18 23 32; 26 21 20; 12 9 26]'; 
semilogx(L,meant,'--r',L,stdt,'--b'); hold on
semilogx(L,mean(meant),'r',L,mean(stdt),'b'); hold on
xlabel('L');ylabel('mean (r) and std (b) norm coal time')
axis([min(L) max(L) 0 100])
title (sprintf('N* d log Nlineage/dt at N=%g, s=0.1, M=3, T =150: mean and std distribution',N))

subplot(2,1,2)
N=[1e3 1e4];  L=500;
meant=[35 28 47; 17 41 51]';
stdt=[22 23 27;  14 19 27]';
semilogx(N,meant,'--r',N,stdt,'--b'); hold on
semilogx(N,mean(meant),'r',N,mean(stdt),'b'); hold on
xlabel('N');ylabel('mean (r) and std (b) norm coal time')
axis([min(N) max(N) 0 100])
title (sprintf('N* d log Nlineage/dt at L=%g, s=0.1, M=3, T =150: mean and std distribution',L))

%% stationary speed with mutation
figure(3)
%  N=1000; L=200; M=3; s=0.1; % f0=0;
% % values of  N*Ub*s^2 and v=<dfw/dt>
% x=[0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 20 200 1000];
% v=[0.0058 0.0072;0.0126 0.0115;0.018 0.0184;0.032 0.0335; 0.044 0.0429; 0.0696 0.0691; 0.084 0.103; 0.129 0.130;0.183 0.202; 0.316 0.332; 1.26 1.22;2.9 2.9];
% v=v';
% %
% N=1000; L=10; M=3; s=0.1;
% x=[0.02 0.05 0.2 1 5 20];
% v=[0.0094 0.0123; 0.018 0.014;   0.021 0.020;    0.0273 0.029; 0.049 0.046;0.1026 0.09755];

N=1000; L=200; M=3; s=0.2;
x=[0.05 0.2 1 5 20 100 1000];
v=[0.032 0.022; 0.0745 0.080 ; 0.175 0.160 ;0.26 0.246 ; 0.49 0.50;0.84 0.86;3.0 3.14];

v=v';
v1=x;
v2=2*log(x/s^2)/log(1/s)^2;
i2=find(v2 > 0 & v2 < v1);
i1=1:min(i2);
loglog(x, mean(v), 'r', x(i1),v1(i1),'m--',x(i2),v2(i2),'b-.')
% hold on
% loglog(x,0.09*sqrt(x),'k')
% hold off

xlabel('Scaled mutation rate, N*Ub*s^2')
ylabel('Adaptation rate, v')
title(sprintf('N=%g L=%g M=%g s=%g f0=0',N,L,M,s))
box off;grid off
