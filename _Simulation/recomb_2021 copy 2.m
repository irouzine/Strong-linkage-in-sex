% Updated 15-Jul-2021
 
% for Fig. 2 in Batorsky et al 2011 
% recomb_0('constant',0.05,0.1,1,500,1000,10,1500,1)

function recomb_2021(distribution_s,r,s0,a,M,L,N,tf,f0,run)
%homedir = '/cluster/shared/rbator01/recomb_model_test/new_model/';
% file location and title for saving figs
homedir = '~/Desktop/Recombination/figs';
filename = sprintf('_%s_%g_%g_%g_%g_%g_%g_%g_%g_%g',...
    distribution_s,r,s0,a, M,L,N,tf,f0,run);

% Arguments:
% distribution_s is one of the four forms below
% r recombination rate per genome
% s0 average selection coefficient
% a index in the power distribution_s (not used)
% L number of loci
% N population number
% M crossover number
% tf full time in genertions
% run the number of the run

% random seed, so can average over runs.
rng(run)

%f0=0.01;                        % preset initial frequency of alleles
muesc=0;                        % standing variation f0(s): deleterious mutation rate during escape
mu=0;% 3e-5;                    % recurrent beneficial mutation rate 
tesc=0;                         % time interval for a kind of initial conditions
T= 0:tf ;                       % times
N_sample=30;                    % Sampling segment for w2

%% Distributions of s and parameter titles of fig 1 

switch distribution_s
    case 'power' 
        % 1: a*s0^a/(s+s0)^(a+1) 
        s=s0*(rand(1,L).^(-1/a)-1); 
        ts=sprintf('1/s^{a+1},s>s0. N=%g,r=%g,L=%g,s0=%g,\n a=%g,f0=%g,T=%g,M=%g',N,r,L,s0,a,f0,max(T),M);
    case 'lognormal'
        % 1/s*exp[-log^2(s/s0)/2/a^2]
        s=s0*exp(a*randn(1,L));
        ts=sprintf('1/s*exp[-log^2(s/s0)/2/a^2]\n N=%g,r=%g,L=%g,s0=%g,a=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n T=%g,M=%g',N,r,L,s0,a,f0, muesc,mu,tesc,max(T),M);
    case 'exponential'
        % exp(-s/s0) %% actually density is, g(s) = (1/s0)*exp(-s/s0)
        s=-s0*log(rand(1,L));
        ts=sprintf('exp(-s/s0)\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n T=%g,M=%g',N,r,L,s0,f0, muesc,mu,tesc,max(T),M);
    case 'const'
        s = s0*ones(1,L);
        ts=sprintf('const, N=%g,r=%g,L=%g,s0=%g\n f0=%g, T=%g, M=%g, run=%g',N,r,L,s0,f0, max(T),M, run); % no mutation
end


%% Initial settings
 
fsite=zeros(length(T),L); Plinkav=fsite; %fanc=fsite;

fav=zeros(length(T),1);  Cav=fav; w2av=fav; 
nlineages=fav; Closs=fav; fpol=fav; Cpol=fav; Vark=fav;
w2=zeros(1,N_sample); Plink=zeros(N_sample,L);  C=zeros(1,L); 
nprog=zeros(N,1);


figure(2);clf
 
Knew=zeros(N,L);  % binary sequence matrix
Anew=zeros(N,L);  % ancestor label matrix

tint=tf/5; % plot time interval 
col='rgbmkrgbmkrgbmkrgbmk';

%% Initial population

% Case 1: randomly distributed good alleles with fixed frequency
if f0~=0
    K=(rand(N,L) < f0); 
    % Matrix K: Each row is a genome, sequence of 0 and 1
    
% Case 2: dynamic distribution based on the antigenic escape compensation model: 
% Alleles acumulate during tesc as deleterious, then change the sign of s. 
% Model valid if 1-site deterministic: f0(s) > 1/Ns <=> mu*N > 1 for s > 1/tesc)
elseif muesc~=0
	f0=muesc*s.^(-1).*(1-exp(-s*tesc));       % row of initial frequencies 
	K=(rand(N,L) < ones(N,1)*f0);
    
% Case 3: no alleles at start
else
    K=zeros(N,L);
end

A=(1:N)'*ones(1,L);  % Initial matrix of ancestor labels

%% Evolution starts...
for t=T
    % beneficial mutation
    K = K | sprand(N,L,mu);
    
    % recombination of some pairs with parent replacement
    npairs=round(r*N/2);
    ii=ceil(rand(npairs,2)*N); i1=ii(:,1); i2=ii(:,2);  % 2 columns of random indices
    for i=1:npairs
        % generating random 1 or 0 for each site, with prob. M/L and 1-M/L,
        % respectively; even/odd xx shows site segments copied from 1st parent,
        % 2nd, 1st, 2nd etc
        xx=cumsum(rand(1,L) < M/L); 
        first=(round(xx/2)==xx/2);    % sites copied from 1st parent
        prog1=K(i1(i),:).*first+K(i2(i),:).*(1-first);    % progeny sequence
        progA=A(i1(i),:).*first+A(i2(i),:).*(1-first);    % progeny ancestor labels
        %prog2=K(i2(i),:).*first+K(i1(i),:).*(1-first);    % throw away complementary progeny
        K(i1(i),:)=prog1; %K(i2(i),:)=prog2;               % only one parent is replaced (1/11)
        A(i1(i),:)=progA;                                  %  replacing ancestor labels
    end
    
%% Random sampling of progeny and natural selection with a  broken stick
    % column of N log-fitnesses 
    % state with 0s only has w=0 (fitness 1) by definition
    w=K*(s'); 

    nprogav=exp(w)/mean(exp(w));
    b2=cumsum(nprogav);
    b1=[0;b2(1:end-1)];
    X=rand(N,1)*N;
    for i=1:N
        nprog(i)=sum( X > b1(i) & X < b2(i));
    end
  
%     %% Random drift: divide in 2 or die 
%     % Same as Poisson, Var[n]=E[n], Neff=N for drift/tree; not for few copies of an allele

%     nprogav=min(2,exp(w)/sum(exp(w)));      % average progeny numbers <= 2
%     nprogav=nprogav/mean(nprogav);          % renormalizing
%     nprog=rand(N,1) < nprogav/2;            % survive or die?
%     % balancing survivals and deaths to keep population constant
%     disbal=1;
%     while disbal ~=0
%         isur=find(nprog); idied=find(1-nprog);
%         disbal=length(isur) - N/2;
%         if disbal > 0
%             iflip=isur(ceil(rand(1,disbal)*length(isur)));
%             nprog(iflip)=zeros(1,disbal); % replacing extra 1 by 0
%         elseif disbal < 0
%             iflip=idied(ceil(rand(1,-disbal)*length(idied)));
%             nprog(iflip)=ones(1,-disbal);
%         end
%     end
%     nprog=nprog*2;
%     
    %% updating population
    is=[0;cumsum(nprog(1:(N-1)))];
    for i=1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*K(i,:);
            Anew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*A(i,:);
        end
    end
    K=Knew;
    A=Anew;
    sK=size(Knew);sK=sK(1);
    if sK  ~= N, disp('N not conserved'),return;end
    
    % Evolution step ended
    
    %% Memorizing observables  
 
    ipol=find(~all(1-K,1));         % polymorphic site  indices
    Closs(t+1)=1-length(ipol)/L;    % fraction of sites with extinct alleles
    fsite(t+1,:)=mean(K);           % allele frequencies at all sites
    fav(t+1)=mean(mean(K));         % average allele frequency per any site or per genome
    fpol(t+1)=mean(mean(K(:,ipol)));% average allele frequency per polymorphic site 
    Vark(t+1)=(std(w)/s0)^2;        % variance of the allele number in a genome
    
    % Sampling random pairs for w2
    ii=ceil(rand(N_sample,2)*N); i1=ii(:,1); i2=ii(:,2);  % 2 columns of N_sample random indices between 1 and N

    for i=1:N_sample
        w2(i)=sum(K(i1(i),:)~=K(i2(i),:))/2; % half-distance
        for j=0:(L-1)
            Plink(i,j+1)=mean(A(i1(i),1:(L-j))==A(i1(i),(j+1):L)); % for LD
            %  probability of 2 sites at distance j in one genome to have the same ancestor 
        end
    end
    % Average over the sample 
    w2av(t+1)= mean(w2); 
    nlineages(t+1)=length(unique(reshape(A,[1,N*L])));  % total dictinct lineage number
    Plinkav(t+1,:)= mean(Plink);  % for LD
    
    % Ancestor spectrum and C
    nanc =[];nmax=zeros(1,L);
    for i=1:1:L
        %Au=unique(A(:,i));
        %spec_anc=sum(ones(N,1)*Au'== A(:,i)*ones(1,length(Au)))/N;
        spec_anc=hist(A(:,i),min(A(:,i)):max(A(:,i)));
        C(i)=sum(spec_anc.^2)/N^2;
        % clone sizes
        if ~all(1-K(:,i))
            nanc =[nanc  spec_anc];
            nmax(i)=max(spec_anc);
        end % non-lost sites only
        %fanc(t+1,1:length(spec_anc))=spec_anc;
    end
  
    Cav(t+1)=mean(C(C>0));
    Cpol(t+1)=mean(C(ipol));
    nmax=nmax(ipol);
    
    %% Plotting the wave and ancestral clone  spectrum at some time points
    if tint*round(t/tint)==t
        c=col(round(t/tint)+1);
        figure(2)
     subplot( tf/tint +2,1,1)
        [nn,xx]=hist( w );                % histogram of fitness among genomes
        plot(xx,nn,c)
        hold on
     subplot( tf/tint +2,2, 2*t/tint+3)
        [nn,xx]=hist( nanc(nanc>0),min(30,length(nanc(nanc>0))) );                % histogram of fitness among genomes
        semilogy(xx/N,nn,[c 'o'])
        title(sprintf('t=%g',t))
        axi=axis;axi(1:2)=[0 1 ];axis(axi);
     subplot( tf/tint +2,2, 2*t/tint+4)
        [nn,xx]=hist( nmax(nmax>0),min(30,length(nmax(nmax>0))) );                % histogram of fitness among genomes
        semilogy(xx/N,nn,[c 'o'])
%         histogram(nanc(nanc>0),30,'FaceColor', c )
        title(sprintf('t=%g',t))
        axi=axis;axi(1:2)=[0 1 ];axis(axi);
    end

end

% Evolution has ended

%% Final plots
figure(2)
% The wave
subplot( tf/tint +2,1,1)
hold off
xlabel('Log fitness,  w  ');
ylabel('Density distribution')
title(ts)   % title  with parameter values
 subplot( tf/tint +2,2,2*tf/tint +3)
 xlabel('Clone size, n_{anc} ');
  subplot( tf/tint +2,2,2*tf/tint +4)
 xlabel('Max clone /site, n_{max} ');
 
%%
figure(1)
% averages
% sliding window averaging
 window=round(tf/20);
 w2av_s=zeros(size(T));
for i=1:length(T)
   w2av_s(i)=mean(w2av(max(i-window,1):min(i+window,length(T))));
end
w2av=w2av_s;
w2av=reshape(w2av, size(T));
fav=reshape(fav, size(T));
Closs=reshape(Closs, size(T));
Cav=reshape(Cav, size(T));
fpol=reshape(fpol, size(T));
Cpol=reshape(Cpol, size(T));
Vark=reshape(Vark, size(T)); % Variance k
xx=length(T);
V=L*(fav(2:xx)-fav(1:(xx-1))); % speed
V=[V(1),V];
Fisher= V./Vark/s0; % checking Fisher theorem
 

f1site=f0*(f0+(1-f0)*exp(-s0*T)).^(-1); % 1-site theory
%size(T),size(w2av),size(Cav),size(Closs),size(f1site) 
subplot(2,1,1)

plot(T,fav,'b',T,fpol,'b+',T,f1site,'b--',T,w2av/L,'g',...
    T,Cav,'m',T,Cpol,'m+',T,1-Closs,'r',T,max(fsite'),'k');

title(sprintf('f_{av} b   f_{pol} b+   f_{1site} b- -   w2av/L g   C m   C_{pol} m+ 1-C_{loss} r'))
xlabel('Time, t');
ylabel(sprintf('All.frequency f f_{pol}, half-distance w^2 \n allele loss C_{loss}, id.descent C C_{pol}'))
axi=axis;axi(3:4)=[0 1 ];axis(axi);
grid
text(tf/2,0.75,ts)

subplot(2,1,2)
plot(T, Vark./w2av, 'co-',T, Fisher, 'c+',...
    T, w2av./fpol./(1-fpol)./(1-Closs) /L, 'k+', T, w2av./fav./(1-fav)./(1-Cav)/L, 'ko',T, w2av./fav./(1-fav-Closs)./(1-Cav).*(1-Closs)/L, 'k')
title(sprintf('Var[k]/w^2  co-   V/Vark/s0 c+   w^2/f_{pol}/(1-f_{pol})/(1-C_{loss})/L  k+ \n  w^2/f_{av}/(1-f_{av})/(1-C_{av})/L  ko  w^2_{av}(1-C_{loss})/f_{av}/(1-f_{av}-C_{loss})/(1-C_{av})/L k-'))
xlabel('Time, t');
ylabel(sprintf('Normalized w^2 and Var[k]'))
axi=axis;axi(3:4)=[0 2 ];axis(axi);
grid

figure(4)
% Histogram fsite
i=0;
for t=0:tf/5:tf
    i=i+1;
    %plot(t*ones(1,L),fsite(t+1,:), [c(i+1),'o'])

    ii=find(fsite(t+1,:) > 0 & fsite(t+1,:) < 1);
%     [nn,xx]=hist(fsite(t+1,ii));                 
%     plot(xx,nn,c(i))
subplot(7,1,i+1)
    histogram(fsite(t+1,ii),'FaceColor',col(i));
    %ir=ceil( length(nn)*rand);
     title( sprintf('t=%g',t))
axi=axis; axi(1:2)=[0 1 ]; axis(axi);
end  
xlabel('fsite ');
%ylabel('# of sites')



subplot(7,1,1)
%     f0=mean(fsite ==0, 2);
%     f1=mean(fsite ==1, 2);
plot(T,fsite')
xlabel('Time');
%ylabel('f0 and f1')
ylabel ('fsite')
axi=axis; axi(3:4)=[0 1 ]; axis(axi);
 title( ts)
grid


% figure(3) 
% % LD
% for t=0:10:tf
%     plot(1:L,Plinkav(t+1,:),col(round(t/10)+1))
%     hold on
% end 
% hold off
% xlabel('Positional distance')
% ylabel('Probability of same ancestor, Plinkav')
% title('t=0:10:tf')


% %% Coalescent probability Pid and its analytic connection to w2 from Eq. 96
% figure(4)
% The upper subplot for Pid is OK: measurable  but noisy
% The lower subplot shows that the old analytic theory is very bad

%  Pid = Cav(2:end)-Cav(1:(end-1));
% 
%  Pid=[Pid(1);Pid]./(1-Cav); 
%  Pid=reshape(Pid, size(T));
%  % sliding window averaging
%  window=round(tf/20);
%  Pid_s=zeros(size(T));
% for i=1:length(T)
%    Pid_s(i)=mean(Pid(max(i-window,1):min(i+window,length(T))));
% end
% Pid=Pid_s;
% Lambda1=log(N*r/2/sqrt(pi));
% Lambda1=log(N*r/2/sqrt(pi*Lambda1)); % one iteration
% kappa=0.5*r^2/s0^3/(2*Lambda1)^2; 
% betap=r/s0*sqrt(2*Lambda1/L); 
% 
% subplot(2,1,1)
%  %plot(T,Pid3/ s0,'r',T,Pid5/ s0,'b',T,Pid/ s0,'k')
%  plot( T,Pid/s0,'k')
% xlabel('Time, t');
% ylabel('Coal. prob., P_{id}/s0')
% title(sprintf('betap=%g kappa=%g Lambda1=%g sl. window=%g',betap,kappa, Lambda1, window))
% grid
% 
% subplot(2,1,2) 
% X=Pid/s0;
% Y=(2*betap)^2*(-log(X/kappa)).^(-2);
% plot(X,Y,'b',X,w2av/L,'r')
% xlabel('P_{id}/s0, simulation');
% ylabel('w_2/L');
% title(' w2_{anal} blue,  w2_{av} red')

% saving figs
figure(1)
saveas(gcf,strcat(homedir,filename,'_1.fig'));
figure(2)
saveas(gcf,strcat(homedir,filename,'_2.fig'));
figure(4)
saveas(gcf,strcat(homedir,filename,'_4.fig'));
figure(1)
% figure(3)
% saveas(gcf,strcat(homedir,filename,'_3.fig'));
% figure(4)
% saveas(gcf,strcat(homedir,filename,'_4.fig'));



end




