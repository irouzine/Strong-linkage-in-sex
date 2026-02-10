function recomb_2021(distribution_s,r,s0,a,M,L,N,tf,f0,run)

% file locations 
homedir = '~/Desktop/Recombination/figs';
% file titles for saving figures
filename = sprintf('_%s_%g_%g_%g_%g_%g_%g_%g_%g_%g',...
    distribution_s,r,s0,a, M,L,N,tf,f0,run);

% Output used in external programs
global ts Nlineage T tcoal Ksample wsample

%% Arguments:
% distribution_s is one of the four forms below
% r recombination rate per genome
% s0 average selection coefficient
% a index in the power distribution_s (not used)
% L number of loci
% N population number
% M crossover number
% tf full time in genertions
% run the number of the run

%% random seed, so can average over runs.
rng(run)

yessave=0; % save figure or not?
yesplot=0; % plotting anything but Fig 1?

%% Mutations are optional and entered here:
NUbs2=0;                      % N*Ub*s^2
mu=NUbs2/N/s0^2/L;
%%

%mu=0;% 3e-5;                    % recurrent beneficial mutation rate 
tesc=0;                         % time interval for a kind of initial conditions
T= 0:tf ;                       % times
N_sample=30;                    % Sampling segment for w2 and LD
muesc=0;                        % standing variation f0(s): deleterious mutation rate during escape

%% Distributions of s and parameter titles of fig 1 

switch distribution_s
    case 'power' 
        % 1: a*s0^a/(s+s0)^(a+1) 
        s=s0*(rand(1,L).^(-1/a)-1); 
        ts=sprintf('1/s^{a+1},s>s0. N=%g,r=%g,L=%g,s0=%g,\n a=%g,f0=%g,T=%g,M=%g',N,r,L,s0,a,f0,max(T),M);
    case 'lognormal'
        % 1/s*exp[-log^2(s/s0)/2/a^2]
        s=s0*exp(a*randn(1,L));
        ts=sprintf('1/s*exp[-log^2(s/s0)/2/a^2]\n N=%g,r=%g,L=%g,s0=%g,a=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n T=%g,M=%g',...
            N,r,L,s0,a,f0, muesc,mu,tesc,max(T),M);
    case 'exponential'
        % exp(-s/s0) %% actually density is, g(s) = (1/s0)*exp(-s/s0)
        s=-s0*log(rand(1,L));
        ts=sprintf('exp(-s/s0)\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n T=%g,M=%g',...
            N,r,L,s0,f0, muesc,mu,tesc,max(T),M);
    case 'const'
        s = s0*ones(1,L);
        ts=sprintf('const, N=%g,r=%g,L=%g,s0=%g\n f0=%g, NUbs2=%g,T=%g, M=%g, run=%g',...
            N,r,L,s0,f0, NUbs2,max(T),M, run); % no mutation
end


%% Initial settings
 
fsite=zeros(length(T),L); Plinkav=fsite; %fanc=fsite;

fav=zeros(length(T),1);  Cav=fav; w2av=fav; 
nlineages=fav; Closs=fav; fpol=fav; Cpol=fav; Vark=fav;
w2=zeros(1,N_sample); Plink=zeros(N_sample,L);  C=zeros(1,L); 
nprog=zeros(N,1);


figure(2);clf
 
Knew=zeros(N,L);  % binary DNA sequences
Anew=zeros(N,L);  % ancestor labels
Pnew=zeros(N,L);  % parent labels

tint=tf/3; % plot time interval 
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

A=(1:N)'*ones(1,L);  % Initial ancestor labels used for C(t)

P1=zeros(N,length(T)); Pm=P1; PL=P1; % Initial parent labels for 3 sites in time, for phylogeny 
W=P1;                                % Initial fitness values
Sig2=zeros(length(T),L);

%% Evolution starts...
for t=T

% beneficial mutation if any
    if mu > 0
        K = K | sprand(N,L,mu);
    end
    
% Initial parent labels for one-time population 
    P=(1:N)'*ones(1,L);  
   
%% Recombination of randomly chosen pairs with one-parent replacement
    npairs=round(r*N/2);
    ii=ceil(rand(npairs,2)*N); i1=ii(:,1); i2=ii(:,2);  % 2 columns of random indices of parents
    for i=1:npairs
        % generating random 1 or 0 for each site, with probabilities M/L and 1-M/L,
        % respectively, to mark crossovers with 1
        % Even/odd xx shows site segments copied from 1st parent, 2nd, 1st, 2nd etc
        xx=cumsum(rand(1,L) < M/L); 
        first=(round(xx/2)==xx/2);                          %  sites copied from 1st parent are marked as 1
        prog1=K(i1(i),:).*first+K(i2(i),:).*(1-first);      % recombinant DNA sequence
        progA=A(i1(i),:).*first+A(i2(i),:).*(1-first);      % recombinant ancestor labels
        progP=P(i1(i),:).*first+P(i2(i),:).*(1-first);      % recombinant parent labels
        K(i1(i),:)=prog1;                                   % 1st parent's DNA replaced   
        A(i1(i),:)=progA;                                   % 1st parent's ancestor labels replaced
        P(i1(i),:)=progP;                                   % 1st parent's parent labels replaced
    end
   
    
%% Random sampling of progeny and natural selection with a broken stick

% column of N log-fitnesses ; state with 0s only has w=0 (fitness 1) by definition
    w=K*(s'); 
% add epistasis here! Input parameters Tij, Eij are set separately in the beginning.
    
    nprogav=exp(w)/mean(exp(w));    % average progeny number
    b2=cumsum(nprogav);
    b1=[0;b2(1:end-1)];             % broken stick
    X=rand(N,1)*N;
    for i=1:N
        nprog(i)=sum( X > b1(i) & X < b2(i)); % actual progeny number
    end
  
%% Updating population
    is=[0;cumsum(nprog(1:(N-1)))];
    for i=1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*K(i,:);  % DNA sequences
            Anew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*A(i,:);  % ancestor labels
            Pnew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*P(i,:);  % parent labels
        end
    end
    K=Knew;
    A=Anew; 
    P=Pnew;
    sK=size(Knew);sK=sK(1);
    if sK  ~= N, disp('N not conserved'),return;end
    
    % updating sequences done
    
    %% Memorizing observables 
    W(:,t+1)=K*(s');                % fitness column
    if tint*round(t/tint)==t
    for jj=1:L
        Sig2(t+1,jj)=std(K(:,1:jj)*s(1:jj)')^2;          % Variation of segment fitness 
    end
    Sig2(t+1,:)=Sig2(t+1,:)/Sig2(t+1,L);   % normalized variation of segment fitness
    end
    P1(:,t+1)=P(:,1);               % parent labels for 3 sites  
    Pm(:,t+1)=P(:,L/2); 
    PL(:,t+1)=P(:,L);                 
    ipol=find(~all(1-K,1));         % non-lost allele site  indices
    Closs(t+1)=1-length(ipol)/L;    % fraction of sites with extinct (or non heterozygous) alleles
    fsite(t+1,:)=mean(K);           % 1-site allele frequencies at all sites
    fav(t+1)=mean(mean(K));         % average allele frequency per site per genome
    fpol(t+1)=mean(mean(K(:,ipol)));% average allele frequency per polymorphic site 
    Vark(t+1)=(std(w)/s0)^2;        % variance of the allele number between genomes
    
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
    Plinkav(t+1,:)= mean(Plink);  
    %allowed=fsite(t+1,:) > 0 && fsite(t+1,:) < 1; % binary mask of polymorphous sites
     
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
    if tint*round(t/tint)==t && yesplot
        c=col(round(t/tint)+1);
        figure(2)
     subplot( tf/tint +2,1,1)
        [nn,xx]=hist( w );                % histogram of fitness among genomes
        semilogy(xx,nn,[c '+'])
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
    
    %% Plotting LD (average Pearson r^2) at some time points
    if tint*round(t/tint)==t && yesplot
        c=col(round(t/tint)+1);
        figure(3)
        LD=zeros(1,L );
        for j=0:(L-1)
            yy=0; nsum=0;
            for i=1:(L-j)
                xx1=fsite(t+1,i)*(1-fsite(t+1,i));
                xx2=fsite(t+1,i+j)*(1-fsite(t+1,i+j));
                % if each site is more than 5% diverse, we add LD
                if xx1 > 0.05 && xx2 > 0.05  
                    nsum=nsum+1;
                    yy=yy+(mean(K(:,i).*K(:,i+j)) - fsite(t+1,i)*fsite(t+1,i+j))^2/xx1/xx2;
                end
            end
            LD(j+1)=yy/nsum;
        end
        plot(0:(L-1),LD,c)
        hold on
    end
    
        %% Testing block independence at some points
    if tint*round(t/tint)==t && yesplot
        c=col(round(t/tint)+1);
        figure(10)
        plot(1:L,Sig2(t+1,:),c)
        iii=round(t/tf*(L-1)+1);
        text(iii,Sig2(t+1,iii),sprintf('t=%g',t))
        xlabel('Segment length')
        ylabel('Fitness variance normalized')
        hold on
    end
end

% Evolution has ended

%% Final plots

if yesplot 
    figure(10);hold off
end

if yesplot
figure(2) % Traveling wave
subplot( tf/tint +2,1,1)
hold off
xlabel('Log fitness,  w  ');
ylabel('Density distribution')
title(ts)   % title  with parameter values
 subplot( tf/tint +2,2,2*tf/tint +3)
 xlabel('Family size/N ');
  subplot( tf/tint +2,2,2*tf/tint +4)
 xlabel('Max family size/N ');

%% 

figure(3)  % LD
hold off
xlabel('Locus distance');
ylabel('Average r^2')
title(ts)   % title  with parameter values
end
 
%% figure 1
figure(1)
% average ovservables vs time
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

% adaptation rate
v=(mean(W(:,end))-mean(W(:,round(end/2))))/(T(end)-T(round(end/2)));

subplot(2,1,1)
plot(T,fav,'b',T,fpol,'b+',T,f1site,'b--',T,10*w2av/L,'g',...
    T,Cav,'m',T,Cpol,'m+',T,1-Closs,'r',T,max(fsite'),'k');
title(sprintf('f_{av} b   f_{pol} b+   f_{1site} b- -   10*w2av/L g   C m   C_{pol} m+ \n v=%g r max(fsite) bk',v))
xlabel('Time, t');
ylabel(sprintf('All.frequency f f_{pol}, half-distance w^2 \n allele loss C_{loss}, id.descent C C_{pol}'))
axi=axis;axi(3:4)=[0 1 ];axis(axi);
box off
text(tf/2,0.5,ts)


%% Tests of relationship between width and C and Fisher
% subplot(2,1,2)
% plot(T, Vark./w2av, 'co-',T, Fisher, 'c+',...
%     T, w2av./fpol./(1-fpol)./(1-Closs) /L, 'k+', T, w2av./fav./(1-fav)./(1-Cav)/L, 'ko',T, w2av./fav./(1-fav-Closs)./(1-Cav).*(1-Closs)/L, 'k')
% title(sprintf('Var[k]/w^2  co-   V/Vark/s0 c+   w^2/f_{pol}/(1-f_{pol})/(1-C_{loss})/L  k+ \n  w^2/f_{av}/(1-f_{av})/(1-C_{av})/L  ko  w^2_{av}(1-C_{loss})/f_{av}/(1-f_{av}-C_{loss})/(1-C_{av})/L k-'))
% xlabel('Time, t');
% ylabel(sprintf('Normalized w^2 and Var[k]'))
% axi=axis;axi(3:4)=[0 2 ];axis(axi);
% grid

%% Histograms fsite vs time

if yesplot
figure(4)

i=0;
for t=0:tf/5:tf
    i=i+1;
    %plot(t*ones(1,L),fsite(t+1,:), [c(i+1),'o'])

    ii=find(fsite(t+1,:) > 0 & fsite(t+1,:) < 1);
%     [nn,xx]=hist(fsite(t+1,ii));                 
%     plot(xx,nn,c(i))'

subplot(7,1,i+1)
    histogram(fsite(t+1,ii),'FaceColor',col(i));
    %ir=ceil( length(nn)*rand);
    title( sprintf('t=%g',t))
    axi=axis; axi(1:2)=[0 1 ]; axis(axi);
    box off
end  
xlabel('Allelic frequency for locus, f_{loc} ');
%ylabel('# of sites')

% One site allele frequncies vs time
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
 box off
end %yesplot

%% Phylogeny plots

nsample=1:N/10:N; %  genome indices in the sample below

if yesplot  % whether to do plots apart from the basic dynamics in Fig. 1
    
nlineage1=nsample; % initial lineage labels
nlineagem=nsample; 
nlineageL=nsample;
m=length(nsample);
yesconnect1=ones(size(nsample)); % segment plotting markers
yesconnectm=ones(size(nsample));
yesconnectL=ones(size(nsample));

Nlineage=zeros(1,length(T)); % initial total lineage number at each time

% colors of lineages
co='krbmgcy';
for i=1:8
    co=[co,co];
end
tcoal=[]; % vector of coalescent times

% Loop back in time 
for t=tf:-1:0

    
%%  fitness trajectory  
    figure(7) 
    for i=1:m
        if t > 0
            plot([t t+1],[W(Pm(nlineagem(i),t+1),t),W(nlineagem(i),t+1)],co(i))
            hold on
        end
    end
   
%% trajectory for middle site
    figure(6)
    subplot(4,1,1)
    for i=1:m
        plot([t t+1],[Pm(nlineagem(i),t+1),nlineagem(i)],co(i))
        hold on
    end

%%  3 site trees
    % first site tree
    subplot(4,1,2)
    for i=1:m 
        % find all intersections upward
        jj=find(P1(nlineage1(i),t+1)==P1(nlineage1(i+1:end),t+1));
        if isempty(jj) && yesconnect1(i)
            % if none, plot the horizontal segment
            plot ([t t+1],[i i],co(i)) 
        elseif yesconnect1(i)
            % plot a segment angled up
            plot ([t t+1],[max(jj)+i, i],co(i)) 
            % do not plot the lineage further
            yesconnect1(i)=0; 
            % memorize the time point
            tcoal=[tcoal,t];  
        end 
        hold on  
    end
    % recursive trajectory step
    nlineage1=P1(nlineage1,t+1)'; 
    
    %  middle site tree 
    subplot(4,1,3)
    for i=1:m 
        jj=find(Pm(nlineagem(i),t+1)==Pm(nlineagem(i+1:end),t+1));
        if isempty(jj) && yesconnectm(i)
            plot ([t t+1],[i i],co(i)) 
        elseif yesconnectm(i)
            plot ([t t+1],[max(jj)+i, i],co(i)) 
            yesconnectm(i)=0; 
             tcoal=[tcoal,t];  
        end 
        hold on  
    end
    nlineagem=Pm(nlineagem,t+1)'; 
    
    %  last site tree 
    subplot(4,1,4)
    for i=1:m 
        jj=find(PL(nlineageL(i),t+1)==PL(nlineageL(i+1:end),t+1));
        if isempty(jj) && yesconnectL(i)
            plot ([t t+1],[i i],co(i)) 
        elseif yesconnectL(i)
            plot ([t t+1],[max(jj)+i, i],co(i)) 
            yesconnectL(i)=0; 
             tcoal=[tcoal,t];  
        end 
        hold on  
    end
    nlineageL=PL(nlineageL,t+1)';  
 
 % Number of unique lineages averaged over 3 loci
    Nlineage(t+1)=mean([length(unique(nlineage1)) length(unique(nlineagem)) length(unique(nlineageL))]); 
end % Loop in time ends

figure(7);hold off
title(sprintf(['Fitness trajectory middle locus \n' ts])) 
axi=axis; axi(1:2)=[0 tf+1]; 
axis(axi); box off
ylabel ('Fitness, w');
xlabel('Time')

figure(6)
subplot(4,1,1); hold off; 
title(sprintf(['Trajectory middle locus \n' ts])) 
axis([0 tf+1 0 N]); box off
ylabel ('Individual genome');

subplot(4,1,2); hold off; 
title('Phylogeny of first locus'); box off
axis([0 tf+1 0 length(nsample)+1])
 ylabel('Lineages')
 
subplot(4,1,3); hold off; 
title('Phylogeny of middle locus'); box off
axis([0 tf+1 0 length(nsample)+1])
 ylabel('Lineages')
 
subplot(4,1,4); hold off; 
title('Phylogeny of last locus'); box off
axis([0 tf+1 0 length(nsample)+1])
xlabel('Time'); ylabel('Lineages')


end %yesplot


%% Memorize the last sequences of the sample for MEGA inference and their fitness values
Ksample=K(nsample,:);
wsample=w(nsample);


%% obsolete coal density
% Number of lineages and coal. density vs time 
% subplot(5,1,5) 
% %sfit=0.05; 
% %t50=80; 
% nfit= Nlineage(end)-1;
% 
% %  interpolating Nlineage vs time
% banana=@(x)sum((1 + nfit*(1+exp(-x(1)*(T-x(2)))).^(-1) - Nlineage).^2);
% xout = fminsearch(banana, [0.05 80]);
% sfit=xout(1);   
% t50=xout(2); 
% % Smooth density of coalescent events per lineage, dlogy/dt 
% % invNeff=nfit*sfit*(exp(sfit*(T-t50)/2)+exp(-sfit*(T-t50)/2)).^(-2) ./y;
% y =1 + nfit*(1+exp(-sfit*(T-t50))).^(-1); %interpolation of Nlineage
% invNeff=y(2:end)-y(1:(end-1)); 
% invNeff=[invNeff,invNeff(end)]; 
% % normalizing to neutral and subtr. 1
% invNeff=N*invNeff./(y.*(y-1)/2)-1; 
% 
% % The coal.density without smoothing
% [hi,xx]=hist(tcoal,round(sqrt(3*nfit))); % best number of bins 
% area=(xx(2)-xx(1))*sum(hi); % normalization  
% invNeff_num =(nfit*hi/area); 
% % normalizing to neutral and subtr. 1
% intN=interp1(T,Nlineage,xx);
% invNeff_num = N*invNeff_num./(intN.*(intN-1)/2)-1;
% % mean and std
% meant = sum(xx.*invNeff_num)/sum(invNeff_num);
% stdt=sqrt(sum((xx-meant).^2 .*invNeff_num)/sum(invNeff_num));
% 
% % Plot row and interpolated Nlineage and norm.coal.density
% plot(T,Nlineage*3,'b',T,y*3,'b',T,invNeff,'r',xx,invNeff_num,'ro');
% [yy,ii]=max(invNeff); %  printing best fit parameters
% text(T(ii),yy,sprintf('meant=%.2g stdt=%.2g',meant,stdt))
% ylabel('coal.den./neu. - 1')
% xlabel('Time'); 
% title('Mean lineage # x 3, blue; coal.den./neu. - 1, red'); box off
% axi=axis; axi(3)=0; axis(axi)




%% obsolete LD
% 
% figure(3) 
% 
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


if yessave 
% saving figs
figure(1)
saveas(gcf,strcat(homedir,filename,'_1.fig'));
figure(2)
saveas(gcf,strcat(homedir,filename,'_2.fig'));
figure(3)
saveas(gcf,strcat(homedir,filename,'_3.fig'));
figure(4)
saveas(gcf,strcat(homedir,filename,'_4.fig'));
figure(6)
saveas(gcf,strcat(homedir,filename,'_6.fig'));
figure(7)
saveas(gcf,strcat(homedir,filename,'_7.fig'));

end


end




