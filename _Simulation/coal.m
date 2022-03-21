%% histogram for coalescent density

global ts Nlineage T tcoal

N=1000; % population size 
mruns=10;
Nlin =[]; tco =[];

%% Loop in runs
for run=1:mruns
    recomb_2021('const',1,0.1,1,3,10,N,150,0.02,run)
%   recomb_2021(distribution_s,r,s0,a,M,L,N,tf,f0,run)
    Nlin=[Nlin;Nlineage];   % adding rows over runs
    tco=[tco,tcoal];        % growing string of coal.time points
end
Nlin=mean(Nlin);

%% The coal.density  
[hi,xx]=hist(tco,round(sqrt(length(tco))));  
area=(xx(2)-xx(1))*sum(hi); % normalization
nfit= Nlin(end)-Nlin(1);
invNeff_num =(nfit*hi/area); 
% normalizing to neutral and subtr. 1
intN=interp1(T,Nlin,xx);
invNeff_num = N*invNeff_num./(intN.*(intN-1)/2) - 1;
% mean and std
meant = sum(xx.*invNeff_num)/sum(invNeff_num);
stdt=sqrt(sum((xx-meant).^2 .*invNeff_num)/sum(invNeff_num));

%% Plot norm.coal.density
figure(10)
%plot(xx,invNeff_num,'ro-');
bar(xx,invNeff_num)
%  printing mean and std
ylabel('coal.den./neu. - 1')
xlabel('Time'); 
title( ts); box off
axi=axis; axi(3)=0; axis(axi)
text((axi(1)+axi(2))/2,(axi(3)+axi(4))/2,sprintf('meant=%.2g stdt=%.2g',meant,stdt))