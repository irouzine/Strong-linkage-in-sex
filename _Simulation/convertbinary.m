% Takes matrix of 0 and 1 and converts to A and C for alignment and
% phylogeny and prints to a text file for MEGA software
global Ksample

for run=1:3
recomb_2021('const',1,0.1,1,3,200, 1000,20,0.02,run)

x=size(Ksample); L=x(2);n=x(1);
out=zeros(1,L);
FileID=fopen(sprintf('/Users/igorrouzine/Desktop/Recombination/_Simulation/sequenceprint%g.txt',run),'w');
% FileID1=fopen('/Users/igorrouzine/Desktop/Recombination/_Simulation/sequenceprint1.txt','w');
% FileID2=fopen('/Users/igorrouzine/Desktop/Recombination/_Simulation/sequenceprint2.txt','w');
% FileID3=fopen('/Users/igorrouzine/Desktop/Recombination/_Simulation/sequenceprint3.txt','w');

for i=1:n
    for j=1:L 
        if Ksample(i,j)
            out(j) = 'A';
        else
            out(j) = 'C';
        end
    end

    fprintf(FileID,'%s\n',out);  % full genome
%     fprintf(FileID1,'%s\n',out(1:round(L/3)));  % first third
%     fprintf(FileID2,'%s\n',out(round(L/3)+1:round(2*L/3)));  % middle third
%     fprintf(FileID3,'%s\n',out(round(2*L/3)+1:end));  % last third
end

fclose(FileID); %fclose(FileID1);fclose(FileID2);fclose(FileID3);
    mean(mean(Ksample))
end %in runs