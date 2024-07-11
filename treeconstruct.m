%This function calculates parivise distances between extant species and
%constructs a phylogenetic tree. 

%Further construct phylocom sample file

function [phylotree, newickstr, phylocom_sample]=treeconstruct(Output, mergee,time,nHab,run_string)

%First calculate parivise distance
sp=unique(mergee(:,1));

for i=1:length(sp)
    %List ansestry history
    sp1=sp(i);
    nodehist1=Output(sp1,1);
    while nodehist1(end)>0
        nodehist1(end+1)=Output(nodehist1(end),1);
    end        
    
    for j=1:length(sp)
        %List ansestry history
        sp2=sp(j);
        nodehist2=Output(sp2,1);
        while nodehist2(end)>0
            nodehist2(end+1)=Output(nodehist2(end),1);
        end      
        
        %Find earlyest common node for the two species
        common_node=max(intersect(nodehist1,nodehist2));
        node_index=find(Output(:,1)==common_node);
        dist=time-Output(node_index(1),2);
        
        %dist=(time-Output(sp(i),2))+(time-Output(sp(j),2));
        
        PVDist(i,j)=dist;
    end
end


%Create name vector for leafs in the tree
name={};
namebase1 = 'species';
for i=1:length(sp)
    namebase2 = num2str(sp(i));
    namebase=[namebase1 namebase2];  
    name(i)={num2str(namebase)}; 
end


%Create tree and save newick string in txt file
phylotree=seqlinkage(PVDist,'average',name);

newickstr=getnewickstr(phylotree); %extract newick string from tree
% fil=fopen(['newicktree' run_string],'w');          
% fprintf(fil,newickstr);
% fclose(fil);


%Create phylocom sample file
%Three columns with habitat, abundance, species_id


% %Write to file and save file
% phylocom_sample=[];
% 
% hab=1:nHab;
% 
% fil=fopen(['phylosample' run_string],'w'); %open file for writing
% for i=hab
%     nameHab=['habitat' num2str(i)]; %define name of habitat 
%     
%     %Find species with abundance>1 in habitat
%     tmp=find(mergee(:,3+i)>0);
%     sp_id=unique(mergee(tmp,1));
%     
%     %Calculate abundance for each of the species and write to file
%     for j=sp_id'
%         tmp=find(mergee(:,1)==j);
%         abun=sum(mergee(tmp,3+i));
%         
%         nameSpec=['species' num2str(j)];
%         
%         fprintf(fil,'%s\t%i\t%s\t\n',nameHab,abun,nameSpec);
%     end
% end
% fclose(fil);



%Creat cellarray for the info
phylocom_sample={};

hab=1:nHab; 

for i=hab
    nameHab=['habitat' num2str(i)]; %define name of habitat 
    
    %Find species with abundance>10 in habitat
    tmp=find(mergee(:,3+i)>10);
    sp_id=unique(mergee(tmp,1));
    
    %Calculate abundance for each of the species save in phylocom_sample
    for j=sp_id'
        tmp=find(mergee(:,1)==j);
        abun=sum(mergee(tmp,3+i));
        
        nameSpec=['species' num2str(j)];
        
        phylocom_sample(end+1,1)={nameHab};
        phylocom_sample(end,2)={abun};
        phylocom_sample(end,3)={nameSpec};
        
    end
end
    


    