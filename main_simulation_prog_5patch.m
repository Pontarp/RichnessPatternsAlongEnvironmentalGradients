function lunarcsim(var1,var2,var3,var4,var5)

%var1 controls K (3 values)
%var2 controls sigma (3 values)
%var3 controls dispersal (3 values)
%var4 controls delta temp (5 values)
%var5 controls delta Kopt (5 values)



%Main program
%This program runs simulations of a metacommunity

%Itterate
iter=1:1;
for iii=iter
    
    
run_string=[num2str(var1) num2str(var2) num2str(var3) num2str(var4) num2str(var5) num2str(iter(iii))];
run_num=str2num(run_string);

%Set seeds for randfunctions
rand('seed',fix(sum(1000*clock)+var2))
randn('seed',fix(sum(2000*clock)+var2))

%Create info matrix for this run
    Output = [];
    invasion_analys=[];

%R�knare f�r respektive typ av mutationer 
raknazzuu=1;
raknazz=1;
raknauu=1;


%Parameters

sig_a_vector=[0.1 0.2 0.3];
sig_a= sig_a_vector(var2); %niche width for the resource trait
sig_u=0.2;  %niche width for the temperatur trait

r= 1; %intrinsic growth rate

K0_vector=[500 1500 2500];
K0= K0_vector(var1); %Maximium K-value

sig_K=1; %Width of the resource distribution
sig_U=1; %Width of the temp. distribution

%Different delta R_opt
Ropt_mat=[-0.2 -0.1 0 0.1 0.2;
      -1 -0.5 0 0.5 1;
      -2 -1 0 1 2;
      -3 -1.5 0 1.5 3;
      -4 -2 0 2 4]; 
Ropt=Ropt_mat(var5,:); %Resource peak postition in trait (z) space

%Different delta temp
Topt_mat=[-0.2 -0.1 0 0.1 0.2;
      -1 -0.5 0 0.5 1;
      -2 -1 0 1 2;
      -3 -1.5 0 1.5 3;
      -4 -2 0 2 4]; 
  Topt=Topt_mat(var4,:); %Temperature optimum in trait (u) space

nHab=5; %No. of habitats (must be same as length of Ropt)

P_mut_z = 1e-3; %Mutation probability
P_mut_u = 1e-3;
sig_mut = 0.02; %standard deviation of mutations

P_disp_vector=[1e-3 1e-2 1e-1]; %Dispersal probability
P_disp=P_disp_vector(var3);

%Dispersal regime according to matrix
P_disp_mat=[0 0.5 0 0 0; 
            1 0 0.5 0 0; 
            0 0.5 0 0.5 0;
            0 0 0.5 0 1;
            0 0 0 0.5 0];
P_disp_cum= cumsum(P_disp_mat);        


%Seed the system 
zz=Ropt(1); %Resource traits in a vector named zz
uu=Topt(1); %Resource traits in a vector named uu
NN= [10 0 0 0 0]; %Distribution of first species
spec_id=0; 



%Calculate Carying capacity for each of the zz values
KK= K0*exp(-(zz-Ropt).^2/2/sig_K^2);
%Calculate fitness loss because of temperature mismatch
UU= exp(-(uu-Topt).^2/2/sig_U^2);

time=100000;
%Loop over a number of time steps
for t=1:time

    %New vectors and matrices for the nex generation
    newNN = zeros(size(NN)); 
    newzz= zz; 
    newuu= uu; 
    newKK= KK; 
    newUU=UU; 

    %Calculate fitness for each lineage in each habitat (where it exist)
    
    for i=1:length(zz) %loop over zz vector
        fitness=zeros(1,nHab);
        for h=find(NN(i,:)>0)
            fitness(h)=fitfunc(zz(i),KK(i,h),UU(i,h),zz,NN(:,h),r,sig_a,sig_K);
        end
    
    
        %Calculate the number of offspring
        no_offspring= poissrnd(fitness.*NN(i,:));

        %Mutations:
        %No mut for z trait in the different habitats (output vector 1:length(no.hab))
        no_mutants_z=min(no_offspring,poissrnd(no_offspring*P_mut_z));
        %For u trait
        no_mutants_u=min(no_offspring,poissrnd(no_offspring*P_mut_u));

        %Output which habitat should mutate
        m= find(no_mutants_z>0);
        n= find(no_mutants_u>0);
        
        
        if length([m n])>0  %If mutation should occur in any habitat
            for j=unique([m n]) %loop over the habitat in which mutation will occure
                
                %Determine which offspring should mutate by randomly pulling a UNIQE number
                %between 1 and the number of offspring
                mutant_z=randperm(no_offspring(j));
                mutant_z=mutant_z(1:no_mutants_z(j));
                mutant_u=randperm(no_offspring(j));
                mutant_u=mutant_u(1:no_mutants_u(j));
                
                %Vector for number of mutants, the number vill be extracted
                %from the reproduction of the non mutating population
                mut_no_rep=zeros(1,nHab);
                

                %Find dubble mutatns
                tmp=intersect(mutant_z,mutant_u);
                
                %Create the dubble mutants
                for dumy=1:length(tmp)
                    newzz(end+1,1)=zz(i)+sig_mut*randn;
                    newuu(end+1,1)=uu(i)+sig_mut*randn;
                    newNN(end+1,j)=1; 
                    newKK(end+1,:)=K0*exp(-(newzz(end)-Ropt).^2/2/sig_K^2);
                    newUU(end+1,:)=exp(-(newuu(end)-Topt).^2/2/sig_U^2);
                    mut_no_rep(j)=mut_no_rep(j)+1;                    
                    raknazzuu=raknazzuu+1;
                    spec_id(end+1,1)=spec_id(i);
                end
                %Create single zz mutants 
                for dumy=1:length(mutant_z)-length(tmp)
                    newzz(end+1,1)=zz(i)+sig_mut*randn;
                    newuu(end+1,1)=uu(i);
                    newNN(end+1,j)=1; 
                    newKK(end+1,:)=K0*exp(-(newzz(end)-Ropt).^2/2/sig_K^2);
                    newUU(end+1,:)=UU(i,:);
                    mut_no_rep(j)=mut_no_rep(j)+1;
                    raknazz=raknazz+1;
                    spec_id(end+1,1)=spec_id(i);
                end
                %Create single uu mutants
                for dumy=1:length(mutant_u)-length(tmp)
                    newzz(end+1,1)=zz(i);
                    newuu(end+1,1)=uu(i)+sig_mut*randn;
                    newNN(end+1,j)=1; 
                    newKK(end+1,:)=KK(i,:);
                    newUU(end+1,:)=exp(-(newuu(end)-Topt).^2/2/sig_U^2);
                    mut_no_rep(j)=mut_no_rep(j)+1;
                    raknauu=raknauu+1;
                    spec_id(end+1,1)=spec_id(i);
                end


                newNN(i,:)=no_offspring-mut_no_rep;
            end%end of mutants per habitat and zz value
        else
            newNN(i,:)=no_offspring;
        end        
    end %end of zz loop for reproduction and mutation
    
     %Update infovectors
        NN=newNN; 
        KK=newKK;
        UU=newUU;
        zz=newzz;
        uu=newuu;
        
        %Remove extinct lineages
        alive = find(sum(NN,2)>0); 
        zz=zz(alive);
        uu=uu(alive);
        spec_id=spec_id(alive); 
        NN=NN(alive,:);
        KK=KK(alive,:);
        UU=UU(alive,:);
        
        
    %   Dispersal
    for i=1:length(zz)
        for h = find(NN(i,:)>0)
            number_dispersers = binornd(NN(i,h), P_disp);
            for d = 1:number_dispersers
                new_habitat = find(rand < P_disp_cum(:,h), 1, 'first');
                NN(i,h) = NN(i,h)-1;
                NN(i,new_habitat) = NN(i,new_habitat)+1;
            end
        end
    end
    
    %Cluster phenotypic values into clusters 
    
    if length(zz)>1 %No use to start this analysis if there is no variation
        if rem(t,100)==0 %Alos only do the analysis every 100 time steps
            [Output mergee invasion_analys]=speccheck(t,zz,uu,spec_id,NN,KK,UU,Output,invasion_analys);  
            spec_id=mergee(:,1); 
            zz=mergee(:,2); 
            uu=mergee(:,3);
            NN=mergee(:,4:8);
            KK=mergee(:,9:13);
            UU=mergee(:,14:18);
            
        end        
    end
    
      
        
    %Every now and then plot individuals in zz (y-axis) and uu (z-axis) trait space as a
    %function of time (x-axis)       
    if rem(t,100)==0
        doplot(t,zz,uu,NN,run_num)       
    end
   
    
    
end %end of time step loop



%Construct tree for end community
[phylotree, newickstr, phylocom_sample]=treeconstruct(Output, mergee, time, nHab, run_string);

%Save information for this run
saveas(run_num,['figure' run_string]);
%Save data
save(['outputdata' run_string],'Output','mergee','invasion_analys');
%Save tree and phylosample
save(['treedata' run_string],'phylotree','newickstr','phylocom_sample');
close

end
