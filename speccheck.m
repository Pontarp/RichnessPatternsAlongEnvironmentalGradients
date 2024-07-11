%This function checks the zz and uu vector for gaps in the traitvalues
%If a gap is found two new species are registered

function [Output mergee invasion_analys] = speccheck(t,zz,uu,spec_id,NN,KK,UU,Output,invasion_analys)

%mergee zz, uu and spec_id together
mergee=[spec_id zz uu NN,KK,UU]; 

%First check for gaps in zz
mergee=sortrows(mergee,[1 2]); %Sort according to species id and then according to zz-value 
for i=1:length(zz)-1    %loop over all zz values
    gap=abs(mergee(i,2)-mergee(i+1,2));  %calculate gaps between zz
    if gap>0.1  %if gap exist
        if mergee(i,1)==mergee(i+1,1)  %if gap is within one species
            origin=mergee(i,1);
            tmp=find(mergee(:,1)==mergee(i,1)); %find index for the species which will branch
            
            if sum(sum(mergee(min(tmp):i,4:8)))>20 & sum(sum(mergee(i+1:max(tmp),4:8)))>20 %if abundance in both new clusters is above some value
                mergee(min(tmp):i,1)=max(mergee(:,1))+1;
                mergee(i+1:max(tmp),1)=max(mergee(:,1))+1;
                
                Output(end+1:end+2,1)=origin;
                Output(end-1:end,2)=t;
                Output(end-1,3:7)=sum(mergee(min(tmp):i,4:8),1);
                Output(end,3:7)=sum(mergee(i+1:max(tmp),4:8),1);
                
%                 disp('zz branch')
%                 mergee
%                 Output
                
            end
        end
    end
end

%Now check for gaps in uu
mergee=sortrows(mergee,[1 3]);
for i=1:length(zz)-1
    gap=abs(mergee(i,3)-mergee(i+1,3));
    if gap>0.1
        if mergee(i,1)==mergee(i+1,1)
            origin=mergee(i,1);
            tmp=find(mergee(:,1)==mergee(i,1));
            
            if sum(sum(mergee(min(tmp):i,4:8)))>20 & sum(sum(mergee(i+1:max(tmp),4:8)))>20 %if abundance in both new clusters is above some value

                mergee(min(tmp):i,1)=max(mergee(:,1))+1;
                mergee(i+1:max(tmp),1)=max(mergee(:,1)+1);

               
                Output(end+1:end+2,1)=origin;
                Output(end-1:end,2)=t;
                Output(end-1,3:7)=sum(mergee(min(tmp):i,4:8),1);
                Output(end,3:7)=sum(mergee(i+1:max(tmp),4:8),1);
                
%                  disp('uu branch')
%                 mergee
%                 Output
                 
            end
        end
    end
end

%Sumarize the spatial distribution of species over time 

tmp=unique(mergee(:,1)); 
for i=tmp'
    tmp2=find(mergee(:,1)==i);
    invasion_analys(end+1,1)=t; 
    invasion_analys(end,2)=i;
    invasion_analys(end,3:7)=sum(mergee(tmp2',4:8),1);
end


