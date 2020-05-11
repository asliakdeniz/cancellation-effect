%Shift Process

%Parameter Values
iM              = 50;                   %number of groups
iN              = 15;                   %group size
iT              = 3000000;              %maximum number of time periods
dR              = 0.1;                  %probability of a migration event
dP              = (1-dR)*iN/(iN+1);     %probability of an individual selection event
dQ              = (1-dR)*1/(iN+1);      %probability of a group selection event
dWi             = 0.1;                  %intensity of individual selection
dWg             = 0.1;                  %intensity of group selection
iRepl           = 1000000;              %number of replications
dB              = 0.095;                %benefit of cooperation
dC              = 0.1;                  %cost of cooperation
a               = 0;                    %to check whether any run is finished before reaching either fixation of extinction
error           = 0;                    %to check whether two versions of population structure is the same
   
rng(40);                                %setting the seed

vPayoff         = zeros(iM,1);          %group payoff vector
vCoopFinal      = zeros(iRepl,1);       %final fraction of cooperators
vCoopFinalM     = zeros(iRepl,1);       %final fraction of cooperators (second version)
iFc             = 1-dWi*dC;             %payoff of a Cooperator within the group
iFd             = 1;                    %payoff of a Defector within the group

for     i=1:iRepl       %simulation loop
    i
vGroups         =   zeros(iM,1);    %vector with no. of cooperators in each group at a point in time
vGroups(1,1)    =   1;              %initial mutation
mGroups         =   zeros(iM, 2);   %matrix with no. of cooperators in each group at a point in time 
mGroups(1,1)    =   1;              %initial mutation
mGroups(:,2)    =   mGroups(:,1);
vCoop           =   zeros(iT,1);    %vector to store the total no. of cooperators at each time step
vCoop(1,1)      =   sum(vGroups(:,1));
vCoopM          =   zeros(iT,1);    %vector to store the total no. of cooperators at each time step (second version)
vCoopM(1,1)     =   sum(mGroups(:,1));

t               =   2;
while    t <=iT && vCoop(t-1)>0 
    
    dRand   = rand;
    if mGroups(:,1) == vGroups(:,1)     %check whether the two versions match    
        if dRand < dP   %individual selection event
            iIndexInd  = ceil(iM*rand);     %randomly selecting a group 
            iCo     = mGroups(iIndexInd,1); 
            dBirthC = iCo*iFc/(iCo*iFc+(iN-iCo)*iFd);   %birth rate of cooperators
            dBirthD = (iN-iCo)*iFd/(iCo*iFc+(iN-iCo)*iFd);  %birth rate of defectors
    
            dProb   = rand;

            if dProb < dBirthC*(iN-iCo)/iN     %probability that a cooperator replaces a defector
                vGroups(iIndexInd,1)    = iCo+1;
                mGroups(iIndexInd,2)    = iCo+1;
            elseif  dProb < (dBirthC*((iN-iCo)/iN)+dBirthD*(iCo/iN))    %probability that a defector replaces a cooperator
                vGroups(iIndexInd,1)    = iCo-1;
                mGroups(iIndexInd,2)    = iCo-1;
            end

        elseif dRand < dP+dQ    %group selection event
            vCo         = mGroups(:,1);    
            vPayoff     = ones(iM,1)+dWg*(vCo.*(dB/iN));    %group payoffs
            vP          = vPayoff/sum(vPayoff);     %relative group payoffs (group reproduction probabilities)

            dProbGr       = rand;
            iIndexGr      = 0;

            while dProbGr > 0     %while loop to pick a group for reproduction
                 iIndexGr  = iIndexGr + 1;
                 dProbGr   = dProbGr - vP(iIndexGr);    
            end
            
            iRandom  = ceil(rand*(iM));     %randomly selecting a group to die 
            D = rand;
            if iRandom < iIndexGr
                if D<0.5
                    vGroups(iRandom:iIndexGr-1,1) = vGroups(iRandom+1:iIndexGr,1);
                    mGroups(iRandom:iIndexGr-1,2) = mGroups(iRandom+1:iIndexGr,1);
                else 
                    vGroups(2:iRandom,1) = vGroups(1:iRandom-1,1); 
                    vGroups(1,1) = vGroups(iM,1);
                    vGroups(iIndexGr+1:iM,1) = vGroups(iIndexGr:iM-1,1);
                    mGroups(2:iRandom,2) = mGroups(1:iRandom-1,1); 
                    mGroups(1,2) = mGroups(iM,1);
                    mGroups(iIndexGr+1:iM,2) = mGroups(iIndexGr:iM-1,1);   
                end
            elseif iRandom > iIndexGr
                if D>=0.5 
                    vGroups(iIndexGr+1:iRandom,1) = vGroups(iIndexGr:iRandom-1,1);
                    mGroups(iIndexGr+1:iRandom,2) = mGroups(iIndexGr:iRandom-1,1);
                else 
                    vGroups(iRandom:iM-1,1) = vGroups(iRandom+1:iM,1); 
                    vGroups(iM,1) = vGroups(1,1);
                    vGroups(1:iIndexGr-1,1) = vGroups(2:iIndexGr,1);
                    mGroups(iRandom:iM-1,2) = mGroups(iRandom+1:iM,1); 
                    mGroups(iM,2) = mGroups(1,1);
                    mGroups(1:iIndexGr-1,2) = mGroups(2:iIndexGr,1);   
                end
            else    %iRandom == iIndexGr
                vGroups(:,1) = vGroups(:,1);
                mGroups(:,2) = mGroups(:,1);
            end 
        
        else    %dRand<1 %migration event
            iIndexMig  = ceil(iM*rand);  %randomly selecting a group
            dProbMig = rand; %cumulative prob. for the migration -- not used in the current version of migration
            dDirection = rand;  %direction of migration
        
            iCoMig = mGroups(iIndexMig,1);  %no. of coop in migrating group
            iX = ceil(iN*rand);     %randomly picking the migrating player
            iY = ceil(iN*rand);     %randomly picking the migrating player
            if dDirection<0.5   %migration to the right
                if iIndexMig==iM
                    iPartner = 1;
                else
                    iPartner = iIndexMig+1;
                end
                iCoRec = mGroups(iPartner,1);   %no. of coop in the partner group
                if iIndexMig~=iPartner
                    if (iX<=iCoMig) && (iY>iCoRec)  %migrating gr. sends a cooperator & partner gr. sends a defector
                        mGroups(iIndexMig,2)=iCoMig-1;           
                        vGroups(iIndexMig,1)=iCoMig-1;           
                        mGroups(iPartner,2)=iCoRec+1;            
                        vGroups(iPartner,1)=iCoRec+1;            
                    elseif (iX>iCoMig) && (iY<=iCoRec)  %migrating gr. sends a defector & partner gr. sends a cooperator
                        mGroups(iIndexMig,2)=iCoMig+1;       
                        vGroups(iIndexMig,1)=iCoMig+1;       
                        mGroups(iPartner,2)=iCoRec-1;        
                        vGroups(iPartner,1)=iCoRec-1;        
                    end
                end
            else    %migration to the left
                if iIndexMig==1
                    iPartner = iM;
                else
                    iPartner = iIndexMig-1;
                end
                iCoRec = mGroups(iPartner,1);   %no. of coop in the partner group
                if iIndexMig~=iPartner
                    if (iX<=iCoMig) && (iY>iCoRec)  %migrating gr. sends a cooperator & partner gr. sends a defector
                        mGroups(iIndexMig,2)=iCoMig-1;           
                        vGroups(iIndexMig,1)=iCoMig-1;           
                        mGroups(iPartner,2)=iCoRec+1;            
                        vGroups(iPartner,1)=iCoRec+1;            
                    elseif (iX>iCoMig) && (iY<=iCoRec)  %migrating gr. sends a defector & partner gr. sends a cooperator
                        mGroups(iIndexMig,2)=iCoMig+1;       
                        vGroups(iIndexMig,1)=iCoMig+1;       
                        mGroups(iPartner,2)=iCoRec-1;        
                        vGroups(iPartner,1)=iCoRec-1;        
                    end
                end
            end
        end         
    else
        error=error+1;
    end

vCoop(t)=sum(vGroups(:,1));    %no. of cooperators at t
vCoopM(t) = sum(mGroups(:,2));  %no. of cooperators at t (second version) 
mGroups(:,1) = mGroups(:,2);
t   = t+1; 
if vCoop(t-1)==iM*iN
    t=iT;
   vCoop(iT)=iM*iN; 
   vCoopM(iT)=iM*iN; 
end
end
vCoopFinal(i)=vCoop(iT)/(iM*iN);    %final fraction of cooperators at run i
vCoopFinalM(i) = vCoopM(iT)/(iM*iN);    %final fraction of cooperators at run i (second version)
if (vCoopFinal(i)<1) && (vCoopFinal(i)>0)   %check if the final fraction is either 0 (extinction) or 1 (fixation)
    a = a+1;
end
end
ShiftSum = sum(vCoopFinalM);    %final no. of cases where fixation occured
dShift = sum(vCoopFinal)/iRepl;     %final fraction of cases where fixation occured
vF = [iM; iN; dR; dB; iT; dWi; dWg; dShift; ShiftSum; a; error];    %vector containing model parameters and outcome variables
