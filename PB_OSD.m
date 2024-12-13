function [cEst,checkNum] = PB_OSD(r, p, noisePower)
    
    m = p.m; 
    N = p.N;
    K = p.K;
    G = p.G;
    comSet = p.comSet;
    Tsz = p.Tsz;
    TEPlist = p.TEPlist;
    N0 = 2 * noisePower;
    checkNum = 0;
    sucset = p.sucset;
    minDis = inf;
%% Oder statictic permutation

    a  = abs(r);
    EPini = 1./(1+exp(4*a/N0));
    
    NonEP = prod(1-EPini);
    
    [EPord, P1] = sort(EPini,'ascend');

    
%% Gaussian Elimination

    HD = (r<0);

    G1 = G(:,P1);
    [Gs,P2] = GE(G1);
    P = P1(P2);  
    iP = 1:N;
    iP(P) = 1:N;
    rRo = r(P);
    rRom = rRo(1:K);
    aRo = abs(rRo);
    aR0m = abs(rRom);
    oHD = (rRo <= 0);
    rHD = (rRom <= 0);
    

%% Probability base 

    EPordMRB = EPord(1:K);
    PPPord = sum(EPordMRB)/K;
    
    PbSuc = 0;
    for j = 0:m 
        PbSuc = PbSuc + comSet(j+1)*PPPord^j * (1-PPPord)^(K-j);     
    end  

    baseEt = prod(1 - EPordMRB);
    baseEp = NonEP/baseEt;

%% REPROCESSING
    centercode = mod(rHD*Gs,2);
    tempcenter = abs(centercode - oHD);
    minTEPWeight = inf;
    PbPro = p.proset* sqrt((1-PbSuc)/Tsz);
    EpLRB = mean(aRo(K+1:N));
    ErrPLRB = sum(EPord(K+1:N));
    M1 = ErrPLRB;
    V1 = ErrPLRB * (1 - ErrPLRB/(N-K));
    M2 = (N-K) * 0.5;
    V2 = (N-K) * 0.25;
    
    for i = 1:Tsz
        TEP = TEPlist(i,:);
        TEPdis = aR0m * TEP';
        
        if TEPdis > minTEPWeight ||  TEPdis > minDis
            continue;
        else
            Ptep = exp(-4*TEPdis/N0)*baseEt;
            WHDbound = minDis - TEPdis; 
            EstNum = WHDbound /EpLRB;  
            Promising = Ptep*(1 - 0.5 *erfc((EstNum-M1)/sqrt(2*V1)))  + (1-Ptep)*(1 - 0.5 *erfc((EstNum-M2)/sqrt(2*V2))); 
            if Promising < PbPro         
                minTEPWeight = min(minTEPWeight,TEPdis);
                continue;
            end
        end  
        corCode = rem(TEP*Gs,2);
        diPattern = abs(tempcenter - corCode);
        currentDis = diPattern*aRo';     % calculate the cerrent distance
        checkNum = checkNum + 1;     
        
        if (currentDis < minDis)            % if new codeword is closer, recode it as cerrent best codeword
            minDis = currentDis;
            currentBest = abs(corCode - centercode);  
            Ptep = exp(-4*TEPdis/N0)*baseEt;    
            Ppar = exp(-4*(currentDis - TEPdis)/N0)*baseEp;    
            Psuc = (Ptep.*Ppar)./(Ptep.*Ppar + (1-Ptep).*2^(K-N));  
            if  Psuc > sucset*PbSuc
                cEst = currentBest(iP);
                return;
            end
%             
        end
        
    end   
    
    cEst = currentBest(iP);

end