function p = setp

p.codetype = 'BCH';
p.SNRtype = 'SNR'; %'SNR'

p.N = 64;   %Block length
p.K = 51;   %Information length
p.D = 6;    %Dmin (not used in the algorithm, but used to identify matrix file)

p.m = 2;    %OSD decoding order
rate = p.K/p.N;

p.proset = 0.001;  % Promosing probablilty treshold 
p.sucset = 0.99;   % success probablilty treshold 

% SNR for similation
SNRrage = 3:0.5:8;

% Convert SNR to EbN0 if needed
switch p.SNRtype
    case 'SNR'
        p.SNR = SNRrage;
        p.ebN0 = p.SNR - 10*log10(rate) + 10*log10(1/2) ;
    case 'EbN0'
        p.ebN0 = SNRrage;
        p.SNR = p.ebN0 + 10*log10(rate) - 10*log10(1/2) ;
end


switch p.codetype
    case 'LDPC'
    [p.G] = AlistToMatrix(p.N,p.K,p.D);
    
    case 'BCH'

            Gm_name = ['genmat_' num2str(p.N) '_' num2str(p.K) '_' num2str(p.D) '.txt'];
            p.G = load(Gm_name);
            P = p.G(:,(p.K + 1) : p.N); 
            p.H = [P' eye(p.N - p.K)]; 
            
    case 'Polar'

%           Gm_name = ['Polar_' num2str(p.N) '_' num2str(p.K) '_0_Gen_CRC11.txt'];
            Gm_name = ['purePolarIRHARQ_GenM_n' num2str(p.N) '_k' num2str(p.K) '.txt'];
            p.G = load(Gm_name);
            p.G = GE(p.G);
            P = p.G(:,(p.K + 1) : p.N); 
            p.H = [P' eye(p.N - p.K)];  
      
end


fprintf('Running Simulation for (%d,%d,%d) %s code\n',p.N,p.K,p.D,p.codetype);
fprintf('The maxmum order of OSD decoding is: %d, \n ',p.m); 
fprintf('------------------------------------------\n');

rate = p.K/p.N;
loop = length(p.SNR);
p.PPV = zeros(1,loop);


%% Pre generation of TEP lists
for i = 1: loop
  
    snr = 10^(p.SNR(i)/10);
    %% PPV bound
    C = 0.5*log2(1+snr);
    V = (snr/2) * (snr+2)/((snr+1)^2) * log2(exp(1))^2;
    p.PPV(i) = qfunc(-(rate - C - 0.5*log2(p.N)/p.N)/sqrt(V/p.N));

end

   p.TEPlist = [];
    p.Tsz = 0;

    K = p.K;
    for order = 0:p.m  
        err = VChooseK(1:K,order);
        if order==0
            sz = 1;
            TEPpure = zeros(sz,K);
        else
            sz = Combination(K,order);   
            TEPpure = zeros(sz,K);
            for i = 1:sz
                TEPpure(i,err(i,:)) = 1;
            end
        end    
        p.TEPlist = [p.TEPlist; TEPpure];
        p.Tsz = p.Tsz+sz;
    end
    p.TEPlist = p.TEPlist(:,end:-1:1);

p.comSet = [];
for i = 0:round(p.K/4)
 p.comSet = [p.comSet nchoosek(p.K,i)];
end


end



