%%%% Simulation for OSD Decoding for BCH codes based on MATLAB.
%%%% By Chentao Yue, The University of Sydney, 
%%%% Jan, 2021
clear all;

%% 

%%% Include folder "matrices" in the path!!!

p = setp();   % Set up all coder/and decoder configurations in step();                            
i = 0;
loop = size(p.SNR,2);
canNum = zeros(1,loop);
Operations = zeros(1,loop);
WorstNum = zeros(1,loop);
Q = zeros(1,loop);
ebN0 = zeros(1,loop);
FER = zeros(1,loop);
BER = zeros(1,loop);
time = zeros(1,loop);
SCnum = zeros(1,loop);
yo = 0;
for snr = 1:loop
    tic;
    fprintf('SNR:%2.2f   ',p.SNR(snr));
    fprintf('ebN0:%2.2f   ',p.ebN0(snr));    
    rng(sum(100*clock));
    blockErr = 0;                                   % number of decoding block error
    blockNum = 0;                                   % total number of block
    berr = 0;   
    bitNum = 0;
    pos = [];
    while(blockErr <200 &&   blockNum <1e5)       

        %% Encoder
          mSg = randi([0,1],1,p.K);
%         mSg = zeros(1,p.K);
        [c,G] = Generator_Codewords(p, mSg);

        %% Modulation (BPSK)
        mC = (-1).^c;                 %BPSK

        %% Noise Channel
%       r = awgn(mC, p.SNR(snr) , 'measured');
        
        sigPower = sum(abs(mC.^2))/length(mC);
        sigPower = 10*log10(sigPower);
        noisePower = sigPower - p.SNR(snr);
        noisePower = 10^(noisePower/10);
        r = mC + sqrt(noisePower) * randn(1,p.N);
        
    
        %% Decoder 
              
        [cEst,checkNum] = PB_OSD(r, p, noisePower);  
        % cEst is the OSD output
        % checkNum is the number of TEPs used for this block

        %% count for TEPs
        canNum(snr) = canNum(snr) + checkNum;
        if checkNum > WorstNum(snr)
            WorstNum(snr) = checkNum;
        end     
        %% Error detection
        FLAG =  biterr(cEst,c);
        berr = berr + FLAG;
        if(FLAG ~= 0)
            blockErr = blockErr+1;  
        end    
        bitNum = bitNum + p.N;
        blockNum = blockNum + 1;  
    end
    
    BER(snr) = berr/bitNum;
    FER(snr) = blockErr/blockNum;
    canNum(snr) = canNum(snr)/blockNum;
    time(snr) = toc*1000/(blockNum);
    
    fprintf('Block Error:%d   ', blockErr);
    fprintf('Block Number:%7.d   ', blockNum);
    fprintf('BER:%.5f   ', BER(snr));
    fprintf('FER:%.5f   ', FER(snr));
    fprintf('Ave.TEP:%.3f    ', canNum(snr));
    fprintf('Wst.TEP:%.3f    ', WorstNum(snr));
    fprintf('Time(ms) per block:%.3f   \n ', time(snr));


end

save('Lastrun.mat')
fprintf('------------------------------------------\n');
fprintf('Simulation finished.');


Visualization(p.SNR, FER, 'SNR' , 'BLER' ,'decibel')
Visualization(p.SNR, canNum, 'SNR' , 'Ave. number of TEPs' ,'decibel')
Visualization(p.SNR, time, 'SNR' , 'Ave. time(ms) per block' ,'linear')
