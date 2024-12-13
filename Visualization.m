function Visualization(SNR,BER,xl,yl,scale)

    figure;
    
    switch scale
        case 'decibel'
            semilogy(SNR,BER,'.-')
        case 'linear'
            plot(SNR,BER,'.-')
    end
    xlabel(xl)
    ylabel(yl)

end