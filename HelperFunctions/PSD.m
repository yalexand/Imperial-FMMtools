function [ Omega, PSD ] = PSD( x, Fs )
        datasize=size(x);
        numsample=datasize(1);
        %Windowing
        %H=hann(numsample);
        %W=H.*x;
        W = x;
        %Fourier Transform
        FFTX=fft(x,numsample);
        %Power: magnitude^2
        X=FFTX(1:floor(numsample/2)).*conj(FFTX(1:floor(numsample/2)));
        %Bandwidth
        BW=1.5*Fs*numsample;
        %PSD=magnitude^2/bandwidth
        PSD=X/BW;
        %Computing the corresponding frequency values
        Omega=Fs*(0:numsample-1)/numsample;
        Omega=Omega(1:floor(numsample/2));
end

