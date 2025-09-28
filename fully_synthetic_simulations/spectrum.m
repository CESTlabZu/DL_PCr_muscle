function y = spectrum(amp, width, off)
    offset= -1000:25:1000;
    rffreq=[-2000, -1750, -1500, -1250, offset, 1250, 1500,1750,2000];
    y = (amp)./(1+(((rffreq-off).^2)./((0.5*width).^2)));
end
