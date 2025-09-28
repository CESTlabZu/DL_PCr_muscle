function [y,Rbak] = background_func(x,xdata,FitParam)
    peakoffset = FitParam.PeakOffset;
    Rbak = x(1) + x(2)*(xdata-peakoffset) + x(3)*((xdata-peakoffset).^2)+ x(4)*((xdata-peakoffset).^3);
    Rall = Rbak;
    y = Rall;
end