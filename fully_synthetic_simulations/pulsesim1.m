function  p = pulsesim1(dw, ksw1,ksw2,ksw3, ksw4,ksw5, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2, R2S3,R2S4,  R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle)
    w1=getsatpulse(curve, angle, 1, duration, TR);
    init = pulsesolv1(w1, dw, ksw1, ksw2,ksw3,ksw4,ksw5, kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2, R2S3,R2S4,  R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, init, duration/1, excitewait);
    p = init;
end