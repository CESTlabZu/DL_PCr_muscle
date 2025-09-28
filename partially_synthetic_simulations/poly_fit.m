function background = plof_4(k, input,fitParam)
    input_background = input(5:37);
    k_background = k(5:37);
    input_background(18:28) = [];
    k_background(18:28) = [];
    FitResult.xindex = linspace(min(input),max(input),length(input_background));
    x0_background = [0.1,1,0.5,-190];
    lb=[0,0,0,-1000];
    ub=[100,100,1000,0];
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@background_func,x0_background,k_background',input_background,lb,ub,options,fitParam);
    background = background_func(FitResult.Coefficents,k(5:37),fitParam);

    
end