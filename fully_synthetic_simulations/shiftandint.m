function outputspec = shiftandint(inputspec)    
    max=1000;
    step=25;
    offset= -max:step:max;
    k=[-2000, -1750, -1500, -1250, offset, 1250, 1500,1750,2000];
    k = k(5:85);
    input = inputspec(5:85);
    outputspec = inputspec;
    k_index = k(input == min(input));
    if(k_index ~= k(90/2))
        k_shift = k-(k_index);
        spec = interp1(k_shift,input,k);
    else
        spec = input';
    end
    outputspec(5:85) =  spec;
end

