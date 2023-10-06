function f=benchmark_func(x,func_num)
    bias = [100,1100,700,1900,1700,1600,2100,2200,2400,2500];
    f = cec20_func(x',func_num) - bias(func_num);
    f = f';
end