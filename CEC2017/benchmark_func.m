function f=benchmark_func(x,func_num)
    bias = [100:100:2900];
    f = cec17_func(x',func_num) - bias(func_num);
    f = f';
end