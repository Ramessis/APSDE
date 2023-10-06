function Run()
    clc
    close all
    format long
    format compact
    addpath('Public');
    %% Select CEC Version or UAV
    % 1: cec2017; 2:cec2020; 3: UAV
    run_version = 1;
    switch run_version
        case 1
            addpath('CEC2017');
            Fun_Nums = [1,3:29];
            runmax=51;
        case 2
            addpath('CEC2020');
            Fun_Nums = 1:10;
            runmax=51;
        case 3
            addpath('UAV');
            Fun_Nums = 1;
            runmax = 1;
    end
        
    % population size
    N=100;   

    TestFitness=[];
    TestResult=[];
    TestValue={};
    TestTime=[];
    TestRatio=[];
    TestFES=[];
    TestOptimization={};
    TestParameter={};
    
    %% APSDE
    name = 'APSDE';
    TestFitness=[];
    TestResult=[];
    TestValue={};
    TestTime=[];
    TestRatio=[];
    TestFES=[];
    TestOptimization={};
    TestParameter={};
    for problem=Fun_Nums
        RunResult = [];
        RunValue = [];
        RunTime = [];
        RunFES = [];
        RunOptimization = [];
        RunParameter = {};
        for run = 1:runmax
            [RunResult(run),RunValue(run,:),RunTime(run),RunFES(run),RunOptimization(run,:),RunParameter{1,run}]...
                =APSDE(problem,N,run);
        end
    
        TEV=Error(problem);
        sign=(RunResult<=TEV);
        Ratio=sum(sign)/runmax;
        FES=sign.*RunFES;
        TestFitness=[TestFitness;RunResult];
        TestResult=[TestResult;min(RunResult) max(RunResult) median(RunResult) mean(RunResult) std(RunResult)];
        TestValue=[TestValue;mean(RunValue)];
        TestTime=[TestTime;mean(RunTime)];
        TestRatio=[TestRatio;Ratio];
        TestFES=[TestFES;mean(FES)];
        TestOptimization=[TestOptimization;RunOptimization];
        TestParameter=[TestParameter;RunParameter];
    end
    Test=sprintf('Results/TestFitness/%s.mat',name);
    save(Test,'TestFitness');
    Test=sprintf('Results/TestResult/%s.mat',name);
    save(Test,'TestResult');
    Test=sprintf('Results/TestValue_FES/%s.mat',name);
    save(Test,'TestValue');
    Test=sprintf('Results/TestTime/%s.mat',name);
    save(Test,'TestTime');
    Test=sprintf('Results/TestRatio/%s.mat',name);
    save(Test,'TestRatio');
    Test=sprintf('Results/TestFES/%s.mat',name);
    save(Test,'TestFES');
    Test=sprintf('Results/TestOptimization/%s.mat',name);
    save(Test,'TestOptimization');
    Test=sprintf('Results/TestParameter_FES/%s.mat',name);
    save(Test,'TestParameter');
end

