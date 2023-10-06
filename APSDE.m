%% Adaptive population size based differential evolution (APSDE)
%problem: the serial number of testing function recorded in "Public\benchmark_func.m"
%N: the population size
%runmax: the number of the algorithm runs
%RunResult: the  optimal value produced by each algorithm runs
%RunOptimization: the optimal value produced by reach algorithm runs
%RunValue: the fitness of optimal value produced by each 10000 FES
%RunParameter:the optimal value produced by each 10000 FES
%RunTime: the time spent by each algorithm runs
%RunFES: the FES required to satisfy the conditions
function RunResult=APSDE(problem,~,run)
    'APSDE'
    D=Dim(problem);%13-16行的意思参考CEP
    lu=Boundary(problem,D);
    TEV = Error(problem);
    FESMAX = D*10000;
    
    TimeFlag=0;
    TempFES=FESMAX;
    t1=clock;
    
    %%  parameter settings

    rate_popSize = 4;
    p_best_rate = 0.11;

    
    init_popSize = round(D * rate_popSize);
    min_popSize = 4;
    pop_size = init_popSize;
    

    %% Initialize the main population
    pop = repmat(lu(1, :), pop_size, 1) + rand(pop_size, D) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
    fitness = benchmark_func(pop, problem);
    best = min(fitness);
    previous_popSize = pop_size;
    
    FES = pop_size;
    stagnation = 0;
    
    mu_sf = zeros(pop_size,1);
    mu_cr = zeros(pop_size,1);
    PFIA_size = init_popSize; %H * D;   

    
    p = 0.3;            % similarity rate
    PFIA = cell(5,1);                  % PFIA=(features pair, MF, MCR, pop_size/init_popSize,SR(success rate));
    PFIA{1,1} = ones(PFIA_size,3) .* inf; 
    PFIA{2,1} = ones(PFIA_size,1) .* 0.5;
    PFIA{3,1} = ones(PFIA_size,1) .* 0.5;
    PFIA{4,1} = ones(PFIA_size,1) .* (init_popSize/init_popSize); 
    PFIA{5,1} = ones(PFIA_size,1) .* 0.000; 
    PFIA_pos = 1;

    archive.NP = init_popSize; % the maximum size of the archive
    archive.pop = repmat(lu(1, :), archive.NP, 1) + rand(archive.NP, D) .* (repmat(lu(2, :) - lu(1, :), archive.NP, 1));
    archive.funvalues = benchmark_func(archive.pop, problem);
    archive.flag = 1; % first in, first out.
    archive.index = 1; % the updated position 
    k=1;
    for i=1:archive.NP
        FES=FES+1;
        if FES==10000*0.1||mod(FES,10000)==0
            [kk,ll]=min(fitness);
            RunValue(1,k)=kk;
            Para(k,:)=pop(ll,:);
            k=k+1;
            fprintf('Algorithm:%s problemIndex:%d Run:%d FES:%d Best:%g\n','APSDE',problem,run,FES,kk);
        end
        if TimeFlag==0
            if min(fitness)<=TEV
                TempFES=FES;
                TimeFlag=1;
            end
        end
    end
    %% main loop
    while FES <= FESMAX
        %calculating current population feature
        currentPopStd = sum(std(pop))/pop_size;
        currentFitStd = std(fitness)/pop_size;
        currentState = [currentPopStd, currentFitStd, stagnation];
        % state matching 
        try
            similarity = pdist2(PFIA{1,1},currentState,'mahal'); %Calculating feature similarity by using Mahalanobis distance
        catch ME %#ok<NASGU>
            similarity = sqrt(sum((currentState - PFIA{1,1}).^2,2)); %Calculating feature similarity by using Euclidean distance
        end
        [~, index] = sort(similarity, 'ascend');
        
        %% calculating the size of current population
        popSize_index = index(1:ceil(p * PFIA_size));
        history_popSize = PFIA{4,1}(popSize_index,1);
        history_IR = PFIA{5,1}(popSize_index,1);
        mu_popSize = history_IR' * (history_popSize.^2)/(history_IR'* history_popSize);

        if isnan(mu_popSize)
            if(FES/FESMAX < 0.5)
                rate_popSize = 0.75 * mean(history_popSize); %(min_popSize/init_popSize - 1)/FESMAX * FES + 1; %FES/FESMAX;
            else
                rate_popSize = 1.25 * mean(history_popSize);
            end
        else
            rate_popSize =  normrnd(mu_popSize,0.1); %mu_popSize + 0.1 * tan(pi * (rand - 0.5));
        end

        if FES/FESMAX > 0.75 && rate_popSize > 0.25
            rate_popSize = normrnd(0.15, 0.1);  %20220403-normrnd(0.2, 0.1);
        end
        
        if rate_popSize > 1
            rate_popSize = 1;
        end

        if rate_popSize < (min_popSize/init_popSize)
            rate_popSize = min_popSize/init_popSize;
        end
        
        pop_size = ceil(rate_popSize * init_popSize);

        if pop_size < previous_popSize
            [~,sorted_index] = sort(fitness);
            pop = pop(sorted_index(1:pop_size),:);
            fitness = fitness(sorted_index(1:pop_size),:);
        end

        if pop_size > previous_popSize
           [~,sorted_index_arc] = sort(archive.funvalues);
           pop = [pop;archive.pop(sorted_index_arc(1:(pop_size - previous_popSize)),:)]; 
           fitness = [fitness;archive.funvalues(sorted_index_arc(1:(pop_size - previous_popSize)),:)];
        end
        previous_popSize = pop_size;

        %% each individual select randomly a pair (F, CR) form top p% similarity records
        top_index = ceil(max(1,p * PFIA_size) * rand(pop_size, 1));
        PFIA_rand_index = index(top_index,1);
        mu_sf = PFIA{2,1}(PFIA_rand_index,1);
        mu_cr = PFIA{3,1}(PFIA_rand_index,1);

        %% for generating crossover rate
        CR = normrnd(mu_cr, 0.1);
        CR = min(CR, 1);
        pos = find(CR <= 0);

        while ~ isempty(pos)
            CR(pos) = normrnd(mu_cr(pos), 0.1);
            pos = find(CR <= 0);
        end

        %% for generating scaling factor
        F = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
        pos = find(F <= 0);

        while ~ isempty(pos)
            F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
            pos = find(F <= 0);
        end

        F = min(F, 1);

        r0 = 1 : pop_size;
        popAll = [pop; archive.pop];
        [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);

        [~, sorted_index] = sort(fitness, 'ascend');
        pNP = ceil(p_best_rate * pop_size); 
        randindex = ceil(rand(pop_size, 1) .* pNP); %% select from [1, 2, 3, ..., pNP]
        randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
        pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

        vi = pop + F(:, ones(1, D)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
        vi = boundConstraint(vi,lu);

        mask = rand(pop_size, D) > CR(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
        rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
        jrand = sub2ind([pop_size D], rows, cols); mask(jrand) = false;
        ui = vi; ui(mask) = pop(mask);

        children_fitness = benchmark_func(ui, problem);

        dif = abs(fitness - children_fitness);


        %% I == 1: the parent is better; I == 2: the offspring is better
        I = (fitness > children_fitness);
        goodCR = CR(I == 1);
        goodF = F(I == 1);
        dif_val = dif(I == 1);

        archive = updateArchiveSequentially(archive, pop(I == 1, :), fitness(I == 1));

        [fitness, I] = min([fitness, children_fitness], [], 2);
        pop(I == 2, :) = ui(I == 2, :);
        fitness(I==2,:) = children_fitness(I==2,:);
        num_success_params = numel(goodCR);
        
        currentbest = min(fitness);
        if best ~= 0
            improve_rate = (best - currentbest)/best;
        else
            improve_rate = 1;
        end

        if currentbest < best
            best = currentbest;
            stagnation = 0;
        else
            stagnation = stagnation + 1;
        end

        if num_success_params > 0
            % UPdate the archive
            PFIA{1,1}(PFIA_pos,:) = currentState;%[currentPopStd, currentFitStd, stagnation];
        
            sum_dif = sum(dif_val);
            dif_val = dif_val / sum_dif;
            
            %% for updating the memory of scaling factor
            PFIA{2,1}(PFIA_pos,1) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

            %% for updating the memory of crossover rate
            PFIA{3,1}(PFIA_pos,1) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);

            %% for updating the pop_size and success rate
            PFIA{4,1}(PFIA_pos,1) = rate_popSize;
            PFIA{5,1}(PFIA_pos,1) = improve_rate;    %sum_dif/pop_size;             %num_success_params / pop_size;
        end

        PFIA_pos = PFIA_pos + 1;
        
        if PFIA_pos > PFIA_size 
            PFIA_pos = 1; 
        end

        for i=1:pop_size
            FES=FES+1;
            if FES==10000*0.1||mod(FES,10000)==0
                [kk,ll]=min(fitness);
                RunValue(1,k)=kk;
                Para(k,:)=pop(ll,:);
                k=k+1;
                fprintf('Algorithm:%s problemIndex:%d Run:%d FES:%d Best:%g\n','APSDE',problem,run,FES,kk);
            end
            if TimeFlag==0
                if min(fitness)<=TEV
                    TempFES=FES;
                    TimeFlag=1;
                end
            end
        end
    end
    [kk,ll] = min(fitness);
    gbest = pop(ll,:);
    t2 = clock;
    RunTime = etime(t2,t1);
    RunResult = kk;
    RunFES = TempFES;
    RunOptimization = gbest;
    RunParameter = Para;
end %% end 1 function run


function [r1, r2] = gnR1R2(NP1, NP2, r0)

    % gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
    %    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
    %    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
    %
    % Call:
    %    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
    %    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
    %
    % Version: 2.1  Date: 2008/07/01
    % Written by Jingqiao Zhang (jingqiao@gmail.com)

    NP0 = length(r0);

    r1 = floor(rand(1, NP0) * NP1) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    r2 = floor(rand(1, NP0) * NP2) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos)==0
            break;
        else % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end
end