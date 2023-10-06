function vi=boundConstraint(vi, lu);

    [NP, D] = size(vi);  % the population size and the problem's dimension

    xl = repmat(lu(1, :), NP, 1);
    xu = repmat(lu(2, :), NP, 1);
    %% check the lower bound
    pos = vi < xl;
    vi(pos) = 2 .* xl(pos) - vi(pos);
    
    %% check the upper bound
    pos = vi > xu;
    vi(pos) = 2 .* xu(pos) - vi(pos);
    
    vi = min(xu,vi);
    vi = max(xl,vi);
end