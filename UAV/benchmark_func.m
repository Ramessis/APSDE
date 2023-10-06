 %% Fitness function
 % x is the UAV flight path
 % f is the UAV flight cost (path length + capture value)
function f=benchmark_func(x,~)
    [NP,D] = size(x);
    interval = 100/(D-1);   % interval length
    angleLimitValue = 100;  % Penalty value for UAV turning angle exceeding the limiting angle
    capturedValue = 1000;    % Increased fitness value after being captured
    
    VR = 2; % Ratio of blue UAV speed to red one speed, range (1, infinity)
    warningPoints = [20,00,10;60,80,10;50,45,20; 100,20,20;100,80,20];% Coordinates and radius of the center of the circle
    
    f = zeros(NP,1);

    for i = 1:NP
        watchedPoints = [];
        watchedPoints = Watching([0*interval,x(i,1)],[0*interval,x(i,1)],warningPoints,watchedPoints);
        
        for j = 2:1:D
            % Detect whether the current point is the first entry into the alert area
            watchedPoints = Watching([(j-2)*interval,x(i,j-1)],[(j-1)*interval,x(i,j)],warningPoints,watchedPoints);
            
            step = sqrt(interval^2 + (x(i,j)-x(i,j-1))^2);   % Single-step flight distance            
            f(i,1) = f(i,1) + step + floor(step/(2*interval))*angleLimitValue;

            isCaptured = caputreFunc([(j-1)*interval,x(i,j)],watchedPoints, VR); % Was the UAV captured
            
            if isCaptured
                f(i,1) = f(i,1) + capturedValue;
            end
        end
    end
end

%% Detect whether the current point is the first entry into the alert area
function watchedPoints = Watching(priorPoint, point,warningPoints,watchedPoints)
    [px,~] = size(warningPoints);
    
    for i = 1:px
        % Distance of the current point from the alert point
        dist = sqrt((point(1)-warningPoints(i,1))^2 + (point(2)-warningPoints(i,2))^2);
        
        % Whether or not the UAV enters the cordoned-off point surveillance area
        if dist <= warningPoints(i,3) 
            [wx,~] = size(watchedPoints);
            isAdded = false;
            for j = 1:wx
                % Surveillance has been initiated at this perimeter.
                if warningPoints(i,1) == watchedPoints(j,1) && warningPoints(i,2) == watchedPoints(j,2) 
                    isAdded = true;
                end
            end
            
            if ~isAdded
                if dist < (warningPoints(i,3) - 1)
                    watchedPoints = [watchedPoints;warningPoints(i,1),warningPoints(i,2),priorPoint(1),priorPoint(2)];
                else
                    watchedPoints = [watchedPoints;warningPoints(i,1),warningPoints(i,2),point(1),point(2)];
                end
            end
        end
    end
end

%% Was the UAV captured
function isCaptured = caputreFunc(point,watchedPoints,VR)
    [px,~] = size(watchedPoints);
    isCaptured = false;
    for i = 1:px
        dist1 = sqrt((point(1) - watchedPoints(i,3))^2 + (point(2) - watchedPoints(i,4))^2);% Distance to Trigger Monitoring Points
        dist2 = sqrt((point(1) - watchedPoints(i,1))^2 + (point(2) - watchedPoints(i,2))^2);% Distance to Surveillance Points Drone Distance
        
        % Uncaptured
        if dist1/dist2 < VR 
            isCaptured = false;
        % Captured
        else
            isCaptured = true;
            return;
        end
    end
end

