function [t_tot,s_tot,d_tot,x_tot,y_tot,relevantTimes,relevantArcLengths,relevantXi,relevantYi,relevantDistances] = pp_computeArcLength2(points,xi,yi,velocities,deadTime)

currentX = 0;
currentY = 0;
pastIndex = 1;
s_tot = [];
t_tot = [];
d_tot = [];
x_tot = [];
y_tot = [];
totalDeltaT = 0;
totalDeltaS = 0;
totalDeltaD = 0;
relevantTimes = [];
relevantArcLengths = [];
relevantDistances = [];
relevantXi = [];
relevantYi = [];

closestPointsPath = pp_findClosestPointsPath(points,xi,yi);

for j=2:size(closestPointsPath,1)
    
    firstPointIndex = find(xi==closestPointsPath(j-1,1));
    secondPointIndex = find(xi==closestPointsPath(j,1));

    % Compute DeltaS in terms of steps of the current point from the
    % previous one
    DeltaS = abs(secondPointIndex-firstPointIndex);
    
    % The first sub-trajectory counts 1 sample less
    % because the start point is not counted
    % so we add 1 manually
    if j==2
        DeltaS = DeltaS+1;
    else
        DeltaS = DeltaS;
    end


    DeltaD = norm([xi(firstPointIndex),yi(firstPointIndex)]-[xi(secondPointIndex),yi(secondPointIndex)]);


    % v1d --> velocity in m/s
    % v2d --> velocity in samples/s
    v1d = velocities(j-1);
    v2d = velocities(j);
    
    % Convert the velocity from m/s to samples/s
    v1s = v1d*DeltaS/DeltaD;
    v2s = v2d*DeltaS/DeltaD;
    
    % Compute DeltaT with the given paper formula
    DeltaT = 2*DeltaS/(v1s+v2s);
    aT = (v2s-v1s)/DeltaT;
    % 
    % if v2d==0
    %     DeltaT = DeltaT+deadTime;
    % end
    
    % Compute the time vector needed to cover the whole DeltaS
    if j==2
        t = linspace(0,DeltaT,DeltaS);
    else
        % This is to avoid that the sub-trajectories between two subsequent
        % points overlap. So we jump the first sample
        t = linspace((DeltaT/DeltaS),DeltaT,DeltaS);
    end

    % Compute the arc length
    s = round((1/2)*aT*t.^2+v1s*t);
    
    % Compute the distance covered
    % Equal to the arc length multiplied by the distance between each
    % sample
    d = (DeltaD/DeltaS).*s;
    
    % Offsets
    s = s+totalDeltaS;
    d = d+totalDeltaD;

    % Compute the total distance covered in terms of steps and the total
    % time elapsed
    totalDeltaS = totalDeltaS+DeltaS;
    t = totalDeltaT+t;
    totalDeltaT = totalDeltaT+DeltaT;
    totalDeltaD = totalDeltaD+DeltaD;
    
    % Compute x
    x=[];
    for k=1:size(s,2)
        if s(k)==0
            s(k)=1;
        end
        x = [x, xi(s(k))];
    end

    % Compute y
    y=[];
    for k=1:size(s,2)
        if s(k)==0
            s(k)=1;
        end
        y = [y, yi(s(k))];
    end

    t_tot = [t_tot t];
    s_tot = [s_tot, s];
    d_tot = [d_tot, d];
    x_tot = [x_tot, x];
    y_tot = [y_tot, y];

    % Retrieve the points of the path at their corresponding time
    relevantTimes = [relevantTimes, totalDeltaT];
    relevantArcLengths = [relevantArcLengths, totalDeltaS];
    relevantDistances = [relevantDistances, totalDeltaD];
    relevantXi = [relevantXi, xi(totalDeltaS)];
    relevantYi = [relevantYi, yi(totalDeltaS)];

end

end

