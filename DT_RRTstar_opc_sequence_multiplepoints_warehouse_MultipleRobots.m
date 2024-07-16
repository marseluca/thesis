close all
clear all
clear figure

% load('Workspaces\couples.mat');

mappp.points=[];

mappp.lines=[];

mapp = import_data_catia_vrml('semaphore.wrl');
figure;

for i=1:size(mapp.points,2)
    if mapp.points(3,i)==0.02
        mappp.points=[mappp.points; mapp.points(1,i) mapp.points(2,i)]; % to put in vertical
    end
end

map.points=unique(mappp.points,'rows'); % to select the unique points
map.points=map.points*2000;

for i=1:size(mapp.lines,2)
    if mapp.points(3,mapp.lines(1,i))==0.02 && mapp.points(3,mapp.lines(1,i))==0.02
        x1=mapp.points(1,mapp.lines(1,i)); 
        y1=mapp.points(2,mapp.lines(1,i));
        x2=mapp.points(1,mapp.lines(2,i)); 
        y2=mapp.points(2,mapp.lines(2,i));
        for j=1:size(map.points,1)
            if x1*2000==map.points(j,1) && y1*2000==map.points(j,2) 
                jj=j;
            elseif x1*2000==map.points(j,2) && y1*2000==map.points(j,1)
                jj=j;
            end
        end
        for z=1:size(map.points,1)
            if x2*2000==map.points(z,1) && y2*2000==map.points(z,2) 
                zz=z;
            elseif x2*2000==map.points(z,2) && y2*2000==map.points(z,1)
                zz=z;
            end
        end
        mappp.lines=[mappp.lines; jj zz];

    end
end
% map.lines=unique(mappp.lines, 'rows');
map.lines=mappp.lines;

% Plot the map lines
for i=1:size(map.lines,1)
    x1=map.points(map.lines(i,1),1);
    y1=map.points(map.lines(i,1),2);
    x2=map.points(map.lines(i,2),1);
    y2=map.points(map.lines(i,2),2);

    x=[x1 x2];
    y=[y1 y2];
    hold on;
    plot(x, y, '-r','LineWidth',2);
end

xlabel('[m]');
ylabel('[m]');


max_x=max(map.points(:,1));
min_x=min(map.points(:,1));

max_y=max(map.points(:,2));
min_y=min(map.points(:,2));

%% Triangulation
tic;
triangulation = delaunayTriangulation(map.points);
toc1=toc;

%% Uncomment the following line to show triangulation
% triplot(triangulation);
% hold on;
% plot(map.points(:,1), map.points(:,2), 'ro');

mmap.points=unique(map.points,'rows');

% for i=1:size(triangulation.Points,1)
%     text(triangulation.Points(i,1),triangulation.Points(i,2),num2str(i),'Color','red');
% end


%% Writing the number of the triangle in each triangle
% P=triangulation.Points;
% T=triangulation.ConnectivityList;
% 
% vertices1 = P(T(:,1), :);
% vertices2 = P(T(:,2), :);
% vertices3 = P(T(:,3), :);
% 
% triangle_centers = (vertices1 + vertices2 + vertices3) / 3;
% 
% num_triangles = size(T, 1);
% 
% for i = 1:num_triangles
%     % Get the center of the current triangle
%     center = triangle_centers(i, :);
% 
%     % Display the triangle index at the center
%     text(center(1), center(2), num2str(i), 'HorizontalAlignment', 'center');
% end



conc_edges_=[];
    for i=1:size(triangulation.ConnectivityList,1)
        conc_edges_=[conc_edges_; triangulation.ConnectivityList(i,1) triangulation.ConnectivityList(i,2); triangulation.ConnectivityList(i,2) triangulation.ConnectivityList(i,3)];
    end
    conc_edges=unique(conc_edges_,'rows');



% Nrobots = 1;
% Ngoals = 1;
OutputMatrix = [];

% for i =1:Nrobots
%     [x,y] = ginput(1);
%     robot_x_coords(i) = x;
%     robot_y_coords(i) = y;
%     plot(robot_x_coords(i), robot_y_coords(i), 'o', Color='r');
% end
% 
% for i =1:Ngoals
%     [x,y] = ginput(1);
%     goal_x_coords(i) = x;
%     goal_y_coords(i) = y;
%     plot(goal_x_coords(i), goal_y_coords(i), 'o', Color='r');
% end

maxIter=6000;
    
num_of_sequences=1;
num_paths=1;

% paths_registration = cell((num_paths*num_of_sequences*size(couples,1)), 1);
paths_registration = cell((num_paths*num_of_sequences), 1);

% paths_dtrrt = inf(num_paths*num_of_sequences*size(couples,1), 4);
paths_dtrrt = inf(num_paths*num_of_sequences, 4);
paths_dtrrt(:,1)=1;


maxLength = (max(cellfun(@(x) size(x, 1), paths_registration)));
% dynMatrix = cell(num_paths*num_of_sequences*size(couples,1), maxLength+2);
dynMatrix = cell(num_paths*num_of_sequences, maxLength+2);

toc2=[];

%%
for loopp = 1: 1      %size(couples,1)
    
    %% To get it in input
    % Define start and goal points
    % startNode = [robot_x_coords(loop), robot_y_coords(loop)];
    % goalNode = [goal_x_coords(1), goal_y_coords(1)];
    
    %% For warehouse 1
    startNode = [-200, -100];
    goalNode = [200, 0]; 

    %% For the couples of start and goal nodes
    % startNode=couples(loopp,1:2);
    % goalNode=couples(loopp,3:4);
    min_dist_triang=[];
    for i=1:size(triangulation.ConnectivityList,1)
        p1=triangulation.Points(triangulation.ConnectivityList(i,1),:);
        p2=triangulation.Points(triangulation.ConnectivityList(i,2),:);
        p3=triangulation.Points(triangulation.ConnectivityList(i,3),:);
    
        dist1=sqrt((p1(1)-goalNode(1))^2+(p1(2)-goalNode(2))^2);
        dist2=sqrt((p2(1)-goalNode(1))^2+(p2(2)-goalNode(2))^2);
        dist3=sqrt((p3(1)-goalNode(1))^2+(p3(2)-goalNode(2))^2);
    
        [minimum, index_minimum]=min([dist1;dist2;dist3]);
    
        min_dist_triang=[min_dist_triang; i minimum];
    end

    
    
    %% Iterate for the number of triangle sequences you want to obtain (in this case 1)
    for f=1:num_of_sequences
        tic;
        %% Main loop
        
        nRowStart=whichTriangle(startNode, triangulation); %first input is a coordinate (x and y)
        nRowGoal=whichTriangle(goalNode, triangulation);

        if nRowStart==nRowGoal
            tic;
            total_distance=sqrt((goalNode(1) - startNode(1))^2 + (goalNode(2) - startNode(2))^2);
            toc3=toc;
            for b=1:num_paths
                paths_dtrrt((((loopp-1)*num_paths)+b),2)=loopp;
                paths_dtrrt((((loopp-1)*num_paths)+b),3)=toc3;
                paths_dtrrt((((loopp-1)*num_paths)+b),4)=total_distance;
            end
            break;
        end
        
        queue = nRowStart;
        visitedTriangles = false(size(triangulation.ConnectivityList, 1), 1);
        visitedTriangles(nRowStart) = true;
        visited_order=[nRowStart];
        cont=0;
        
        first_triangle_edges=[triangulation.ConnectivityList(nRowStart,1) triangulation.ConnectivityList(nRowStart,2); triangulation.ConnectivityList(nRowStart,2) triangulation.ConnectivityList(nRowStart,3); triangulation.ConnectivityList(nRowStart,1) triangulation.ConnectivityList(nRowStart,3)];
        
        first_triangle_adj_triangles=[];
        for i=1:size(first_triangle_edges,1)
            m=findTriangleContainingEdge(triangulation.ConnectivityList, first_triangle_edges(i,:), visitedTriangles);
            first_triangle_adj_triangles=[first_triangle_adj_triangles; m];
        end
        
        realEdgesStart=findRealEdges(nRowStart, triangulation, map);
        virtualEdgesStart=findVirtualEdges(nRowStart, triangulation, realEdgesStart);
        
        loop=false;
        
        while ~isempty(queue)
        
                number_of_n=0;
                flag=false;
                currentTriangleIndex = queue(1);
                queue(1) = [];
        
                if currentTriangleIndex == nRowGoal
                    disp('Goal reached');
                    reconstructPath(triangulation, visited_order);
                    break;  
                end
        
                realEdges=findRealEdges(currentTriangleIndex, triangulation, map);
                virtualEdges=findVirtualEdges(currentTriangleIndex, triangulation, realEdges);
        
                adjacentTriangles=[];
        
                bueno=false;
                chosen=[];
                all_vertices=0;
                
                n_triangle=[];
                previous_vertices=triangulation.ConnectivityList(visited_order(end),:);
    
                for i=1:size(virtualEdges,1)
                        n_triangle=[n_triangle; findTriangleContainingEdge(triangulation.ConnectivityList, virtualEdges(i,:), visitedTriangles)];
                end
                remove=find(n_triangle==0);
                n_triangle(remove)=[];
    
                n_triangle_=n_triangle; %in order to keep the variable. Useful to find random_number
                
                n_row_triangle=[];
                for i=1:size(n_triangle_,1)
                    n_row_triangle=[n_row_triangle; find(min_dist_triang(:,1)==n_triangle_(i))];
                end
                min_dist_this_triang=min_dist_triang(n_row_triangle,:);
    
    
                while bueno==false %until I get a valid next triangle
                    all_vertices=all_vertices+1; %keeps the count of how many adjacencies I have already checked
                    
                    [min_next_tri, min_next_tri_index]=min(min_dist_this_triang(:,2));
                    if size(min_next_tri,1)>1 && size(min_next_tri_index,1)>1
                        min_next_tri=min_next_tri(1);
                        min_next_tri_index=min_next_tri_index(1);
                    end
                    n_row=find(min_dist_this_triang(:,2)==min_next_tri);
                    if size(n_row,1)>1
                        n_row=n_row(1);
                    end
                    n=min_dist_this_triang(n_row,1);
                    index_remove=find(n_triangle==n);
                    n_triangle(index_remove)=[];
                    
                    if n~=0
                        size_n=size(findRealEdges(n,triangulation, map),1);
                    else 
                        size_n=0;
                    end
                
                    if size(n,1)>0 && size(findRealEdges(n,triangulation, map),1)<=1
                        triang_n_edges=[triangulation.ConnectivityList(n,1) triangulation.ConnectivityList(n,2); triangulation.ConnectivityList(n,2) triangulation.ConnectivityList(n,3); triangulation.ConnectivityList(n,1) triangulation.ConnectivityList(n,3)];
                        for i=1:size(map.lines,1) %intersection check
                            for j=1:size(triang_n_edges,1)
                                if segmentsIntersect(map.points(map.lines(i,1),:), map.points(map.lines(i,2),:), triangulation.Points(triang_n_edges(j,1),:), triangulation.Points(triang_n_edges(j,2),:))
                                    bueno=false;
                                    flag=true;
                                    break;
                                end
                            end
                            if flag==true
                                break;
                            end
                        end %end of intersection check
    
                        if any(ismember(previous_vertices, virtualEdges(index_remove,1))) && any(ismember(previous_vertices, virtualEdges(index_remove,2)))
                            bueno=false;
                        else
                            bueno=true;
                            adjacentTriangle=n;
                            if adjacentTriangle>0
                                adjacentTriangles = [adjacentTriangles; adjacentTriangle];
                                visited_order=[visited_order; adjacentTriangle];
                            end
                        end
                    elseif size(n,1)>0 && size_n>1
                         number_of_n=number_of_n+1;
                    else 
                        number_of_n=number_of_n+1;
                    end
    
                    if number_of_n>1
                        % disp('Entering exitLoop');
                        % disp(visited_order);
                        [currentTriangleIndex, visited_order, adjacentTriangles] =  exitLoop(visited_order, triangulation, map, visitedTriangles, adjacentTriangles);
                        % disp('Exiting exitLoop');
                        % disp(visited_order);
                        bueno=true;
                    end
                end
               
                queue = [queue; adjacentTriangles];
                visitedTriangles(adjacentTriangles) = true;
        
                cont=cont+1;
        end

        
        %% Initialization for equidistant samples

        num_samples = 30; 
        samples = struct('rows', [], 'points', [], 'points_', []);
        
        %% Initialization for Gaussian distribution
    
        % num_samples = 30; 
        % samples = struct('rows', [], 'points', [], 'points_', []);
        % mean_val = 0.5;    
        % std_dev = 0.5;     
        % 
        %% Trolley parameters
    
        circle_radius = 0.003*2000; 
        safety_margin = circle_radius + circle_radius*0.05; %+il 5%
    
    
    
        %% Placing the equidistant samples on the edges
        for i = 1:numel(visited_order) - 1

            triangle1_idx = visited_order(i);
            triangle2_idx = visited_order(i + 1);

            % Get the vertices (points) of each triangle
            vertices1 = triangulation.Points(triangulation.ConnectivityList(triangle1_idx, :), :);
            vertices2 = triangulation.Points(triangulation.ConnectivityList(triangle2_idx, :), :);

            % Identify common vertices between the two triangles
            common_vertices = intersect(triangulation.ConnectivityList(triangle1_idx, :), triangulation.ConnectivityList(triangle2_idx, :));

            % Check if there are exactly two common vertices (one shared edge)
            if numel(common_vertices) ~= 2
                error('Unexpected number of common vertices between triangles.');
            end

            % Determine the indices of the common vertices within each triangle
            common_vertex_indices1 = find(ismember(triangulation.ConnectivityList(triangle1_idx, :), common_vertices));
            common_vertex_indices2 = find(ismember(triangulation.ConnectivityList(triangle2_idx, :), common_vertices));

            % Get the common vertices from each triangle
            common_vertex1 = vertices1(common_vertex_indices1, :);
            common_vertex2 = vertices2(common_vertex_indices2, :);

            % Calculate equidistant samples along the edge connecting common vertices
            t = linspace(0, 1, num_samples + 2)';  % Equally spaced parameter values
            t = t(2:end-1);  % Exclude the endpoints (0 and 1)

            % Check dimensions of common_vertex1 and common_vertex2
            if size(common_vertex1, 1) ~= size(common_vertex2, 1)
                error('Common vertices have incompatible dimensions.');
            end

            t = t(:);  % Force t to be a column vector

            % Calculate edge_samples using element-wise operations
            edge_samples = common_vertex1(1,:) + t.* (common_vertex1(2,:) - common_vertex1(1,:));

            % Get the other vertices (not part of the common edge) for each triangle
            other_vertices1 = vertices1(~ismember(1:size(vertices1, 1), common_vertex_indices1), :);
            other_vertices2 = vertices2(~ismember(1:size(vertices2, 1), common_vertex_indices2), :);

            % Append the data to the samples structure
            edge_indices = repmat(common_vertex_indices1, numel(t), 1);  % Repeated indices for the common edge
            edge_points = [edge_samples; other_vertices1; other_vertices2];

            samples.rows = [samples.rows; common_vertices];
            samples.points_ = [samples.points_; edge_points];
        end

        [~, idx_to_remove] = ismember(samples.points_, triangulation.Points, "rows");

        filtered_points = samples.points_(all(idx_to_remove == 0, 2), :);

        samples.points=filtered_points;

        %% Plot the equidistant points
        % for i = 1:size(samples.points, 1)
        %     % Edge's center
        %     x_center = samples.points(i, 1);
        %     y_center = samples.points(i, 2);
        % 
        %     % Draw circle on the edge
        %     theta = linspace(0, 2*pi, 100); 
        %     x_circle = x_center + circle_radius * cos(theta);
        %     y_circle = y_center + circle_radius * sin(theta);
        %     %plot(x_circle, y_circle, 'r');
        % end
        

        %% Placing the Gaussian distribution of samples on the edges
    
        % for i = 1:numel(visited_order) - 1
        % 
        %     if nRowStart==nRowGoal
        %         break;
        %     end
        % 
        %     triangle1_idx = visited_order(i);
        %     triangle2_idx = visited_order(i + 1);
        % 
        %     % Get the vertices (points) of each triangle
        %     vertices1 = triangulation.Points(triangulation.ConnectivityList(triangle1_idx, :), :);
        %     vertices2 = triangulation.Points(triangulation.ConnectivityList(triangle2_idx, :), :);
        % 
        %     % Identify common vertices between the two triangles
        %     common_vertices = intersect(triangulation.ConnectivityList(triangle1_idx, :), triangulation.ConnectivityList(triangle2_idx, :));
        % 
        %     % Check if there are exactly two common vertices (one shared edge)
        %     if numel(common_vertices) ~= 2
        %         error('Unexpected number of common vertices between triangles.');
        %     end
        % 
        %     % Determine the indices of the common vertices within each triangle
        %     common_vertex_indices1 = find(ismember(triangulation.ConnectivityList(triangle1_idx, :), common_vertices));
        %     common_vertex_indices2 = find(ismember(triangulation.ConnectivityList(triangle2_idx, :), common_vertices));
        % 
        %     % Get the common vertices from each triangle
        %     common_vertex1 = vertices1(common_vertex_indices1, :);
        %     common_vertex2 = vertices2(common_vertex_indices2, :);
        % 
        %     t = normrnd(mean_val, std_dev, num_samples, 1);
        % 
        %     % Clip the samples to the range [0, 1]
        %     t = min(max(t, 0), 1);
        % 
        %     % Check dimensions of common_vertex1 and common_vertex2
        %     if size(common_vertex1, 1) ~= size(common_vertex2, 1)
        %         error('Common vertices have incompatible dimensions.');
        %     end
        % 
        %     % Ensure t is a column vector matching the number of rows in common_vertex1
        %     t = t(:);  % Force t to be a column vector
        % 
        %     % Calculate edge_samples using element-wise operations
        %     edge_samples = common_vertex1(1,:) + t.* (common_vertex1(2,:) - common_vertex1(1,:));
        % 
        %     % Get the other vertices (not part of the common edge) for each triangle
        %     other_vertices1 = vertices1(~ismember(1:size(vertices1, 1), common_vertex_indices1), :);
        %     other_vertices2 = vertices2(~ismember(1:size(vertices2, 1), common_vertex_indices2), :);
        % 
        %     % Append the data to the samples structure
        %     edge_indices = repmat(common_vertex_indices1, numel(t), 1);  % Repeated indices for the common edge
        %     edge_points = [edge_samples; other_vertices1; other_vertices2];
        % 
        %     samples.rows = [samples.rows; edge_indices];
        %     samples.points_ = [samples.points_; edge_points];
        % end
        % 
        % if nRowStart~=nRowGoal
        %     [~, idx_to_remove] = ismember(samples.points_, triangulation.Points, "rows");
        % 
        %     filtered_points = samples.points_(all(idx_to_remove == 0, 2), :);
        % 
        %     samples.points=filtered_points;
        % end 

        %% Plot the Gaussian distribution points
        % hold on;
        % plot(samples.points_(:, 1), samples.points_(:, 2), 'ro'); 
        
        toc2_=toc;
        toc2=[toc2; toc2_];
    
        %% Generating paths
        for b=1:num_paths
    
            tic;
            total_distance=0;
            % disp('Path number');
            % disp(b);
        
            hold on;
            points_sequence = [startNode];
            current_point = startNode;
        
        
            stop=false;
            for i = 1:num_samples:size(samples.points, 1)-num_samples+1
                % disp('i:');
                % disp(i);
                % disp('---------------------');
                random_index_taken=[];
                contin=true;
                another_random=false;
        
                while contin
                    if size(random_index_taken,1)>=num_samples
                        disp('Failed to find a collision-free connection :(');
                        stop=true;
                        % contin=false;
                        break;
                    end
                    random_index = randi([i, i + num_samples - 1]);
                    % disp('Random index:');
                    % disp(random_index);
                    if ~ismember(random_index, random_index_taken)
                        random_index_taken=[random_index_taken; random_index];
                        new_point = samples.points(random_index, :);
                        another_random=false;
                    else 
                        disp('Random index already taken, trying with another one');
                        another_random=true;
                    end
        
                    % intersects_map_line=collisionCheckAlong(current_point, new_point, safety_margin, map);
        
                    if another_random==false
                        intersects_map_line = checkCircleSegmentIntersection(current_point, new_point, safety_margin, map);
                        % disp('Intersection check');
                        % disp(intersects_map_line);
        
                        circle_ok=false;
                        for j = 1:size(map.lines, 1)
                            segment_start = map.points(map.lines(j, 1),:);
                            segment_end = map.points(map.lines(j, 2),:);
        
                            % Check intersection between circle and line segment
                            if circleLineIntersect(safety_margin, new_point, circle_radius, segment_start, segment_end)
                                circle_ok = false;
                                break; 
                            end
                         end
        
                        if ~intersects_map_line
        
                            contin=false;
                            total_distance=total_distance + sqrt((new_point(1) - current_point(1))^2 + (new_point(2) - current_point(2))^2);
                            current_point = new_point;
                            points_sequence = [points_sequence; new_point];
        
                        else
                            contin=true;
                        end
                       
                    end
                end
                if stop==true
                    % disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                    break;
                end
        
            end
    
            total_distance=total_distance+sqrt((goalNode(1) - points_sequence(end,1))^2 + (goalNode(2) - points_sequence(end,2))^2);
            % if f==1
            %     if b==1
            %         total_distance11=total_distance;
            %     elseif b==2
            %         total_distance21=total_distance;
            %     elseif b==3
            %         total_distance31=total_distance;
            %     end
            % elseif f==2
            %     if b==1
            %         total_distance12=total_distance;
            %     elseif b==2
            %         total_distance22=total_distance;
            %     elseif b==3
            %         total_distance32=total_distance;
            %     end
            % elseif f==3
            %     if b==1
            %         total_distance13=total_distance;
            %     elseif b==2
            %         total_distance23=total_distance;
            %     elseif b==3
            %         total_distance33=total_distance;
            %     end
            % end
        
            % if isOnSegment(triangulation.ConnectivityList(nRowGoal,1), triangulation.ConnectivityList(nRowGoal,2) ,points_sequence(end)) || isOnSegment(triangulation.ConnectivityList(nRowGoal,2), triangulation.ConnectivityList(nRowGoal,3) ,points_sequence(end)) || isOnSegment(triangulation.ConnectivityList(nRowGoal,1), triangulation.ConnectivityList(nRowGoal,3) ,points_sequence(end))
            points_sequence=[points_sequence; goalNode];
            % end
        
            x = points_sequence(:, 1); % Extract all rows from the first column
            y = points_sequence(:, 2); % Extract all rows from the second column
        
            hold on;
            %% To plot the paths in different colours
            % if f==1
            %     if b==1 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1); % 'bo-' specifies blue circles connected by lines
            %     elseif b==2 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1);
            %     elseif b==3 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1);
            %     elseif b==4 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1);
            %     elseif b==5 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1);
            %     elseif b==6 && stop==false
            %         plot(x, y, 'bo-', 'LineWidth', 1);
            %     end
            % end
            % if f==2
            %     if b==1 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1); % 'bo-' specifies blue circles connected by lines
            %     elseif b==2 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1);
            %     elseif b==3 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1);
            %     elseif b==4 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1);
            %     elseif b==5 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1);
            %     elseif b==6 && stop==false
            %         plot(x, y, 'mo-', 'LineWidth', 1);
            %     end
            % end
            % if f==3
            %     if b==1 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1); % 'bo-' specifies blue circles connected by lines
            %     elseif b==2 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1);
            %     elseif b==3 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1);
            %     elseif b==4 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1);
            %     elseif b==5 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1);
            %     elseif b==6 && stop==false
            %         plot(x, y, 'go-', 'LineWidth', 1);
            %     end
            % end
            % disp('/////////////////////////////////////////////////////////////////////////////////////////////////');
    
    
            toc3=toc;
            %% Filling the "paths" variable with 1) n of map 2) sequence number 3) exec. time 4) path length
            % if loopp==1
            %     paths_dtrrt(loopp*b,2)=loopp;
            %     paths_dtrrt(loopp*b,3)=toc3;
            %     paths_dtrrt(loopp*b,4)=total_distance;
            % else

            paths_dtrrt((((loopp-1)*num_paths)+b),2)=loopp;
            paths_dtrrt((((loopp-1)*num_paths)+b),3)=toc3;
            paths_dtrrt((((loopp-1)*num_paths)+b),4)=total_distance;

            % end
            %% Registering paths
            if f==1
                paths_registration{b}=points_sequence;
            else
                paths_registration{((f-1)*num_paths)+b}=points_sequence;
            end

    
    
                %% Filling the dynMatrix
            %     if f==1
            %         numPoints = size(paths_registration{b}, 1);
            %         dynMatrix{b, 1} = f; 
            %         dynMatrix{b, 2} = total_distance;  
            % 
            %         for j = 1:numPoints
            %             dynMatrix{b, 2 + j} = paths_registration{b}(j, :); % [x, y]
            %         end
            %     else
            %         numPoints = size(paths_registration{((f-1)*num_paths)+b}, 1);
            %         dynMatrix{((f-1)*num_paths)+b, 1} = f; 
            %         dynMatrix{((f-1)*num_paths)+b, 2} = total_distance;
            % 
            %         for j = 1:numPoints
            %             dynMatrix{((f-1)*num_paths)+b, 2 + j} = paths_registration{((f-1)*num_paths)+b}(j, :); % [x, y]
            %         end
            %     end
            %     maxLength = max(cellfun(@(x) size(x, 1), paths_registration));
            %     for i = 1:((f-1)*num_paths)+b
            %         for j = size(paths_registration{i}, 1) + 3:maxLength+2
            %             dynMatrix{i, j} = [Inf, Inf];
            %         end
            %     end
            % 
            % 
        end

    end
    
% 
%     dynMatrix=dynMatrix(:,1:maxLength+2);
%     if loop ==1
%         OutputMatrix = dynMatrix;
%     else
%         if size(OutputMatrix,2) < size(dynMatrix,2)
%             AddColumns = size(dynMatrix,2) - size(OutputMatrix,2);
%             NewCol = cell(size(OutputMatrix,1),AddColumns);
%             NewColPop = Inf(1,2);
%             for n=1:AddColumns
%                 for m=1:size(OutputMatrix,1)
%                     NewCol{m,n} = NewColPop;
%                 end
%             end
%             %NewCol = Inf(size(OutputMatrix,1),AddColumns);
%             %NewCol = num2cell(NewCol); 
%             OutputMatrix = [OutputMatrix NewCol];
%             OutputMatrix = [OutputMatrix; dynMatrix];
%         elseif size(OutputMatrix,2) > size(dynMatrix,2)
%             AddColumns = size(OutputMatrix,2) - size(dynMatrix,2);
%             NewCol = cell(size(dynMatrix,1),AddColumns);
%             NewColPop = Inf(1,2);
%             for n=1:AddColumns
%                 for m=1:size(dynMatrix,1)
%                     NewCol{m,n} = NewColPop;
%                 end
%             end
%             %NewCol = Inf(size(dynMatrix,1),AddColumns);
%             %NewCol = num2cell(NewCol);
%             dynMatrix = [dynMatrix NewCol];
%             OutputMatrix = [OutputMatrix; dynMatrix];
% 
%         end
%     end
% 
% 
end
% 
%     updatedMatrix=cell(size(OutputMatrix,1),size(OutputMatrix,2));
%     % oldDyn=dynMatrix;
% 
% 
%     [numRows, ~] = size(OutputMatrix);
% 
%     % Iterate on each row of the dynMatrix
%     for i = 1:numRows
%         % Extract coordinates of points of the current row
%         rowCoords= OutputMatrix(i, 3:end);
% 
%         % Initialize an array to store the valid points
%         points = [];
%         k = 1;
%         cont_inf=0; %number of [Inf, Inf] in the row
% 
%         % Collect the couples of coordinates [x, y] until [Inf, Inf]
%         while k + 1 <= length(rowCoords)
%             if ~isequal(rowCoords(k), {[Inf, Inf]}) && ~isequal(rowCoords(k + 1), {[Inf, Inf]})
%                 points = [points; rowCoords{k}(1) rowCoords{k}(2); rowCoords{k + 1}(1) rowCoords{k + 1}(2)];
%             elseif ~isequal(rowCoords(k), {[Inf, Inf]}) && isequal(rowCoords(k + 1), {[Inf, Inf]})
%                 points = [points; rowCoords{k}(1) rowCoords{k}(2)];
%                 cont_inf=cont_inf+1;
%             elseif isequal(rowCoords(k), {[Inf, Inf]}) && ~isequal(rowCoords(k + 1), {[Inf, Inf]})
%                 points = [points; rowCoords{k + 1}(1) rowCoords{k + 1}(2)];
%                 cont_inf=cont_inf+1;
%             elseif isequal(rowCoords(k), {[Inf, Inf]}) && isequal(rowCoords(k + 1), {[Inf, Inf]})
%                 cont_inf=cont_inf+2;
%             end
%             k=k+2;
%         end
% 
%         if mod(length(rowCoords), 2) ~= 0
%             if isequal(rowCoords(end), {[Inf, Inf]})
%                 cont_inf=cont_inf+1;
%             else
%                 points=[points; rowCoords(end)];
%             end
%         end
% %% Replace all the [Inf, Inf] with middle points until there are no more [Inf, Inf]
%         newPoints_=points;
%         j=1;
%         if cont_inf>0
%             while(cont_inf~=0)
%                 point1 = newPoints_(j, :);
%                 point2 =  newPoints_(j + 1, :);
%                 midPoint = (point1 + point2) / 2;
%                 newPoints_ = [newPoints_(1:j, :); midPoint; newPoints_(j+1:end, :)];
%                 cont_inf=cont_inf-1;
%                 j = j + 2;
%                 if j>=size(newPoints_,1)
%                     j=1;
%                 end
%             end
%             newRowCoords = cell(1, size(newPoints_, 1)); %create a cell array from the newPoints matrix
%             for z = 1:size(newPoints_, 1)
%                 newRowCoords{z} = newPoints_(z, :);
%             end
%             updatedMatrix(i, 1)=OutputMatrix(i, 1);
%             updatedMatrix(i, 2)=OutputMatrix(i, 2);
%             updatedMatrix(i, 3:2+length(newRowCoords)) = newRowCoords;
%         else
%             updatedMatrix(i, :)=OutputMatrix(i, :);
%         end
% 
%     end
% 
% 
% excel_file = 'updatedMatrix.xlsx';
% writecell(updatedMatrix, excel_file);







