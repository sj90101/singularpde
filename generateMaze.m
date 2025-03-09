function maze = generateMaze(M, N)
    % Initialize maze with walls
    % maze(y, x, direction): logical array indicating if wall exists (true) or not (false)
    maze = true(M, N, 4); % Directions: N=1, E=2, S=3, W=4

    % Initialize visited array
    visited = false(M, N);

    % Start carving from (1,1)
    [maze, ~] = carve_passages_from(1, 1, maze, visited);
%    drawMaze(maze);
end

function [maze, visited] = carve_passages_from(cx, cy, maze, visited)
    visited(cy, cx) = true; % Mark current cell as visited

    % Directions
    dx = [0, 1, 0, -1];
    dy = [-1, 0, 1, 0];
    opposite = [3, 4, 1, 2];

    directions = randperm(4); % Random order of directions

    for i = 1:4
        direction = directions(i);
        nx = cx + dx(direction);
        ny = cy + dy(direction);

        % Check if nx and ny are within bounds and not yet visited
        if ny >= 1 && ny <= size(maze, 1) && nx >= 1 && nx <= size(maze, 2) && ~visited(ny, nx)
            % Remove walls between current cell and next cell
            maze(cy, cx, direction) = false;
            maze(ny, nx, opposite(direction)) = false;

            % Recursive call
            [maze, visited] = carve_passages_from(nx, ny, maze, visited);
        end
    end
end

function drawMaze(maze)
    [M, N, ~] = size(maze);

    % Initialize figure
    figure; hold all; axis equal off;

    for y = 1:M
        for x = 1:N
            % Coordinates
            x1 = x - 1;
            x2 = x;
            y1 = M - y;
            y2 = M - y + 1;

            if maze(y, x, 1) % Top wall
                plot([x1, x2], [y2, y2],'LineWidth',2.0);
            end
            if maze(y, x, 2) % Right wall
                plot([x2, x2], [y1, y2],'LineWidth',2.0);
            end
            if maze(y, x, 3) % Bottom wall
                plot([x1, x2], [y1, y1],'LineWidth',2.0);
            end
            if maze(y, x, 4) % Left wall
                plot([x1, x1], [y1, y2],'LineWidth',2.0);
            end
        end
    end

    % Adjust axis limits
    xlim([0, N]);
    ylim([0, M]);
    hold off;
end

