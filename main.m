clear all
close all
clc
load("data.mat")

% Plot formation and error
plot_formation(z)

% Position vector in 2D in time per drone/agent
z_pos = zeros(K,N,2);
z_pos(1,:,:) = z;
z = reshape(z_pos(1,:,:), size(z));

% Initialization
U = zeros(K,N,2);
dist = zeros(K,N,2);

% Positional error from the optimum location
pos_err = zeros(K,1);

% Array with the trajectory of one node in 2D
trajectory = zeros(K,N,2);

%% Formation control noiseless case
for k = 1:K
    for i = 4:N
        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Caluclate the current input
        U(k,i,:) = L(i,:)*(z_i-z);

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + 10*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
    trajectory(k,:,:) = z_pos(k,:,:);
end
plot_formation(z);

figure
plot(pos_err)
disp("Final error noiseless case")
disp(pos_err(end))

%% Formation control Noise case
% Initialization
for k = 1:K
    for i = 4:N
        % Generate noise
        v = randn(size(z))*R;

        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Caluclate the current input
        U(k,i,:) = L(i,:)*(z_i-z+v);

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + 10*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
    trajectory(k,:,:) = z_pos(k,:,:);
end
plot_formation(z);

figure
plot(pos_err)
disp("Final error noise case")
disp(pos_err(end))


%% Formation control Noise case with averaging over 10 samples estimator
% Initialization
for k = 1:K
    for i = 4:N
        % Generate noise
        v = randn(size(z))*R;

        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Calculate the current distance
        dist(k,:,:) = z_i-z+v;

        if k < 10
            distance = sum(dist(k:k-k,:,:),1) / 10;
        else
            distance = sum(dist(k:k-10,:,:),1) / 10;
        end
        distance = reshape(distance, 7,2);
        
        % Caluclate the current input
        U(k,i,:) = L(i,:)*distance;

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + 10*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
    trajectory(k,:,:) = z_pos(k,:,:);
end
plot_formation(z);

figure
plot(pos_err)
disp("Final error noise case")
disp(pos_err(end))