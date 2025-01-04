clear all
close all
clc
load("data.mat")
z_org = z;

% Amount of iterations
dt = 0.01;
K = 10000;
noise_power = 1;
control_law_speed = 10;

% Plot formation and error
plot_formation(z, "Initial state")
plot_formation(z_star, "Desired state")

% Position vector in 2D in time per drone/agent
z_pos = zeros(K,N,2);
z_pos(1,:,:) = z_org;
z = reshape(z_pos(1,:,:), size(z));

% Initialization
U = zeros(K,N,2);
dist = zeros(K,N,2);

% Positional error from the optimum location
pos_err = zeros(K,1);

% Formation control noiseless case
for k = 1:K
    for i = 4:N
        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Caluclate the current input
        U(k,i,:) = L(i,:)*(z_i-z);

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + control_law_speed*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
end
% plot_formation(z, "Final state, noiseless case");

z_pos(end,1:3,:) = z_org(1:3,:); %Fill the start positions of the first 3 nodes
plot_formation_trajectory(z_pos,"Final state, noiseless case");

figure
plot(pos_err)
yscale("log")
grid("on")
ylabel("Procrutes error")
xlabel("Time")
title("Error, noiseless case")

disp("Final error, noiseless case")
disp(pos_err(end))

%% Formation control Noise case

% Position vector in 2D in time per drone/agent
z_pos = zeros(K,N,2);
z_pos(1,:,:) = z_org;
z = reshape(z_pos(1,:,:), size(z));

% Initialization
U = zeros(K,N,2);
dist = zeros(K,N,2);

% Positional error from the optimum location
pos_err = zeros(K,1);

for k = 1:K
    for i = 4:N
        % Generate noise
        v = noise_power*randn(size(z))*R;

        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Caluclate the current input
        U(k,i,:) = L(i,:)*(z_i-z+v);

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + control_law_speed*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
end
% plot_formation(z, "Final state, noise case");

z_pos(end,1:3,:) = z_org(1:3,:); %Fill the start positions of the first 3 nodes
plot_formation_trajectory(z_pos,"Final state, noise case");

figure
plot(pos_err)
yscale("log")
grid("on")
ylabel("Procrutes error")
xlabel("Time")
title("Error, noise case")

disp("Final error, noise case")
disp(pos_err(end))

%% Formation control Noise case with averaging over 10 samples estimator

% Position vector in 2D in time per drone/agent
z_pos = zeros(K,N,2);
z_pos(1,:,:) = z_org;
z = reshape(z_pos(1,:,:), size(z));

% Initialization
U = zeros(K,N,2);
dist = zeros(K,N,2);

% Positional error from the optimum location
pos_err = zeros(K,1);
MA_size = 10;

for k = 1:K
    for i = 4:N
        % Generate noise
        v = noise_power*randn(size(z))*R;

        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));

        % Calculate the current distance
        dist(k,:,:) = z_i-z+v;

        if k <= MA_size
            distance = sum(dist(1:k,:,:),1) / k;
        else
            distance = sum(dist(k-MA_size:k,:,:),1) / MA_size;
        end
        distance_reshaped = reshape(distance,7,2);

        % Caluclate the current input
        U(k,i,:) = L(i,:)*distance_reshaped;

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + control_law_speed*dt*U(k,i,:);

        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);
end
% plot_formation(z, "Final state, noise case with MA estimator");

z_pos(end,1:3,:) = z_org(1:3,:); %Fill the start positions of the first 3 nodes
plot_formation_trajectory(z_pos, "Final state, noise case with MA estimator");

figure
plot(pos_err)
yscale("log")
grid("on")
ylabel("Procrutes error")
xlabel("Time")
title("Error, noise case with MA estimator")

disp("Final error, noise case with MA estimator")
disp(pos_err(end))

%% Formation control Noise case with MLE, Gaussian noise and linear model estimator

% Position vector in 2D in time per drone/agent
z_pos = zeros(K,N,2);
z_pos(1,:,:) = z_org;
z = reshape(z_pos(1,:,:), size(z));

% Initialization
U = zeros(K,N,2);
dist = zeros(K,N,2);

% Positional error from the optimum location
pos_err = zeros(K,1);

T = 100;
D = 2;
t_vec = linspace(1,1,T)';
zHat = zeros(N,D);

% Kronicker product of 1_T and I_D
H = kron(t_vec,eye(D));

% Kronicker product of 1_D and R
Rtilde = kron(eye(T),R);

for k = 1:K
    for i = 4:N
        % Reshape z_pos per node a 2D matrix
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));
        
        for j = 1:N %Can be faster by only using neighbours
            if L(i,j) ~= 0
                % Generate noise
                v = noise_power*Rtilde*randn(D*T,1);

                y = H*(z_i - z(j,:))' + v;
                zHat(j,:) = (1/(T))*H'*y;
            end
        end

        % Caluclate the current input
        U(k,i,:) = L(i,:)*zHat;

        % Change position according to input
        z_pos(k+1,i,:) = z_pos(k,i,:) + control_law_speed*dt*U(k,i,:);
        
        % Reshape 2D z_pos per node to fill into z with all nodes
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err(k) = norm(z-z_star,2);

end
% plot_formation(z, "Final state, noise case with MLE estimator");

z_pos(end,1:3,:) = z_org(1:3,:); %Fill the start positions of the first 3 nodes
plot_formation_trajectory(z_pos, "Final state, noise case with MLE estimator");

figure
plot(pos_err)
yscale("log")
title("Error, noise case with MLE estimator")

disp("Final error, MLE noise case")
disp(pos_err(end))