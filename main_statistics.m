clear all
close all
clc
load("data.mat")
z_org = z;

% Parameters
dt = 0.01;
K = 50000;
noise_power = 1;
control_law_speed = 1;
maxTime = K * dt * control_law_speed; %in sec

% Start the statistics routine here
M = 5;
Pro_error_struct = zeros(5,M,K);
for m = 1:M
    %% Formation control Noiseless case
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
    
    % Add the error curve
    Pro_error_struct(1,m,:) = pos_err;

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
    
    % Add the error curve
    Pro_error_struct(2,m,:) = pos_err;

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
     
    % Add the error curve
    Pro_error_struct(3,m,:) = pos_err;

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
    
    T = 10;
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
    
    % Add the error curve
    Pro_error_struct(4,m,:) = pos_err;

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
    
    % Add the error curve
    Pro_error_struct(5,m,:) = pos_err;
    
    % Loading bar progress
    m
end

%% Plotting
plot_noiseless_case = reshape(mean(Pro_error_struct(1,:,:),2),K,1);
plot_noise_case = reshape(mean(Pro_error_struct(2,:,:),2),K,1);
plot_MA = reshape(mean(Pro_error_struct(3,:,:),2),K,1);
plot_MLE = reshape(mean(Pro_error_struct(4,:,:),2),K,1);
plot_MLE_100 = reshape(mean(Pro_error_struct(5,:,:),2),K,1);

time = linspace(0,maxTime,K);
figure

plot(time,plot_noiseless_case, 'DisplayName','Noiseless case');
hold on;
plot(time,plot_noise_case, 'DisplayName','Noise case');
hold on;
plot(time,plot_MA, 'DisplayName','MA (LP size=10)');
hold on;
plot(time,plot_MLE, 'DisplayName','MLE (T=10)');
hold on;
plot(time,plot_MLE_100, 'DisplayName','MLE (T=100)');
hold on;

yscale("log")
grid("on")
ylabel("Procrutes error")
xlabel("Time [s]")
title("Mean convergence over 50 runs")
legend("AutoUpdate","on")