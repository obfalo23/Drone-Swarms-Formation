clear all
close all
load("data.mat")

% Plot formation and error
plot_formation(z)

z_pos = zeros(K,N,2);
z_pos(1,:,:) = z;
z = reshape(z_pos(1,:,:), size(z));

pos_err = double.empty;
pos_err = [pos_err; norm(z-z_star,2)];

trajectory = double.empty;
trajectory = [trajectory; z_pos(1,7,:)];

% Formation control noiseless case
% Initialization
U = zeros(K,N,2);
for k = 1:K
    for i = 4:N
        v = randn(size(z))*R;
        z_i = reshape(z_pos(k,i,:), size(z(i,:)));
        U(k,i,:) = (L(i,:)*(z_i-z+v));
        z_pos(k+1,i,:) = z_pos(k,i,:) + 10*dt*U(k,i,:);
        z(i,:) = reshape(z_pos(k+1,i,:), size(z(i,:)));
    end
    pos_err = [pos_err; norm(z-z_star,2)];
    trajectory = [trajectory; z_pos(k,7,:)];
end
plot_formation(z);

figure
plot(pos_err)
disp("Final error")
disp(pos_err(end))