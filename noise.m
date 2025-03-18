n1 = randn(1000,1); % Gaussian noise with mean=0 and std=1

n2 = 5 +randn(1000,1); % Gaussian noise with mean=5 and std=1

n3 = 0.2 * randn(1000,1);% Gaussian noise with mean=0 and std=0.2
mean(n1)

std(n1)
hist(n1, 50)
x = zeros(1000, 1);

y = zeros(1000, 1);

% noise vectors

nx = 0.5 * randn(1000, 1);

ny = 0.5 * randn(1000, 1);

% generate trajectory
function [x, y] = brownianModel()
    %Brownian motion model
% initial position vectors
x = zeros(1000, 1);
y = zeros(1000, 1);
% noise vectors
nx = 0.5 * randn(1000, 1);
ny = 0.5 * randn(1000, 1);
% generate trajectory
for i = 2:1000
x(i) = x(i-1) + nx(i);
y(i) = y(i-1) + ny(i);
end
end
    
  function [x, vx, y, vy] = cvModel()
% generate the 2D trajectory and velocity of a point moving
% according to a Constant Velocity motion model
dt = 0.033; % time interval
% initial position and velocity vectors
x = zeros(1000, 1); % x position
vx = zeros(1000, 1); % x velocity
y = zeros(1000, 1); % y position
vy = zeros(1000, 1); % y velocity
% noise vectors
nx = 0.5 * randn(1000, 1);
nvx = 0.1 * randn(1000, 1);
ny = 0.5 * randn(1000, 1);
nvy = 0.1 * randn(1000, 1);
% generate trajectory
for i = 2:1000
x(i) = x(i-1) + vx(i-1)*dt + nx(i);
vx(i) = vx(i-1) + nvx(i);
y(i) = y(i-1) + vy(i-1)*dt + ny(i);
vy(i) = vy(i-1) + nvy(i);
end
end

[x, y] = brownianModel;
figure(2)
plot(x,y)
figure(3)
plot(x)
figure(4)
plot(y)
figure(5)
plot(x,y)
% Load the data from CSV files
xp = csvread('xp.csv');
yp = csvread('yp.csv');
u = csvread('u.csv');
v = csvread('v.csv');

% Compute the noise components
nu = u - xp;
nv = v - yp;

% Compute the mean and standard deviation of the noise
mean_nu = mean(nu);
std_nu = std(nu);
mean_nv = mean(nv);
std_nv = std(nv);

% Display results
fprintf('Mean of nu: %.4f, Standard Deviation of nu: %.4f\n', mean_nu, std_nu);
fprintf('Mean of nv: %.4f, Standard Deviation of nv: %.4f\n', mean_nv, std_nv);

% Plot histograms of the noise components
figure;
subplot(2,1,1);
histogram(nu, 5);
title('Histogram of nu');
xlabel('Noise in x direction');
ylabel('Frequency');

subplot(2,1,2);
histogram(nv, 5);
title('Histogram of nv');
xlabel('Noise in y direction');
ylabel('Frequency');

% Check if noise follows a zero-mean Gaussian
if abs(mean_nu) < 0.1 * std_nu && abs(mean_nv) < 0.1 * std_nv
    disp('The noise appears to be approximately zero-mean Gaussian.');
else
    disp('The noise does not appear to be zero-mean Gaussian.');
end
