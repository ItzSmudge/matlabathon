clc; close all; 
load("ExampleData.mat");
%%
lat = Position.latitude;
lon = Position.longitude;
positionDatetime = Position.Timestamp;

accelDatetime = Acceleration.Timestamp;

% Convert datetime to seconds
positionTime = timeElapsed(positionDatetime);
accelTime = timeElapsed(accelDatetime);

% Earth's circumference in meters
earthCirc = 40075000;  % Equatorial circumference in meters
totaldis = 0; 

% Calculate total distance in meters
for i = 1:(length(lat)-1)
    lat1 = lat(i);     % The first latitude
    lat2 = lat(i+1);   % The second latitude
    lon1 = lon(i);     % The first longitude
    lon2 = lon(i+1);   % The second longitude
    degDis = distance(lat1, lon1, lat2, lon2); % Distance in degrees
    dis = (degDis/360) * earthCirc; % Convert degrees to meters
    totaldis = totaldis + dis;
end

% Stride length in meters (0.762m is ~2.5 ft, a typical average)
stride = 0.762; 

% Calculate steps directly from meters
steps = totaldis / stride;

% Display results in METERS
disp(['The total distance traveled is: ', num2str(totaldis), ' meters'])
disp(['You took ', num2str(steps), ' steps'])
%%
% convert position to x,y,z
latitude = Position.latitude;
longitude = Position.longitude; 
altitude = Position.altitude;
position_time = Position.Timestamp; 

origin = [latitude(1), longitude(1), altitude(1)];
[x,y,z] = latlon2local(latitude, longitude, altitude, origin);

figure(1);
subplot(3,1,1);
plot(position_time, x, "LineWidth",1.5);
xlabel("time");
ylabel("x(m)"); 
ax = gca;
ax.LineWidth = 1.5; 
ax.FontSize = 10; 

subplot(3,1,2);
plot(position_time, y, "LineWidth",1.5);
xlabel("time");
ylabel("y(m)"); 
ax = gca; 
ax.LineWidth = 1.5; 
ax.FontSize = 10; 

subplot(3,1,3);
plot(position_time, z, "LineWidth",1.5);
xlabel("time");
ylabel("z(m)"); 
ax = gca; 
ax.LineWidth = 1.5; 
ax.FontSize = 10; 


xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
ax = gca; 
ax.LineWidth = 1.5; 
ax.FontSize = 10; 
title('3D Trajectory');
grid on;
%axis equal;
hold off;

figure(2);
plot3(x, y, z, 'b-', "LineWidth",1.5);
hold on;
plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot3(x(end), y(end), z(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');


%% calculate total distancce 
dx = diff(x);
dy = diff(y); 
dz = diff(z);


distances = sqrt(dx.^2+dy.^2+dz.^2);

total_distance = sum(distances);


%% velocity calculation
dt = seconds(diff(position_time));
vx = diff(x);
vy = diff(y);
vz = diff(z); 

% 0 < alpha <= q
alpha = 0.3; % higher alpha = less smoothing, lower alpha = more smoothing


vx_filtered = zeros(size(vx));
vy_filtered = zeros(size(vy));
vz_filtered = zeros(size(vz));

% Apply Exponential Moving Average filter
vx_filtered(1) = vx(1); % Initialize with first value
vy_filtered(1) = vy(1);
vz_filtered(1) = vz(1);

for i = 2:length(vx)
    vx_filtered(i) = alpha * vx(i) + (1 - alpha) * vx_filtered(i-1);
    vy_filtered(i) = alpha * vy(i) + (1 - alpha) * vy_filtered(i-1);
    vz_filtered(i) = alpha * vz(i) + (1 - alpha) * vz_filtered(i-1);
end

velocities = sqrt(vx.^2 + vy.^2 + vz.^2);
total_velocity = mean(velocities);

velocities_filtered = sqrt(vx_filtered.^2 + vy_filtered.^2 + vz_filtered.^2);
total_velocity_filtered = mean(velocities_filtered);



figure(3);
%x velocity
subplot(3,1,1);
plot(position_time(2:end), vx, 'b','DisplayName', 'Raw');
hold on;
plot(position_time(2:end), vx_filtered, 'r', 'DisplayName', 'Filtered', "LineWidth",1);
hold off;
ylabel('V_x (m/s)');
title('Velocity Components and Magnitude');
ax = gca; %
ax.LineWidth = 1.5; 
ax.FontSize = 10; 

legend('Location', 'best');
grid on;

% y velocity
subplot(3,1,2);
plot(position_time(2:end), vy, 'b','DisplayName', 'Raw',"LineWidth",1.5);
hold on;
plot(position_time(2:end), vy_filtered, 'r', 'DisplayName', 'Filtered',"LineWidth",1.5);
hold off;
ylabel('V_y (m/s)');
ax = gca;
ax.LineWidth = 1.5; 
ax.FontSize = 10; 

legend('Location', 'best');
grid on;

%z velocity
subplot(3,1,3);
plot(position_time(2:end), vz, 'b','DisplayName', 'Raw',"LineWidth",1.5);
hold on;
plot(position_time(2:end), vz_filtered, 'r', 'DisplayName', 'Filtered',"LineWidth",1.5);
hold off;
ylabel('V_z (m/s)');
ax = gca; 
ax.LineWidth = 1.5;
ax.FontSize = 10;

legend('Location', 'best');
grid on;

figure(4);
plot(position_time(2:end), sqrt(vx.^2 + vy.^2 + vz.^2), 'b',"LineWidth",1.5,'DisplayName', 'Raw');
hold on;
plot(position_time(2:end), sqrt(vx_filtered.^2 + vy_filtered.^2 + vz_filtered.^2), 'r',"LineWidth",1.5,'DisplayName', 'Filtered');
hold off; 
%plot(position_time(2:end), total_velocity, 'r');
ylim([0 100]);
ylabel('Speed (m/s)');
xlabel('Time');
ax = gca; % Get current axes
ax.LineWidth = 1.5; % Make axis lines thicker
ax.FontSize = 10; % Optional: Increase axis label font size
legend('Location', 'best');
grid on;



%% anomaly detection 

% Check for NaN values
nan_values = isnan(velocities_filtered);
if any(nan_values)
    fprintf('Found %d NaN values at indices: %s\n', ...
            sum(nan_values), mat2str(find(nan_values)));
    %
    velocities_filtered(nan_values) = interp1(find(~nan_values), ...
                                          velocities_filtered(~nan_values), ...
                                          find(nan_values));
end

% Detect velocity jumps (>15 m/s then change)
jump_threshold = 30;
velocity_jumps = find(abs(diff(velocities_filtered)) > jump_threshold);
if ~isempty(velocity_jumps)
    fprintf('Found %d velocity jumps exceeding %d m/s at:\n', ...
            length(velocity_jumps), jump_threshold);
    disp(arrayfun(@(x) position_time(x+1), velocity_jumps, 'UniformOutput', false));
    
    %median filtering
    for jump = velocity_jumps'
        window = max(1,jump-3):min(length(velocities_filtered),jump+3);
        velocities_filtered(jump) = median(velocities_filtered(window));
    end
end

% Detect gps anomalies
gps_jump_threshold = 30; 
gps_jumps = find(distances > gps_jump_threshold);
if ~isempty(gps_jumps)
    fprintf('Found %d GPS position jumps > %.1f meters at:\n', ...
            length(gps_jumps), gps_jump_threshold);
    disp(arrayfun(@(x) position_time(x+1), gps_jumps, 'UniformOutput', false));
    
end

%%
justAcc = timetable2table(Acceleration, "ConvertRowTimes", false);
yfit = trainedModelNN.predictFcn(justAcc);

figure;
plot(accelTime,a_x);
hold on;
plot(accelTime,a_y); 
plot(accelTime,a_z-9);
hold off

figure;
plot(accelTime,a_x);
hold on;
plot(accelTime,a_y); 
plot(accelTime,a_z-9);
xlim([0 50])
hold off

figure;
plot(accelTime,a_x);
hold on;
plot(accelTime,a_y); 
plot(accelTime,a_z-9);
xlim([0 50])
legend('X Acceleration','Y Acceleration','Z Aceeleration');
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)');
title('Acceleration Data Vs. Time');
hold off

yfitcat = categorical(cellstr(yfit));
hist(yfitcat)
pie(yfitcat)


% Plot filtered acceleration for sitting
figure;
subplot(3,1,1);
plot(sitAcceleration.X);
title('Filtered X - Sitting');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,2);
plot(sitAcceleration.Y);
title('Filtered Y - Sitting');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,3);
plot(sitAcceleration.Z);
title('Filtered Z - Sitting');
xlabel('Sample');
ylabel('Acceleration');

% Plot filtered acceleration for walking
figure;
subplot(3,1,1);
plot(walkAcceleration.X);
title('Filtered X - Walking');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,2);
plot(walkAcceleration.Y);
title('Filtered Y - Walking');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,3);
plot(walkAcceleration.Z);
title('Filtered Z - Walking');
xlabel('Sample');
ylabel('Acceleration');

% Plot filtered acceleration for running
figure;
subplot(3,1,1);
plot(runAcceleration.X);
title('Filtered X - Running');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,2);
plot(runAcceleration.Y);
title('Filtered Y - Running');
xlabel('Sample');
ylabel('Acceleration');

subplot(3,1,3);
plot(runAcceleration.Z);
title('Filtered Z - Running');
xlabel('Sample');
ylabel('Acceleration');


%% Calorie Calculation for Each Activity

% Define your weight in kilograms
weightKg = 60;  % <- you can change this to match your weight

% Get the number of samples and sampling rate (assuming consistent rate)
Fs = 1;  % 1 Hz sample rate
sampleDuration = 1 / Fs;  % seconds per sample

% Compute durations in seconds
dur_sit = height(sitAcceleration) * sampleDuration;
dur_walk = height(walkAcceleration) * sampleDuration;
dur_run = height(runAcceleration) * sampleDuration;

% Convert durations to hours
dur_sit_hr = dur_sit / 3600;
dur_walk_hr = dur_walk / 3600;
dur_run_hr = dur_run / 3600;

% MET values (approximate)
MET_sit = 1.3;
MET_walk = 3.5;
MET_run = 7.0;

% Calorie burn formula: Calories = MET × Weight (kg) × Duration (hr)
cal_sit = MET_sit * weightKg * dur_sit_hr;
cal_walk = MET_walk * weightKg * dur_walk_hr;
cal_run = MET_run * weightKg * dur_run_hr;

% Display results
fprintf("Calories burned while sitting: %.2f kcal\n", cal_sit);
fprintf("Calories burned while walking: %.2f kcal\n", cal_walk);
fprintf("Calories burned while running: %.2f kcal\n", cal_run);

% Optional pie chart
figure;
cal_data = [cal_sit, cal_walk, cal_run];
labels = {'Sitting', 'Walking', 'Running'};
pie(cal_data, labels);
title('Calorie Distribution by Activity');

%% Movement Classification Based on Velocity thresholds
% Define velocity thresholds for different activities
thresholds.sitting = 0.3;    % < 0.3 m/s
thresholds.walking = 1.5;     % 0.3-1.5 m/s
thresholds.running = 4.0;     % 1.5-4.0 m/s
thresholds.fast_run = 10.0;   % > 4.0 m/s (for potential extension)

% Apply additional smoothing to velocity data
smooth_vel = movmean(velocities_filtered, 15); % 15-point moving average

% Classify movement types
modes = strings(length(smooth_vel), 1);
for i = 1:length(smooth_vel)
    v = smooth_vel(i);
    if v < thresholds.sitting
        modes(i) = "sitting";
    elseif v < thresholds.walking
        modes(i) = "walking";
    elseif v < thresholds.running
        modes(i) = "running";
    else
        modes(i) = "fast_run"; % Could be extended for other activities
    end
end


%% Velocity Classification with thresholds (predefined)

% Define thresholds based on estimated velocity ranges (in m/s)
% These may need tuning depending on actual data
thresholds.sitting = 0.5;
thresholds.walking = 2.5;
thresholds.running = 6;

% Smooth velocity using moving average
smooth_vel = movmean(velocities_filtered, 15);

% Classify each time point into activity modes
modes = strings(length(smooth_vel), 1);
for i = 1:length(smooth_vel)
    v = smooth_vel(i);
    if v < thresholds.sitting
        modes(i) = "sitting";
    elseif v < thresholds.walking
        modes(i) = "walking";
    else
        modes(i) = "running";
    end
end

% Define color scheme for plotting
color_map = containers.Map(...
    ["sitting", "walking", "running"], ...
    {[1 0.8 0], [0.2 0.8 0.2], [0 0 0.8]});  % Orange, Green, Blue

% Find transitions between modes
transitions = find([true; ~strcmp(modes(1:end-1), modes(2:end))]);
durations = [diff(transitions); length(modes) - transitions(end) + 1];

% Plot segments by movement mode
figure;
hold on;
for i = 1:length(transitions)
    seg_idx = transitions(i):(transitions(i) + durations(i) - 1);
    plot(position_time(seg_idx+1), velocities_filtered(seg_idx), ...
        'Color', color_map(modes(transitions(i))), 'LineWidth', 2);
end

% Create legend handles
legend_items = [];
mode_order = ["sitting", "walking", "running"];
for i = 1:length(mode_order)
    legend_items(i) = plot(NaN, NaN, 'Color', color_map(mode_order(i)), ...
        'LineWidth', 2, 'DisplayName', mode_order(i));
end

xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend(legend_items, 'Location', 'best');
grid on;

% Format time ticks
ax = gca;
ax.XTick = position_time(1:round(end/10):end);
ax.XTickLabel = datestr(ax.XTick, 'HH:MM');
ax.XTickLabelRotation = 45;
ax.LineWidth = 1.5;
ax.FontSize = 10;

hold off;
