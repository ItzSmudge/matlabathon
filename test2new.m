clc; close all; 
load('ExampleData.mat');
lat = Position.latitude;
lon = Position.longitude;
positionDatetime = Position.Timestamp;

Xacc = Acceleration.X;
Yacc = Acceleration.Y;
Zacc = Acceleration.Z;
accelDatetime = Acceleration.Timestamp;

% Convert datetime to seconds
positionTime = timeElapsed(positionDatetime);
accelTime = timeElapsed(accelDatetime);

% Earth's circumference in METERS (not miles)
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
% Define filter coefficients for a basic low-pass filter
alpha = 0.25;  % 0 < alpha < 1 (lower = smoother)
b = alpha;
a = [1, -(1 - alpha)];

% Apply filtering to sitAcceleration
sitAcceleration.X = filter(b, a, sitAcceleration.X);
sitAcceleration.Y = filter(b, a, sitAcceleration.Y);
sitAcceleration.Z = filter(b, a, sitAcceleration.Z);

sitLabel = 'sitting';
sitLabel = repmat(sitLabel, size(sitAcceleration, 1), 1);
sitAcceleration.Activity = sitLabel;

% Apply filtering to walkAcceleration
walkAcceleration.X = filter(b, a, walkAcceleration.X);
walkAcceleration.Y = filter(b, a, walkAcceleration.Y);
walkAcceleration.Z = filter(b, a, walkAcceleration.Z);

walkLabel = 'walking';
walkLabel = repmat(walkLabel, size(walkAcceleration, 1), 1);
walkAcceleration.Activity = walkLabel;

% Apply filtering to runAcceleration
runAcceleration.X = filter(b, a, runAcceleration.X);
runAcceleration.Y = filter(b, a, runAcceleration.Y);
runAcceleration.Z = filter(b, a, runAcceleration.Z);

runLabel = 'running';
runLabel = repmat(runLabel, size(runAcceleration, 1), 1);
runAcceleration.Activity = runLabel;

allAcceleration = [sitAcceleration; walkAcceleration; runAcceleration];
allAcceleration = timetable2table(allAcceleration, "ConvertRowTimes", false);

%%
unknownAcceleration % preview data to ensure the correct format

%%
justAcc = timetable2table(unknownAcceleration, "ConvertRowTimes", false);
yfit = trainedModel5.predictFcn(justAcc)

figure;
plot(accelTime,Xacc);
hold on;
plot(accelTime,Yacc); 
plot(accelTime,Zacc-9);
hold off

figure;
plot(accelTime,Xacc);
hold on;
plot(accelTime,Yacc); 
plot(accelTime,Zacc-9);
xlim([0 50])
hold off

figure;
plot(accelTime,Xacc);
hold on;
plot(accelTime,Yacc); 
plot(accelTime,Zacc-9);
xlim([0 50])
legend('X Acceleration','Y Acceleration','Z Aceeleration');
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)');
title('Acceleration Data Vs. Time');
hold off

yfitcat = categorical(cellstr(yfit)); % convert to correct data type for histogram
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
