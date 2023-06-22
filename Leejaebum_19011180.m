load("nav.mat");
r_location=zeros(1,3);
r_location(1,1)=input('longitude (degree): ');
r_location(2,1)=input('latitude (degree): ');
r_location(3,1)=input('height (km): ');
n=input("choose satellite (gps:1  QZSS:2  BDS:3: ");

if n == 1
    satellite = nav.GPS;
elseif n == 2
    satellite = nav.QZSS;
else
    satellite = nav.BDS;
end

a = satellite.a/10^3;
e = satellite.e;
i = satellite.i;
omega = satellite.omega;
M0 = satellite.M0;
OMEGA = satellite.OMEGA;

if M0 < 0
    M0 = M0 + 2*pi;
end

t = 0:60:86399; % 0부터 86399까지 1초 간격으로 시간 범위 설정

r_perifocal = zeros(3, length(t));
r_eci = zeros(3, length(t));
r_ecef = zeros(3, length(t));
r_enu = zeros(3, length(t));
r_longitude = zeros(1, length(t));
r_latitude = zeros(1, length(t));
r_azmuth = zeros(1, length(t));
r_elevation = zeros(1, length(t));

for p = 1:1:length(t)

    M = M0 + mean_motion(a) * t(p);
    M = rad2deg(M);
    M = mod(M, 360);
    
    second = rem(t(p), 60);
    minute = floor(t(p)/60);
    minute = mod(minute, 60);
    hour = floor(t(p)/3600);
    toc = satellite.toc;
    toc(1,4) = hour;
    toc(1,5) = minute;
    toc(1,6) = second;

    r_perifocal(:,p) = solveRangeInPerifocalFrame(a, e, mean2true(M,e));
    r_eci(:,p) = PQW2ECI(omega, i, OMEGA)*r_perifocal(:,p);
    r_ecef(:,p) = ECI2ECEF_DCM(toc)*r_eci(:,p);
    r_location_ecef = [(6371+r_location(3,1))*cosd(r_location(1,1))*sind(r_location(2,1)); (6371+r_location(3,1))*cosd(r_location(2,1))*sind(r_location(1,1)); (6371+r_location(3,1))*sind(r_location(2,1))];
    r_enu(:,p) = eci2enu(r_location(1,1), r_location(2,1), toc)*(r_eci(:,p)-ECI2ECEF_DCM(toc)'*r_location_ecef);

    r_longitude(1,p) = rad2deg(lon(r_ecef(:,p)));
    r_latitude(1,p) = rad2deg(lat(r_ecef(:,p)));
    r_azmuth(1,p) = azimuth(r_enu(:,p));
    r_elevation(1,p) = elevation(r_enu(:,p),r_location(2,1));
end   
    subplot(3,2,1), plot3(r_perifocal(1,:), r_perifocal(2,:), r_perifocal(3,:)),title('perifocal frame');
    subplot(3,2,2), plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:)),title('eci frame');
    subplot(3,2,3), plot3(r_ecef(1,:), r_ecef(2,:), r_ecef(3,:)),title('ecef frame');
    subplot(3,2,4), plot3(r_enu(1,:), r_enu(2,:), r_enu(3,:)),title('enu frame');
    subplot(3,2,5), geoplot(r_latitude, r_longitude),title('ground track');
    subplot(3,2,6), polarplot(deg2rad(r_azmuth), r_elevation),title('sky view');

%% PQW to ECI 000
function [rotation_matrix]=PQW2ECI(arg_prg,inc_angle,RAAN) 

arg_prg=deg2rad(arg_prg);
inc_angle=deg2rad(inc_angle);
RAAN=deg2rad(RAAN);

rotation_matrix=[cos(RAAN)*cos(arg_prg)-sin(RAAN)*cos(inc_angle)*sin(arg_prg), -cos(RAAN)*sin(arg_prg)-sin(RAAN)*cos(inc_angle)*cos(arg_prg), sin(RAAN)*sin(inc_angle);
                 sin(RAAN)*cos(arg_prg)+cos(RAAN)*cos(inc_angle)*sin(arg_prg), -sin(RAAN)*sin(arg_prg)+cos(RAAN)*cos(inc_angle)*cos(arg_prg), -cos(RAAN)*sin(inc_angle);
                 sin(inc_angle)*sin(arg_prg), sin(inc_angle)*cos(arg_prg), cos(inc_angle);];

end

%% r_PQW  000
function [rangeInPQW]=solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly)  % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly=deg2rad(true_anomaly);
p=semimajor_axis*(1-eccentricity^2);
r=p/(1+eccentricity*cos(true_anomaly));

rangeInPQW=[r*cos(true_anomaly); r*sin(true_anomaly); 0;];

end

%% eci to ecef
function [DCM] = ECI2ECEF_DCM(time)

standard_time = [2000,1,1,12,0,0];
standard_time = datetime(standard_time);
standard_jd = juliandate(standard_time);

time = datetime(time);
jd = juliandate(time);
GMST = mod(280.4606 + 360.9856473*(jd - standard_jd), 360);
GMST_rad = deg2rad(GMST);

DCM = [cos(GMST_rad), sin(GMST_rad), 0;
       -sin(GMST_rad), cos(GMST_rad), 0;
       0, 0, 1];
end

%% eci to enu
function [rotation_matrix] = eci2enu(longitude, latitude, toc)
    standard_time = [2000,1,1,12,0,0];
    standard_time = datetime(standard_time);
    standard_jd = juliandate(standard_time);

    toc = datetime(toc);
    jd = juliandate(toc);
    GMST = mod(280.4606 + 360.9856473*(jd - standard_jd), 360);
    GMST_rad = deg2rad(GMST);

    rotation_matrix = [-sin(GMST_rad-longitude), cos(GMST_rad-longitude), 0;
                       -sin(latitude)*cos(GMST_rad-longitude), -sin(latitude)*sin(GMST_rad-longitude), cos(latitude);
                       cos(latitude)*cos(GMST_rad-longitude), cos(latitude)*sin(GMST_rad-longitude), sin(latitude)];
end


%% latitude
function [latitude] = lat(ecef)

magnitude_xy = sqrt(ecef(1,:).^2 + ecef(2,:).^2);
latitude = zeros(size(ecef(1,:)));
zero_idx = abs(magnitude_xy) < 1e-8;
latitude(zero_idx) = sign(ecef(3, zero_idx)) * pi/2;
nonzero_idx = ~zero_idx;
latitude(nonzero_idx) = 50*atan2(ecef(3, nonzero_idx), magnitude_xy(nonzero_idx));
end


%% longitude

function [longitude] = lon(ecef)
longitude = atan2(ecef(2,:), ecef(1,:));
end


%% az
function az = azimuth(enu)

az = rad2deg(atan2(enu(2,:), enu(1,:)));
az = mod(az, 360);

end

%% el
function el = elevation(r_enu, el_mask)
el = rad2deg(atan2(r_enu(3,:), sqrt(r_enu(1,:).^2 + r_enu(2,:).^2)));
el(el < el_mask) = NaN;
end

%% mean anomaly to true anomaly
function [true_anomaly] = mean2true(mean_anomaly, eccentricity)
mean_anomaly = deg2rad(mean_anomaly);
eccentricity_anomaly = mean_anomaly;
delta = 1;

while abs(delta) > 10^(-8)
    delta = (eccentricity_anomaly - eccentricity * sin(eccentricity_anomaly) - mean_anomaly) / (1 - eccentricity * cos(eccentricity_anomaly));
    eccentricity_anomaly = eccentricity_anomaly - delta;
end

true_anomaly = atan2(sqrt(1 - eccentricity^2) * sin(eccentricity_anomaly), cos(eccentricity_anomaly) - eccentricity);
true_anomaly = rad2deg(true_anomaly);
end

%% mean motion
function [mean_motion] = mean_motion(semimajor_axis)

mu = 3.986004418*10^5;
mean_motion = sqrt(mu / semimajor_axis^3);

end
