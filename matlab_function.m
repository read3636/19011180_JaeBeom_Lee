
%{
HW11

arg_prg = input('Argument of perigee (degree): ');
inc_angle = input('Inclination (degree): ');
RAAN = input('Right Ascension of the Ascending Node (degree): ');

disp(PQW2ECI(arg_prg,inc_angle,RAAN));
%}

%{
semimajor_axis=input('semimajor_axis (Km): ');
eccentricity=input('ecentricity : ');
true_anomaly=input('true_anomaly (degree): ');

disp('rangeInPQW =');
disp(solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly));
disp('velocityInPQW =');
disp(solveVelocityInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly));
%}

%% HW12

YYYY = input('YYYY: ');
MM = input('MM: ');
DD = input('DD: ');
hh = input('hh: ');
mm = input('mm: ');
ss = input('ss: ');
time = datetime(YYYY, MM, DD, hh, mm, ss);

a=[1000,2000,3000;
    4000,4000,4000;];

disp(ECI2ECEF_DCM(time));
disp(azimuth(a));
disp(elevation(a,1));


%% HW11
function [rotation_matrix]=PQW2ECI(arg_prg,inc_angle,RAAN)

arg_prg=deg2rad(arg_prg);
inc_angle=deg2rad(inc_angle);
RAAN=deg2rad(RAAN);

rotation_matrix=[cos(RAAN)*cos(arg_prg)-sin(RAAN)*cos(inc_angle)*sin(arg_prg), -cos(RAAN)*sin(arg_prg)-sin(RAAN)*cos(inc_angle)*cos(arg_prg), sin(RAAN)*sin(inc_angle);
                 sin(RAAN)*cos(arg_prg)+cos(RAAN)*cos(inc_angle)*sin(arg_prg), -sin(RAAN)*sin(arg_prg)+cos(RAAN)*cos(inc_angle)*cos(arg_prg), -cos(RAAN)*sin(inc_angle);
                 sin(inc_angle)*sin(arg_prg), sin(inc_angle)*cos(arg_prg), cos(inc_angle);];

end
%% HW12
function [rangeInPQW]=solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly)  % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly=deg2rad(true_anomaly);
p=semimajor_axis*(1-eccentricity^2);
r=p/(1+eccentricity*cos(true_anomaly));

rangeInPQW=[r*cos(true_anomaly); r*sin(true_anomaly); 0;];


end

function [velocityInPQW]=solveVelocityInPerifocalFrame(semiajor_axis,eccentricity,true_anomaly) % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly=deg2rad(true_anomaly);
p=semiajor_axis*(1-eccentricity^2);
GM = 3.986004418*10^5; % μ = 3.986004418 × 10^14 [m^3 s^(−2)] --> μ = 3.986004418 × 10^5 [km^3 s^(−2)]

velocityInPQW=sqrt(GM/p)*[-sin(true_anomaly); eccentricity+cos(true_anomaly); 0;];

end

%% HW13
function [DCM]=ECI2ECEF_DCM(time)

time=datetime(time);
jd= juliandate(time);
GMST = mod(280.4606 + 360.9856473*(jd - 2451545), 360);
GMST_rad = deg2rad(GMST);

    DCM = [cos(GMST_rad), sin(GMST_rad), 0;
           -sin(GMST_rad), cos(GMST_rad), 0;
           0, 0, 1];
end

function az = azimuth(ENU)

    az = rad2deg(acos(ENU(:,2)./sqrt(ENU(:,1).^2 + ENU(:,3).^2)));
    az = mod(az, 360);
end

function el = elevation(ENU, el_mask)

    el = rad2deg(asin(ENU(:,3)./sqrt(ENU(:,1).^2 + ENU(:,2).^2 + ENU(:,3).^2)));
    
    el(el < el_mask) = NaN;
end




