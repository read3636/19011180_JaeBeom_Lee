
%{
HW11

arg_prg = input('Argument of perigee (degree): ');
inc_angle = input('Inclination (degree): ');
RAAN = input('Right Ascension of the Ascending Node (degree): ');

disp(PQW2ECI(arg_prg,inc_angle,RAAN));
%}

semimajor_axis=input('semimajor_axis (Km): ');
eccentricity=input('ecentricity : ');
true_anomaly=input('true_anomaly (degree): ');

fprintf('rangeInPQW = %m', solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly));
fprintf('velocityInPQW = %m', solveVelocityInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly));


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
function [rangeinPQW]=solveRangeInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly)  % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly=deg2rad(true_anomaly);
p=semimajor_axis*(1-eccentricity^2);
r=p/(1+eccentricity*cos(true_anomaly));

rangeinPQW=[r*cos(true_anomaly); r*sin(true_anomaly); 0;];


end

function [velocityInPQW]=solveVelocityInPerifocalFrame(semiajor_axis,eccentricity,true_anomaly) % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly=deg2rad(true_anomaly);
p=semiajor_axis*(1-eccentricity^2);
GM = 3.986004418*10^5; % μ = 3.986004418 × 10^14 [m^3 s^(−2)] --> μ = 3.986004418 × 10^5 [km^3 s^(−2)]

velocityInPQW=sqrt(GM/p)*[-sin(true_anomaly); eccentricity+cos(true_anomaly); 0;];

end