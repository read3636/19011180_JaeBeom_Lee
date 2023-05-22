
arg_prg = input('Argument of perigee (rad): ');
inc_angle = input('Inclination (rad): ');
RAAN = input('Right Ascension of the Ascending Node (rad): ');

disp(PQW2ECI(arg_prg,inc_angle,RAAN));


function [rotation_matrix]=PQW2ECI(arg_prg,inc_angle,RAAN)

rotation_matrix=[cos(RAAN)*cos(arg_prg)-sin(RAAN)*cos(inc_angle)*sin(arg_prg), -cos(RAAN)*sin(arg_prg)-sin(RAAN)*cos(inc_angle)*cos(arg_prg), sin(RAAN)*sin(inc_angle);
    sin(RAAN)*cos(arg_prg)+cos(RAAN)*cos(inc_angle)*sin(arg_prg), -sin(RAAN)*sin(arg_prg)+cos(RAAN)*cos(inc_angle)*cos(arg_prg), -cos(RAAN)*sin(inc_angle);
    sin(inc_angle)*sin(arg_prg), sin(inc_angle)*cos(arg_prg), cos(inc_angle);];

end