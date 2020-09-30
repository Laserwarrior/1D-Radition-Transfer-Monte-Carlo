function [Spatial_Material_Parameter, Strength_Probablity_kind, Strength_Probablity_kind_effective, Zeroth_spatial_mement_scalar_flux, First_angular_mement, Current_on_the_right_surface, First_spatial_moment, stand_deviation_Zeroth_spatial_mement_scalar_flux, stand_deviation_First_angular_moment] = Solver(Macroscopic_Cross_Sections, Surface_Sources, Surface_Sources_type, Spatial_distribution, Number_of_histories, Tallies)

%% Determine the probablity to produce the surface source or volumetric source by the strength associated

% Determine the strength of the surface Source
if ~ strcmp(Surface_Sources_type, 'Vacuum')
    M = Surface_Sources(Surface_Sources_type,:);

    Range_of_Omega_z = char( M.(1) ); % str2func(char(Surface_Sources{'Isotropic','Function'}))
    phi = str2func( char ( M.(2) ) ); % str2func(char(Surface_Sources{'Isotropic','Function'}))

    if strcmp(Range_of_Omega_z, 'None')
        Omega_z_min = -1;
        Omega_z_max = 1;
    else
        larger_position  = strfind(Range_of_Omega_z,'>'); % can be modified later to be more general
        Omega_z_min = str2double( Range_of_Omega_z(larger_position + 1:end) );
        Omega_z_max = 1;
    end

    if strcmp(Surface_Sources_type, 'Discrete') % Can be modified later. Replace it by integrating the dirac function
        Strength_Surface_Source = 1;
    else
        Strength_Surface_Source = integral (phi, Omega_z_min, Omega_z_max,'ArrayValued',true);
    end
end
% Determine the strength of the volumertric strength
Spatial_Material  = Spatial_distribution.Material;
Spatial_Thickness = Spatial_distribution.Thickness;
Spatial_Material_Parameter = Macroscopic_Cross_Sections(Spatial_Material,:);

Spatial_Material_Parameter.Thickness = Spatial_Thickness; % build the material system with material, thickness, crosssection, and source
Material_right_boundary = cumsum(Spatial_Material_Parameter.Thickness);
Spatial_Material_Parameter.Material_right_boundary = Material_right_boundary; % build the coordinate
Material_End = Spatial_Material_Parameter.Material_right_boundary(end);
Material_Kind = length(Spatial_Thickness);

if ~ strcmp(Surface_Sources_type, 'Vacuum')
    Total_Strength = Strength_Surface_Source + sum(Spatial_Material_Parameter.Q_0);
    Source_Strength_Probablity = [Strength_Surface_Source;Spatial_Material_Parameter.Q_0]/Total_Strength;
    Source_type = [{Surface_Sources_type};Spatial_Material_Parameter.Row];
    Surface_source_flag = 1;
else
    Total_Strength = sum(Spatial_Material_Parameter.Q_0);
    Source_Strength_Probablity = Spatial_Material_Parameter.Q_0 / Total_Strength;
    Source_type = Spatial_Material_Parameter.Row;
    Surface_source_flag = 0;
end

Strength_Probablity_kind = table(Source_Strength_Probablity, 'RowNames', Source_type); % Probablity for each type of source
Probablity_larger_than_zero_index = Strength_Probablity_kind.(1) > 0;
Strength_Probablity_kind_effective = Strength_Probablity_kind(Probablity_larger_than_zero_index,:);
Strength_Probablity_effective = Strength_Probablity_kind_effective.(1);

Travel_Distance_magnitude_N_times_Simulation = zeros(Material_Kind, Number_of_histories);
Travel_Z_direction_N_times_Simulation = zeros(Material_Kind, Number_of_histories);
First_spatial_moment_N_times_Simulation = zeros(Material_Kind, Number_of_histories);
Current_on_the_right_N_times_Simulation = zeros(1, Number_of_histories);

% Travel_Distance_magnitude_N_times_Simulation = zeros(Material_Kind, 1);
% Travel_Z_direction_N_times_Simulation = zeros(Material_Kind, 1);
% First_spatial_moment_N_times_Simulation = zeros(Material_Kind, 1);
% Current_on_the_right_N_times_Simulation = zeros(1, 1);


% stand_deviation_Zeroth_spatial_mement_scalar_flux = zeros(Material_Kind, 1);
% stand_deviation_fist_angle_moment = zeros(Material_Kind, 1);
% stand_deviation_First_spatial_moment_scalar_flux = zeros(Material_Kind, 1);
% stand_deviation_Current_on_right_surface = zeros(Material_Kind, 1);
%% Simulate N number of particle history

for i = 1: Number_of_histories % Strength is included in the one particle history
    
    [Travel_Distance_magnitude_total_Multiple_Run, Travel_Z_direction_total_total_Multiple_Run, Current_on_right_total, First_spatial_moment_Multiple_Run] = ParticleHistory(Surface_source_flag, Material_Kind, Surface_Sources_type, Surface_Sources, Spatial_Material_Parameter, Material_End, Strength_Probablity_kind_effective, Strength_Probablity_effective);
    
    Travel_Distance_magnitude_N_times_Simulation(:, i) = Travel_Distance_magnitude_total_Multiple_Run; % check right to be column
    Travel_Z_direction_N_times_Simulation(:, i) = Travel_Z_direction_total_total_Multiple_Run;
%     Current_on_the_right_N_times_Simulation(i) = Current_on_right_total;
%     First_spatial_moment_N_times_Simulation(:, i) = First_spatial_moment_Multiple_Run; % only work for single material

%       Travel_Distance_magnitude_N_times_Simulation = Travel_Distance_magnitude_N_times_Simulation + Travel_Distance_magnitude_total_Multiple_Run;
%       Travel_Z_direction_N_times_Simulation = Travel_Z_direction_N_times_Simulation + Travel_Z_direction_total_total_Multiple_Run;
    
    
%     figure(1);scatter(i,Travel_Distance_magnitude_total_Multiple_Run/100); hold on;
%     figure(2);scatter(i,Current_on_right_total); hold on;
%     figure(3);scatter(i,First_spatial_moment_Multiple_Run); hold on;
%     
end

Zeroth_spatial_mement_scalar_flux = mean(Travel_Distance_magnitude_N_times_Simulation, 2) ./ Spatial_Material_Parameter.Thickness;
First_angular_mement = mean(Travel_Z_direction_N_times_Simulation, 2)  ./ Spatial_Material_Parameter.Thickness;
Current_on_the_right_surface = mean(Current_on_the_right_N_times_Simulation);
First_spatial_moment = First_spatial_moment_N_times_Simulation / Number_of_histories; % only work for single material

stand_deviation_Zeroth_spatial_mement_scalar_flux = std(Travel_Distance_magnitude_N_times_Simulation./Spatial_Material_Parameter.Thickness, 0, 2);
stand_deviation_First_angular_moment = std(Travel_Z_direction_N_times_Simulation./Spatial_Material_Parameter.Thickness, 0, 2);

% Zeroth_spatial_mement_scalar_flux = Travel_Distance_magnitude_N_times_Simulation / Number_of_histories  ./Spatial_Material_Parameter.Thickness;
% First_angular_mement = Travel_Z_direction_N_times_Simulation / Number_of_histories  ./Spatial_Material_Parameter.Thickness;
% Current_on_the_right_surface = 0; % no meaning, place holdeer
% First_spatial_moment = 0; % not meaning, place holder
% 
% stand_deviation_Zeroth_spatial_mement_scalar_flux = 0; %place holder
% stand_deviation_First_angular_moment = 0; %place holder
% stand_deviation_Current_on_right_surface = stand_deviation_Current_on_right_surface + (Current_on_right_total - 0).^2;

end