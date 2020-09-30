function [Angle_addition, Strength,flag_2] = ReactionGeneration(Traveled_one_time_Material_Parameter, Strength)


Muliplication_Probablity = Traveled_one_time_Material_Parameter.Sigma_m  /  Traveled_one_time_Material_Parameter.Sigma_t; 
Scattering_Probablity    = Traveled_one_time_Material_Parameter.Sigma_s_0  /  Traveled_one_time_Material_Parameter.Sigma_t; 
Absorption_Probablity = 1 - Muliplication_Probablity - Scattering_Probablity; 

Reaction_Probablity = [Muliplication_Probablity; Scattering_Probablity; Absorption_Probablity];
Reaction_Probablity_region = cumsum(Reaction_Probablity);

% Reaction_Name = {'Muliplication_Probablity'; 'Scattering_Probablity'; 'Absorption_Probablity'};
% Reaction_Parameter = table(Reaction_Probablity, Reaction_Probablity_region, 'RowNames', Reaction_Name);

Right_Reaction_Index = find(  Reaction_Probablity_region >= rand   );
Reaction_Index = Right_Reaction_Index(1);

switch Reaction_Index
    case 1 % multiplication
%         Reaction_Parameter = Traveled_one_time_Material_Parameter(:,{'mu','Sigma_m'});
        Angle_addition = -1 + 2*rand; % Isotropic multiplication (-1, 1)
        Strength = Strength * Traveled_one_time_Material_Parameter.mu;
        flag_2 = 1;
        
    case 2 % scatter
        
        Scattering_Reaction_Parameter = Traveled_one_time_Material_Parameter(:,{'Sigma_s_0','Sigma_s_1'});
        if Scattering_Reaction_Parameter.Sigma_s_1 == 0
            Angle_addition = -1 + 2*rand;
            Strength = Strength * 1;
            flag_2 =1;
        else
        Angle_addition = RandomNumber(Scattering_Reaction_Parameter);
        Strength = Strength * 1;
        flag_2 = 1;
        end
        
    case 3 % absorption
        flag_2 = 0;
        Angle_addition = 0; % dones't have any meaning, because it disappears. Just a holder for output.
        Strength = Strength * 1;
        
    otherwise
        disp('Something wrong in the ReactionGeneration!')
end

end