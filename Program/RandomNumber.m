function [Omega_z] = RandomNumber(CrossSection_Parameter)


%% Probablity_Density = 1/2 +  Sigma_s_0 * Omega_z / ( 2* Sigma_s_1) 

CrossSection_coef = CrossSection_Parameter{1,2}/CrossSection_Parameter{1,1};
Omega_z = (-0.5 + sqrt( 1/4 - CrossSection_coef * (0.5 - 0.25*CrossSection_coef - rand ) )) / (0.5 * CrossSection_coef)  ; % inverse transform sampling

end