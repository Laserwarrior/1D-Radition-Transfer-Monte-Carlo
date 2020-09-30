clear
clc

[CrossSection_DataSource, Surface_Sources_type, Spatial_distribution, Number_of_histories_test, Tallies] = Input(); %% can be designed by manually input
[Macroscopic_Cross_Sections, Surface_Sources] = CrossSection(CrossSection_DataSource); 

Zeroth_spatial_mement_scalar_different_NumOfHis = zeros(length(Spatial_distribution.Thickness), length(Number_of_histories_test) );
First_angular_mement_different_NumOfHis = zeros(length(Spatial_distribution.Thickness), length(Number_of_histories_test) );

for i = 1:length(Number_of_histories_test) 
    
    Number_of_histories = Number_of_histories_test(i);
    [Spatial_Material_Parameter, Strength_Probablity_kind, Strength_Probablity_kind_effective, Zeroth_spatial_mement_scalar_flux, First_angular_mement, Current_on_the_right_surface, First_spatial_moment, stand_deviation_Zeroth_spatial_mement_scalar_flux, stand_deviation_First_angular_moment] = Solver(Macroscopic_Cross_Sections, Surface_Sources, Surface_Sources_type, Spatial_distribution, Number_of_histories, Tallies);
    
    Zeroth_spatial_mement_scalar_different_NumOfHis(:,i) =  Zeroth_spatial_mement_scalar_flux;
    First_angular_mement_different_NumOfHis(:,i) = First_angular_mement;
    
end

% % value of the moment vs. number of simulation times
for i = 1:length(Spatial_distribution.Thickness)
    figure(1);plot(Number_of_histories_test, Zeroth_spatial_mement_scalar_different_NumOfHis(i,:),'LineWidth',3);hold on
    figure(2);plot(Number_of_histories_test, abs(First_angular_mement_different_NumOfHis(i,:)),'LineWidth',3);hold on
end
figure(1);legend(Spatial_distribution.Material)
figure(2);legend(Spatial_distribution.Material)

% %% plot the the results for the last entry of the Number_of_histories_test
% % x label is the material type
% X = categorical(Spatial_Material_Parameter.Row);
% X = reordercats(X,Spatial_Material_Parameter.Row);
% figure(1);bar(X, Zeroth_spatial_mement_scalar_different_NumOfHis(:,length(Number_of_histories_test)))
% figure(2);bar(X, First_angular_mement_different_NumOfHis(:,length(Number_of_histories_test) ))
% figure(3);bar(X, stand_deviation_Zeroth_spatial_mement_scalar_flux)
% figure(4);bar(X, stand_deviation_First_angular_moment)

% % x label is the x coordinate
% figure(1);bar(Spatial_Material_Parameter.Material_right_boundary, Zeroth_spatial_mement_scalar_different_NumOfHis(:,length(Number_of_histories_test)))
% figure(2);bar(Spatial_Material_Parameter.Material_right_boundary, First_angular_mement_different_NumOfHis(:,length(Number_of_histories_test) ))
% 
% figure(3);bar(Spatial_Material_Parameter.Material_right_boundary, stand_deviation_Zeroth_spatial_mement_scalar_flux)
% figure(4);bar(Spatial_Material_Parameter.Material_right_boundary, stand_deviation_First_angular_moment)