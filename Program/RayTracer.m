function [Travel_Distance_magnitude, Position] = RayTracer(Spatial_Material_Parameter, Current_Material_Parameter, Position, Angle)

Total_CrossSection = Current_Material_Parameter.Sigma_t;
Travel_Distance_magnitude = -1/Total_CrossSection*log(rand);
Travel_Distance_Z_direction = Travel_Distance_magnitude * Angle;
Position = Position + Travel_Distance_Z_direction;


end