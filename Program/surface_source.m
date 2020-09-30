function psi = surface_source(Omega_z, Source_Type)

switch Source_Type
    case 'Isotropic' 
        if Omega_z > 0
            psi = 1;
        else
            warning('Unexpected Omega_z as Isotropic surface sources')
        end
    case 'Linear'
        psi = 1 + Omega_z;
    case 'Discrete'
        psi = dirac(Omega_z - 1);
    otherwise
        warning('Unexpected surface source type.')
end

end

%% the boundary condition at z = 0