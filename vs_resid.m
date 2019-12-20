function resid = vs_resid(omega, k, params)
% compute residual of rayleigh eqn at vortex sheet
% following Towne et al. (JFM-2017) paper
    Mj = params.Mj;   % jet Mach number
    gam = params.gam; % ratio of specific heats
    S = params.S;     % density ratio rho_i/rho_o (= To/Ti if ideally expanded)
    m = params.m;     % azimuthal wavenumber
    Rj = params.Rj;   % jet radius

    u_j = 1;            % jet velocity
    a_j = u_j / Mj;     % speed of sound inside jet
    a_o = a_j*sqrt(params.S); % speed of sound outside jet

    gama_o = sqrt(k.^2 - omega.^2);
    gama_i = sqrt(k.^2 - 1./S*(omega - Mj*k).^2);
    if (angle(gama_i) < 0)
        gama_i = -gama_i;
    end;

    Im = besseli(m,gama_i*Rj);
    Imm = besseli(m-1,gama_i*Rj);
    Km = besselk(m,gama_o*Rj);
    Kmm = besselk(m-1,gama_o*Rj);
    resid = 1./(1.-k./omega*Mj).^2  + 1/S*Im./Km*(gama_o*Rj*Kmm + m.*Km)./(gama_i*Rj*Imm - m.*Im);
    %resid = -1*resid;
end
