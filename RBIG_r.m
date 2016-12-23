function [ans] = RBIG_r(dat, N_lay, transformation, porc, precision)
[MI, MIs] = MI_RBIG_2016(dat,N_lay,transformation,porc,precision);
% Save results
ans = struct('MI', MI);

