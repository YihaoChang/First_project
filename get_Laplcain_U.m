function A = get_Laplcain_U(nx,ny,hx,hyc)

Dp = spdiags(hx(2:end),0,nx-1,nx-1)\spdiags([-1,0;ones(nx-3,1)*[-1,1];0,1],[0,1],nx-1,nx-1);
Dn = spdiags(hx(1:end-1),0,nx-1,nx-1)\spdiags([-1,1;ones(nx-3,1)*[-1,1];0,1],[-1,0],nx-1,nx-1);

u_Kx = spdiags(hx(2:end)+hx(1:end-1),0,nx-1,nx-1)\ (2*(Dp-Dn));

expand_hyc = [hyc(1);hyc;hyc(end)];

Dp_stagg = spdiags(expand_hyc(2:end),0,ny,ny) \ spdiags([-1,0;ones(ny-2,1)*[-1,1];0,1],[0,1],ny,ny);
Dn_stagg = spdiags(expand_hyc(1:end-1),0,ny,ny) \ spdiags([-1,0;ones(ny-2,1)*[-1,1];0,1],[-1,0],ny,ny);

u_Ky = spdiags(expand_hyc(2:end)+expand_hyc(1:end-1),0,ny,ny)\ (2*(Dp_stagg-Dn_stagg));
A = kron(speye(ny),u_Kx) + kron(u_Ky,speye(nx-1));
end