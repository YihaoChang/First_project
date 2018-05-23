function A = get_Laplcain_V(nx,ny,hxc,hy)

expand_hxc = [hxc(1);hxc;hxc(end)];


 

Dp_stagg = spdiags(expand_hxc(2:end),0,nx,nx)\spdiags([-1,0;ones(nx-2,1)*[-1,1];0,1],[0,1],nx,nx);
asymmetric_boundary = spdiags([-1,2;ones(nx-2,1)*[-1,1];0,1],[-1,0],nx,nx);
asymmetric_boundary(1,:) = [3,(1/3),zeros(1,nx-2)]; 
Dn_stagg = spdiags(expand_hxc(1:end-1),0,nx,nx)\asymmetric_boundary;

v_Kx = spdiags(expand_hxc(2:end)+expand_hxc(1:end-1),0,nx,nx)\ (2*(Dp_stagg-Dn_stagg));


Dp = spdiags(hy(2:end),0,ny-1,ny-1)\spdiags([-1,0;ones(ny-3,1)*[-1,1];-1,1],[0,1],ny-1,ny-1);
Dn = spdiags(hy(1:end-1),0,ny-1,ny-1)\spdiags([-1,1;ones(ny-3,1)*[-1,1];0,1],[-1,0],ny-1,ny-1);

v_Ky = spdiags(hy(2:end)+hy(1:end-1),0,ny-1,ny-1)\ (2*(Dp-Dn));

A = kron(speye(ny-1),v_Kx) + kron(v_Ky,speye(nx));
end