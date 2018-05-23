function [Dux, Duy] = getDu(nx,ny,hx,hy)

Dux = spdiags(hx,0,nx,nx)\spdiags([ones(nx-2,1)*[1 -1];1,0],[0,-1],nx,nx-1);
Duy = spdiags(hy,0,ny,ny)\spdiags(ones(ny-1,1)*[1 -1],[0,-1],ny,ny-1);

Dux = kron(speye(ny),Dux);
Duy = kron(Duy,speye(nx));

end