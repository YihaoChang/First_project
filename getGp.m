function [Gpx,Gpy] = getGp(nx,ny,hx,hy)

Gpx = spdiags(hx,0,nx-1,nx-1)\spdiags(ones(nx-1,1)*[1,-1],[0,1],nx-1,nx);
Gpx = kron(speye(ny),Gpx);
Gpy = spdiags(hy,0,ny-1,ny-1)\spdiags(ones(ny-1,1)*[1,-1],[0,1],ny-1,ny);
Gpy = kron(Gpy,speye(nx));
end