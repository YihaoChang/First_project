function [U,V] = treat_nonlinesr(U,V,uall,vall,nx,ny,hx,hy,hxc,hyc)

expand_hyc = [hyc(1);hyc;hyc(end)];
dxU1 = spdiags((hx(2:end)+hx(1:end-1)),0,nx-1,nx-1)\(uall(3:end,:) - uall(1:end-2,:));
dyU1 = (spdiags((expand_hyc(2:end)+expand_hyc(1:end-1)),0,ny,ny)\(uall(:,3:end) - uall(:,1:end-2))')';

Vi = (vall(1:end-1,:) + vall(2:end,:)) / 2;
Vi = (Vi(:,1:end-1) + Vi(:,2:end)) / 2;

expand_hxc = [hxc(1);hxc;hxc(end)];
dxU2 = spdiags((expand_hxc(2:end)+expand_hxc(1:end-1)),0,nx,nx)\(vall(3:end,:) - vall(1:end-2,:));
dyU2 = (spdiags((hy(2:end)+hy(1:end-1)),0,ny-1,ny-1)\(vall(:,3:end) - vall(:,1:end-2))')';

Ui = (uall(:,1:end-1) + uall(:,2:end)) / 2;
Ui = (Ui(1:end-1,:) + Ui(2:end,:)) / 2;

U = U.*dxU1(:,2:end-1) + Vi(2:end-1,:).*dyU1(2:end-1,:);
V = Ui(:,2:end-1).*dxU2(:,2:end-1) + V.*dyU2(2:end-1,:);
end