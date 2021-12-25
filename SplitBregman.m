function R = SplitBregman(s,beta,gamma)
[rows,cols] = size(s); 

    dx = zeros(rows,cols);
    dy = zeros(rows,cols);
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
  
for iter = 1:2
    R = FFTsolutionofL2norm(dx-bx,dy-by,s,beta/gamma);
x = Dx(R);
y = Dy(R);
[dx,dy] = shrink2isotropic(x+bx,y+by,1/gamma);
  bx = bx-(dx-x);
  by = by-(dy-y);
end
    
       
function d = Dx(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dy(u)
[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return


function [xs,ys] = shrink2isotropic(x,y,lambda) 
s = sqrt(x.^2 + y.^2);
xs = max(0, s-lambda).*x./s;
ys = max(0, s-lambda).*y./s;


function f = FFTsolutionofL2norm(dx,dy,v,lambda)
%the solution of min{1/2int(gradient u-d)^2+lambda/2int(u-v)^2}
sizev = size(v);
otfDx = psf2otf([-1,1],sizev);
    otfDy = psf2otf([-1;1],sizev);
     Nomin = lambda*fft2(v)+otfDx.*fft2(dx)+otfDy.*fft2(dy);
     Denom = abs(otfDx).^2 + abs(otfDy ).^2+lambda;
     s = Nomin./Denom;
     f = real(ifft2(s));
     
     
     
     
     
     
     
     
     
