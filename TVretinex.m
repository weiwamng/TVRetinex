
%TV Retinex
%By Wang Wei

clear all

%input
exact = double(imread('test_image.jpg'));
initial = cast(exact, 'uint8');
figure;imshow(initial);
[m, n, c] = size(exact(:, :, 1));

%map into HSV
H0 = rgb2hsv(exact);
H = rgb2hsv(exact);
S0 = H(:, :, 3);
 maxS = max(max(S0));
 minS = min(min(S0));
 V_S = (255/(maxS-minS))*(S0-minS);

%convert into the logarithmic domain
s0 = log(V_S);
s = log(V_S+1);

%initialization
l = s;
z = zeros(m, n);
l_old = 0;

%parameters 
alpha = 1; beta = 0.1; mu = 1e-5; lambda = 1; gamma = 5;

%mean loop
tic;
for iter = 1:100
r = SplitBregman(l-s, beta, lambda);
r = max(r, 0);
l = FFTsolution(z, z, r+s, beta/alpha, mu/alpha);
l = max(l, s);

crit = norm(l-l_old,'fro')/norm(l,'fro');
    if crit < 1e-4
        fiter = iter
        break; 
    end;
    l_old=l;
end
toc;

%gamma correction
L = exp(l);
r0 = s-l;
Ts = log(255)+(1/gamma)*(l-log(255))+r0;
S = exp(Ts);
H(:,:,3) = S;
Final = hsv2rgb(H);

figure;imshow(uint8(Final));







