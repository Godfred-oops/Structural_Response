clear
clc

%mass matrix
d1 = 100; %dead load from slab unit in psf
d2 = 25;  %dead load from superimposed load unit in psf

A = 0.5*80*125;  %tributary area

dl = (d1+d2)/1000; %total dead load unit in ksf

g = 386.2205; %unit in in/s^2

fdl = dl*A;

m = fdl/g;

%stiffness matrix
bw = 18; %unit in inches
hw = 20*12; %unit in inches

Iw = 1/12*(bw*(hw^3));
Ew = 57*sqrt(6000); %unit in ksi
Ec = 57*sqrt(6000); %unit in ksi

h = 10.5*12; %height from 2-5

stiff_coef = 0.5; %stiffness coefficient for cracking in column

h_all = 54; %height of the building

M(1:5,1:5) = diag(m*ones(1,5));  %mass matrix

k1 = stiff_coef *  (12*Ew*Iw)/((1.142871*h)^3); %total stiffness for the first floor (12EI/h^3)
k = stiff_coef *(12*Ew*Iw)/(h^3); %total stiffness for 2-5 floor (12EI/h^3)
k2 = stiff_coef *(6*Ew*Iw)/((1.142871*h)^2); %contribution from the moments first floor (6EI/h^2)
k3 = stiff_coef *(6*Ew*Iw)/(h^2); %contribution from the moments other floors (6EI/h^2)
k4 = stiff_coef *(4*Ew*Iw)/(h); %contribution from the moments other floors (4EI/h)
k5 = stiff_coef *(4*Ew*Iw)/(1.142871*h); %contribution from the moments first floors (4EI/h)
k6 = stiff_coef *(2*Ew*Iw)/(h); %contribution from the moments other floors (2EI/h)


K = [k+k1, -k, 0, 0 ,0, -k2+k3, k3, 0, 0,0;...
    -k, 2*k, -k, 0, 0, -k3, 0,k3, 0, 0;...
    0, -k, 2*k, -k, 0, 0, -k3, 0, k3, 0;...
    0, 0 , -k, 2*k, -k, 0, 0, -k3, 0, k3;...
    0, 0, 0, -k, k,0,0,0,-k3,-k3;...
    -k2+k3,-k3,0,0,0,k5+k4, k6,0,0,0;...
    k3,0,-k3,0,0,k6,2*k4,k6,0,0;...
    0,k3,0,-k3,0,0,k6,2*k4,k6,0;...
    0,0,k3,0,-k3,0,0,k6,2*k4,k6;...
    0,0,0,k3,-k3,0,0,0,k6,k4]; %stiffness matrix

Ktt = K(1:5,1:5); %stiffness matrix for both translation
Kto = K(1:5,6:10); %stiffness matrix for translation and rotation
Kot = K(6:10,1:5); %stiffness matrix for rotation and translation
Koo = K(6:10,6:10); %stiffness matrix for both rotation

ktr = Ktt - Kto*inv(Koo)*Kot;
Mtr = M(1:5,1:5);

[V,D] = eig(ktr,Mtr);

for i = 1:size(V, 2)
    subplot(1, length(V), i);
    X = (0:1:length(V));
    Y = [0,V(:,i)'];
    plot( Y,X, 'o-')
    % Set custom y-axis tick locations
    ytick_values = [0, 1, 2, 3, 4, 5, 6];
    yticks(ytick_values);
    % Remove x-axis ticks
    xticks([]);
    title(["Mode", num2str(6-i)])
end

%damping matrix
w1 = (diag(D)'); %w^2 (omega square) for the various modes
w = zeros(1,5);
for i = 1:length(w1)
    w(i) = w1(6-i);
end
Period = zeros(1,5);

%period
for j = 1:length(w)
    Period(1,j) = (2*pi)/sqrt(w(j));
end

%% Equivalent Lateral Force Procedure (ELF)

A1 = 0.125*80*125;

Wx = [A1,A1,A1,A1,A1];
n_story = 5;
h1 = 12;
h_all = 10.5;
hx = zeros(n_story,1);
for i = 1:n_story
    if i == 1
        hx(i) = 12;
    else
        hx(i) = hx(i-1) + h_all;
    end
end

% determine exponent k
if Period(1) < 0.5
    k = 1;
else
    k = 2;
end
% parameters for shear wall
Ct = 0.02;
x = 0.75;
Cu = 1.4;
Sds = 1.518;
Sd1 = 0.615;
S1 = 0.659;
Cd = 5;
R = 5;
Ie = 1;
T = Period(1);
Ta = max(0.1*n_story, Ct * sum(hx)^x);

% find Cs
Cs1 = min(Sds/R , Sd1/(T*(R/Ie)));
Cs2 = Sds*0.044*Ie;
Cs3 = 0.5*S1/R;
Cs = max(max(Cs1,Cs2),Cs3);

% solve for forces
hxk = hx.^k;
Wxhxk = Wx' .* hxk;
Cvx = Wxhxk / sum(Wxhxk);
V1 = Cs * sum(Wx);
Fx = (Cvx * V1)/2; % divide forces by two to get forces on each shear wall

% base shear (kips)
V_elf = sum(Fx);
% base moment (kip-in)
M_elf = 0;
for i = 1:length(hx)
    M_elf = M_elf + Fx(i) * hx(i) * 12;
end
% lateral floor displacements (in)
x_elf = inv(ktr) * Fx * Cd;

%% Response Spectrum Analysis (RSA)
% mode shapes and participation factors
phi = zeros(length(V),length(V));
% reverse mode order of phi
for i = 1:length(phi)
    phi(:,i) = V(:,6-i);
end
% mass normalize phi
alpha = zeros(length(phi),1);
for i = 1:length(phi)
    alpha(i) = sqrt(inv(phi(:,i)'*M*(phi(:,i))));
end 
phi_norm = phi;
for i = 1:length(phi)
    phi_norm(:,i) = alpha(i)*phi(:,i);
end 

% solve for gamma
r = zeros(n_story,1) + 1;
gamma_n = zeros(n_story,1);
for i = 1:length(phi)
    gamma_n(i) = (phi_norm(:,i)' * M * r) / (phi_norm(:,i)'* M *phi_norm(:,i));
end

% solve for Modal Mass Participation (MMP)
MMP = zeros(length(phi),1);
for i = 1:length(MMP)
    MMP(i) = gamma_n(i)^2 /(5*m);
end

% solve for static story forces
s_n = zeros(length(phi),length(phi));
for i = 1:length(phi)
    s_n(:,i) = gamma_n(i) * M * phi_norm(:,i);
end

% find pseudo-accelerations
Ts = Sd1/Sds;
T0 = 0.2*Sd1/Sds;
% two modes required 90% mass participation per ASCE 7-16, but use all
PSa_n = zeros(length(phi),1); % in/s^2
for i = 1:length(PSa_n)
    if Period(i) < T0
        PSa_n(i) = 0.6 * (Sds/T0) * Period(i) + 0.4*Sds;
    elseif Period(i) < Ts
        PSa_n(i) = Sds;
    else
        PSa_n(i) = Sd1 /Period(i);
    end
end

% calculate dyanmic story forces
fx_rsa = zeros(length(phi),length(phi));
for i = 1:length(phi)
    fx_rsa(:,i) = PSa_n(i) * g * s_n(:,i)/R;
end

% calculate base shear for each mode
V_rsa_modal = zeros(length(phi),1);
for i = 1:length(phi)
    V_rsa_modal(i) = sum(fx_rsa(:,i));
end
% compute SRSS combination for base shear and compute scaled V_rsa
V_rsa = sqrt(sum(V_rsa_modal.^2));
scale_factor = V_elf/V_rsa;
V_rsa_scaled = V_rsa * scale_factor;
fx_rsa_scaled = fx_rsa * scale_factor;

% compute scaled SRSS combination for base moment
M_rsa_modal = zeros(length(phi),1);
for i = 1:length(phi)
    for j = 1:length(hx)
        M_rsa_modal(i) = M_rsa_modal(i) + fx_rsa_scaled(j,i) * hx(j) * 12;
    end
end
M_rsa_scaled = sqrt(sum(M_rsa_modal.^2));

% compute SRSS combination for dynamic story forces
Fx_rsa = zeros(length(phi),1);
for i = 1:length(phi)
    Fx_rsa(i) = sqrt(sum(fx_rsa(i,:).^2));
end

% calculate modal story displacements and SRSS combination
x_rsa_modal = zeros(length(phi),length(phi));
for i = 1:length(phi)
    x_rsa_modal(:,i) = inv(ktr)*fx_rsa(:,i) * Cd;
end
x_rsa = zeros(length(phi),1);
for i = 1:length(phi)
    x_rsa(i) = sqrt(sum(x_rsa_modal(i,:).^2));
end

%% Response History Analysis (RHA)
% same as RSA but given three ground motions, 2 directions => 6 ground
% motions





%% Response History Analysis (RHA)
% same as RSA but given three ground motions, 2 directions => 6 ground
% motions

p = dir("*.txt");
%number of data points
l = zeros(1, numel(p));
for j = 1:numel(p)
    l(1,j) = (length((table2array(readtable(p(j).name))))-1)*5;
end
data = zeros(numel(p),max(l)); %storing the data from the txt files
dt2 = zeros(1,numel(p)); %matrix of the time steps

%READING TIME STEPS

for i = 1:numel(p)
    %reading the text file
    fileID = fopen(p(i).name, 'r');
    d = readtable(p(i).name, 'HeaderLines',3);

    %time step
    dt = d(1,1:5);
    dt1 = table2array(dt); %convert the table to array
    dt2(1,i) = dt1(:,4);   
       
end

%READING DATA POINTS
for i = 1:numel(p)
    %reading the text file
    fileID = fopen(p(i).name, 'r');
    d = readtable(p(i).name, 'HeaderLines',3);

    %data points
    pt = d(2:end, 1:5);
    pt1 = table2array(pt);

    n = 1;
    for j = 1:length(pt1)
    data(i,n:n+4) = pt1(j,:);
    n = n+5;
    end

end

%number of data points
npt = zeros(1, numel(p));
for i = 1:size(data,1)
    np = data(i,:);
    np1 = np(~isnan(np));
    npt(1,i) = length(np1(np1 ~= 0));
end
%reading each ground motion

kocaeli_00 = data(1,1:npt(1));
kocaeli_90 = data(2,1:npt(2));
niigata_ew = data(3,1:npt(3));
niigata_ns = data(4,1:npt(4));
northr_00 = data(5,1:npt(5));
northr_90 = data(6,1:npt(6));

%modal mass
mn = zeros(1, n_story);
for i = 1:n_story
    mn(1, i) = sum(s_n(:,i));
end

% w2 = sort(w1);
% ngamma = 1/2; beta=1/4; zeta = 0.05;
% D_kocaeli = zeros(n_story, length(kocaeli_00));
% A_kocaeli = zeros(n_story, length(kocaeli_00));
% 
% 
% response = data(1,1:npt(1))';
% for i = 1:n_story
%         m= mn(i);
%         P = -g*m*response;
%         k = w2(i)*m;
%         w = sqrt(w2(i));
%         c = 2*m*w*zeta;
%         [a, v, u] = NewmarkIntegrator(ngamma, beta, m,  c,k, P, dt2(1));
%         D_kocaeli(i,:) = u;
%         A_kocaeli(i,:) = w^2*D_kocaeli(i,:);
% end
% 


% modal response 
% acceleration and displacement for the first gm
w2 = sort(w1);
ngamma = 1/2; beta=1/4; zeta = 0.05;
D_kocaeli = zeros(2*n_story, length(kocaeli_00));
A_kocaeli = zeros(2*n_story, length(kocaeli_00));
h = 0;
for j = 1:2
    response = data(j,1:npt(j))';
    for i = 1:n_story
        m= 1;
        P = -g*m*response;
        k = w2(i)*m;
        w = sqrt(w2(i));
        c = 2*m*w*zeta;
        [a, v, u] = NewmarkIntegrator(ngamma, beta, m,  c,k, P, dt2(j));
        D_kocaeli(h+i,:) = u;
        A_kocaeli(h+i,:) = w^2*D_kocaeli(h+i,:);
    end
    h = h+5; 

end

%acceleration and displacement for the second gm
D_niigata = zeros(2*n_story, length(niigata_ew));
A_niigata = zeros(2*n_story, length(niigata_ns));
h = 0;
for j = 1:2
    response = data(j+2,1:npt(j+2))';
    for i = 1:n_story
        m= 1;
        P = -g*m*response;
        k = w2(i)*m;
        w = sqrt(w2(i));
        c = 2*m*w*zeta;
        [a, v, u] = NewmarkIntegrator(ngamma, beta, m,  c,k, P, dt2(j+2));
        D_niigata(h+i,:) = u;
        A_niigata(h+i,:) = w^2*D_niigata(h+i,:);
    end
    h = h+5; 

end

%acceleration and displacement for the third gm
D_northr = zeros(2*n_story, length(northr_00));
A_northr = zeros(2*n_story, length(northr_90));
h = 0;

for j = 1:2
    response = data(j+4,1:npt(j+4))';
    for i = 1:n_story
        m= 1;
        P = -g*m*response;
        k = w2(i)*m;
        w = sqrt(w2(i));
        c = 2*m*w*zeta;
        [a, v, u] = NewmarkIntegrator(ngamma, beta, m,  c,k, P, dt2(j+4));
        D_northr(h+i,:) = u;
        A_northr(h+i,:) = w^2*D_northr(h+i,:);
    end
    h = h+5; 

end

%total response for acceleration and displacement in both directions for
%GM1
A1_kocaeli = zeros(2*n_story, length(kocaeli_90));
D1_kocaeli = zeros(2*n_story, length(kocaeli_90));
d = 1;
for i = 1:2
    A1_kocaeli(d:d+4, :) = phi*diag(gamma_n)*A_kocaeli(d:d+4, 1:length(kocaeli_90));
    D1_kocaeli(d:d+4, :) = phi*diag(gamma_n)*D_kocaeli(d:d+4, 1:length(kocaeli_90));
    d = d+5;
end

%total response for acceleration and displacement in both directions for
%GM2
A1_niigata = zeros(2*n_story, length(niigata_ew));
D1_niigata = zeros(2*n_story, length(niigata_ns));
d = 1;
for i = 1:2
    A1_niigata(d:d+4, :) = phi*diag(gamma_n)*A_niigata(d:d+4, 1:length(niigata_ew));
    D1_niigata(d:d+4, :) = phi*diag(gamma_n)*D_niigata(d:d+4, 1:length(niigata_ns));
    d = d+5;
end

%total response for acceleration and displacement in both directions for
%GM3
A1_northr = zeros(2*n_story, length(northr_00));
D1_northr = zeros(2*n_story, length(northr_00));
d = 1;
for i = 1:2
    A1_northr(d:d+4, :) = phi*diag(gamma_n)*A_northr(d:d+4, 1:length(northr_00));
    D1_northr(d:d+4, :) = phi*diag(gamma_n)*D_northr(d:d+4, 1:length(northr_00));
    d = d+5;
end

pseudo_acc_kocaeli = zeros(n_story, 1);
dis_acc_kocaeli = zeros(n_story, 1);
%pseudo_acceleration
for i = 1:n_story 
    pseudo_acc_kocaeli(i,1) = max(max(A1_kocaeli(i,:)), max(A1_kocaeli(i+5,:)));
    dis_acc_kocaeli(i,1) = max(max(D1_kocaeli(i,:)), max(D1_kocaeli(i+5,:)));
end

pseudo_acc_niigata = zeros(n_story, 1);
 dis_acc_niigata = zeros(n_story, 1);
%pseudo_acceleration
for i = 1:n_story 
    pseudo_acc_niigata(i,1) = max(max(A1_niigata(i,:)), max(A1_niigata(i+5,:)));
    dis_acc_niigata(i,1) = max(max(D1_niigata(i,:)), max(D1_niigata(i+5,:)));

end

pseudo_acc_northr = zeros(n_story, 1);
dis_acc_northr = zeros(n_story, 1);
%pseudo_acceleration
for i = 1:n_story
    pseudo_acc_northr(i,1) = max(max(A1_northr(i,:)), max(A1_northr(i+5,:)));
    dis_acc_northr(i,1) = max(max(D1_northr(i,:)), max(D1_northr(i+5,:)));
end

% %Dn calculation
% D_n = zeros(length(gamma_n),1);
% for i = 1:length(gamma_n)
%     D_n(i,1) = dis_acc_niigata(i)/gamma_n(i);
% end
% 
% rr = ones(5,1);
% v_base = zeros(length(gamma_n),1);
% for i = 1:length(gamma_n)
%     v_base(i,1) = phi_norm(:,i)'*Mtr*rr*gamma_n(i)*w2(i)*D_n(i);
% end
%base shear

% calculate dyanmic story forces rha calculation
fx_rha = zeros(length(phi),length(phi));
for i = 1:length(phi)
    fx_rha(:,i) = pseudo_acc_niigata(i) * s_n(:,i)/R;
end

% calculate base shear for each mode rha
V_rha_modal = zeros(length(phi),1);
for i = 1:length(phi)
    V_rha_modal(i) = sum(fx_rha(:,i));
end

V_rha = sqrt(sum(V_rha_modal.^2));
scale_factor = V_elf/V_rha;
V_rha_scaled = V_rha * scale_factor;
fx_rha_scaled = fx_rha * scale_factor;

% compute SRSS combination for dynamic story forces rha
Fx_rha = zeros(length(phi),1);
for i = 1:length(phi)
    Fx_rha(i) = sqrt(sum(fx_rha(i,:).^2));
end

% calculate modal story displacements and SRSS combination rha
x_rha_modal = zeros(length(phi),length(phi));
for i = 1:length(phi)
    x_rha_modal(:,i) = inv(ktr)*fx_rha(:,i) * Cd;
end
x_rha = zeros(length(phi),1);
for i = 1:length(phi)
    x_rha(i) = sqrt(sum(x_rha_modal(i,:).^2));
end

%modal height
%hn = (phi_norm' * M * hx )/(phi_norm' * M * r);

% compute scaled SRSS combination for base moment
M_rha_modal = zeros(length(phi),1);
for i = 1:length(phi)
    for j = 1:length(hx)
        M_rha_modal(i) = M_rha_modal(i) + fx_rha_scaled(j,i) * hx(j) * 12;
    end
end
M_rha_scaled = sqrt(sum(M_rha_modal.^2));



%% Compare values
fprintf('ELF vs RSA vs RHA Results\nLateral floor level displacement for Stories 1-5 (in)\n')
displacement = [x_elf, x_rsa, x_rha];
disp(displacement)
fprintf('Base moment (kip-in)\n')
moment = [M_elf, M_rsa_scaled,M_rha_scaled];
disp(moment)
fprintf('Base Shear (kips)\n')
base_shear = [V_elf, V_rsa_scaled, V_rha_scaled];
disp(base_shear)


% figure; 
% xx = 0:dt2(1):(npt(1)-1)*dt2(1);
% yy = A(:,2);
% plot(xx,yy)