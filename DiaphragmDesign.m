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
    Fx_rsa(i) = sqrt(sum(fx_rsa_scaled(i,:).^2));
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

%% DIAPHRAGM DESIGN 
fx_elf_sort = sort(Fx, 'descend');
fx_rsa_sort = sort(Fx_rsa, 'descend');

%cumulative sum of elf forces
sum_elf = cumsum(fx_elf_sort);
sum_rsa = zeros(5,1);
sum_rsa(1) = sqrt(sum(fx_rsa_scaled(5,:).^2));

for i = 1:length(Fx)-1
    sum_rsa(i+1) = sqrt(sum(sum(fx_rsa_scaled(5:-1:5-i,:)).^2));
end
sum_wpx = (1:1:5)' * fdl; 

%maximum and minimum fpx 

Fpx_min = 0.2 * Sds * fdl; 
Fpx_max = 0.4 * Sds * fdl; 

for i = 1:length(Fx)
    fpx_elf(i) = (sum_elf(i) / sum_wpx(i)) * fdl; 
    fpx_rsa(i) = (sum_rsa(i) / sum_wpx(i)) * fdl; 
end

%data for elf and rsa for diaphragm design 
elf_data = zeros(5,7);
elf_data(:,1) = (5:-1:1)';
elf_data(:,2) = ones(5,1)*fdl;
elf_data(:,3) = fx_elf_sort;
elf_data(:,4) = sum_elf;
elf_data(:,5) = sum_wpx;
elf_data(:,6) = fpx_elf'; 


rsa_data = zeros(5,7);
rsa_data(:,1) = (5:-1:1)';
rsa_data(:,2) = ones(5,1)*fdl;
rsa_data(:,3) = fx_rsa_sort;
rsa_data(:,4) = sum_rsa;
rsa_data(:,5) = sum_wpx;
rsa_data(:,6) = fpx_rsa'; 


min_fx = ones(6,1)*Fpx_min;
max_fx = ones(6,1)*Fpx_max;

for i = 1:length(Fx)
    if rsa_data(i,3) > rsa_data(i,6) || elf_data(i,3) > elf_data(i,6)
        rsa_data(i,7) = max(rsa_data(i,3), Fpx_min);
        elf_data(i,7) = max(elf_data(i,3), Fpx_min);
    else
        rsa_data(i,7) = max(rsa_data(i,6), Fpx_min);
        elf_data(i,7) = max(elf_data(i,6), Fpx_min);
    end
end


figure
%diaphragm forces elf
x = [5,4,3,2,1,0];

subplot(1,3,1)
plot([elf_data(:,3)',0], x, '--')
hold on 
plot(min_fx, x)
hold on
plot(max_fx, x)
hold on 
plot([elf_data(:, 7)',0], x)
plot([elf_data(:,6)', 0], x, '--')
yticks([0,1,2,3,4,5])
xticks([0,100,200,300,400])
hold off
xlabel('Forces_{ELF}')
ylabel('Story')
legend('Fx_{ELF}', 'Fpx_{min}','Fpx_{max}', 'Design_{ELF} Forces', 'F_{px}', 'Location', 'northwest')
title('Diaphragm Forces ELF')

subplot(1,3,2)
%diaphragm rsa
plot([rsa_data(:,3)',0], x, '--')
hold on 
plot(min_fx, x)
hold on
plot(max_fx, x)
hold on 
plot([rsa_data(:, 7)',0], x)
plot([rsa_data(:,6)', 0], x, '--')
yticks([0,1,2,3,4,5])
xticks([0,100,200,300,400])
hold off
xlabel('Forces_{RSA}')
ylabel('Story')
legend('Fx_{RSA}', 'Fpx_{min}','Fpx_{max}', 'Design_{RSA} Forces', 'F_{px}', 'Location','northwest')
title('Diaphragm Forces RSA')

subplot(1,3,3)
%diaphragm rsa
plot([rsa_data(:,7)',0], x)
hold on 
plot([elf_data(:,7)',0],x)
yticks([0,1,2,3,4,5])
xticks([0,100,200,300,400])
hold off
xlabel('Diaphragm Forces')
ylabel('Story')
legend( 'Design_{RSA} Forces', 'Design_{ELF} Forces', 'Location','northwest')
title('Diaphragm Forces ELF vs RSA')