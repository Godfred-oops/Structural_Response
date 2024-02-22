clear
clc

%mass matrix
d1 = 100; %dead load from slab unit in psf
d2 = 25;  %dead load from superimposed load unit in psf

A = 0.25* 80*125;  %tributary area

dl = (d1+d2)/1000; %total dead load unit in ksf

g = 386.2205; %unit in in/s^2

fdl = dl*A;

m = fdl/g;

M = zeros(15,15);
M(1:5,1:5) = diag(m*ones(1,5));  %mass matrix

%stiffness matrix
bc = 20; %unit in inches
hc = 30; %unit in inches
bb = 20; %unit in inches
hb = 32; %unit in inches


Ic = 1/12*(bc*((hc)^3));
Ec = 57*sqrt(6000); %unit in ksi
Ib = 105.207e3;
%Ib = 1/12*(bb*(hb^3)); 
Eb = 57*sqrt(6000);


h = 10.5*12; %height from 2-4 in inches
Li = 25*12; %length for the beam in the moment-frame
stiff_coef = 0.5; %stiffness coefficient for cracking in column

k1 = stiff_coef * 2*(12*Ec*Ic)/((1.142871*h)^3); %total stiffness for the first floor (24EI/1.25h^3)
k = stiff_coef *(24*Ec*Ic)/(h^3); %total stiffness for 2-5 floor (24EI/h^3)
k2 = stiff_coef *(6*Ec*Ic)/(h^2); %total stiffness for 2-5 floor (6EI/h^2)
k3 = stiff_coef *(6*Ec*Ic)/((1.142871*h)^2); %total stiffness for first floor (6EI/1.25h^2)
k4 = stiff_coef *(4*Ec*Ic)/(1.142871*h); %total stiffness for first floor (4EI/1.25h)
k5 = stiff_coef *(4*Ec*Ic)/(h); %total stiffness for 2-5 floor (4EI/h)
k6 = stiff_coef *(4*Eb*Ib)/(Li); %total stiffness for moment first floor (4EI/L)
k7 = stiff_coef *(2*Eb*Ib)/(Li); %total stiffness for first floor (2EI/L)
k8 = stiff_coef *(2*Ec*Ic)/(h); %total stiffness for first floor (2EI/h)



K = [k+k1, -k, 0, 0 ,0, -k3+k2, k2, 0,0,0,-k3+k2, k2, 0,0,0;...
    -k, 2*k, -k, 0, 0, -k2, 0, k2, 0,0, -k2, 0, k2, 0,0;...
    0, -k, 2*k, -k, 0, 0, -k2, 0, k2, 0, 0, -k2, 0, k2, 0;...
    0, 0 , -k, 2*k, -k, 0,0, -k2, 0, k2, 0, 0, -k2, 0, k2;...
    0, 0, 0, -k, k, 0,0,0, -k2, -k2, 0,0,0, -k2,-k2;...
    -k3+k2, -k2, 0, 0, 0, k4+k5+k6, k8, 0,0,0,k7,0,0,0,0;...
    k2, 0, -k2, 0, 0, k8, 2*k5 + k6, k8, 0,0,0,k7,0,0,0;...
    0, k2, 0, -k2, 0, 0, k8, 2*k5+k6, k8, 0,0,0,k7,0,0;...
    0, 0, k2, 0, -k2, 0 , 0, k8, 2*k5+k6, k8, 0,0,0,k7,0;...
    0, 0, 0, k2, -k2, 0, 0, 0,k8, k5, 0, 0,0,0,k7;...
    -k3+k2, -k2, 0,0,0, k7, 0, 0, 0, 0, 2*k5+k6, k8, 0,0,0;...
    k2, 0, -k2, 0, 0, 0, k7, 0, 0, 0, k8, 2*k5+k6, k8, 0, 0;...
    0, k2, 0, -k2, 0, 0, 0, k7, 0, 0, 0, k8, 2*k5+k6, k8, 0;...
    0, 0, k2, 0, -k2, 0, 0,0,0,k7,0,0,k8, 2*k5+k6, k8;...
    0,0,0,k2,-k2,0,0,0,0,k7,0,0,0,k8,k5+k6];


Ktt = K(1:5,1:5); %stiffness matrix for both translation
Kto = K(1:5,6:15); %stiffness matrix for translation and rotation
Kot = K(6:15,1:5); %stiffness matrix for rotation and translation
Koo = K(6:15,6:15); %stiffness matrix for both rotation

ktr = Ktt - Kto*inv(Koo)*Kot;
Mtr = M(1:5,1:5);
[V,D] = eig(ktr,Mtr);


%mode shapes
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

sgtitle('Mode Shapes for moment frame');
%damping matrix
w1 = (diag(D)'); %w^2 (omega square) for the various modes
w = zeros(1,5);
for i = 1:length(w1)
    w(i) = w1(6-i);
end
dr1 = 0.05; %damping ratio for first mode
dr4 = 0.05; %damping ratio for fourth mode

b = (2*dr1*sqrt(w(1))-2*dr4*sqrt(w(4)))/(w(1)-w(4));
a = 2*dr1*sqrt(w(1))-b*w(1);

Dm = a*M + b*K;

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
Ct = 0.016;
x = 0.9;
Cu = 1.4;
Sds = 1.518;
Sd1 = 0.615;
S1 = 0.659;
Cd = 5.5;
R = 8;
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
Fx = (Cvx * V1)/2; % divide forces by two to get forces on each level of the moment frame
% base shear (kips)
V_elf = sum(Fx);
% % base moment (kip-in)
% M_elf = 0;
% for i = 1:length(hx)
%     M_elf = M_elf + Fx(i) * hx(i) * 12;
% end
% lateral floor displacements (in)
x_elf = inv(ktr) * Fx * Cd;

%% member forces for the ELF

%displacement for the drift calculations
dx = ktr\Fx;

fx1 = [Fx;0;0;0;0;0;0;0;0;0;0];
dx1 = K\fx1;

nele = 10;
L = zeros(nele, 4);
for i = 1:5
    if i == 1
        L(i,1) = 0;
        L(i,2) = 1;
        L(i,3) = 0;
        L(i,4) = 1 + nele/2;
    else
        L(i,1) = i - 1;
        L(i,2) = i;
        L(i,3) = nele/2 + i - 1;
        L(i,4) = nele/2 + i; 
    end
end
for j = 6:10
    L(j,1) = L(j-5,1);
    L(j,2) = L(j-5,2);
    if j == 6
        L(j,3) = 0;
    else
        L(j,3) = L(j-5,3) + 5;
    end
    L(j,4) = L(j-5,4) + 5;
end

length = [12,10.5,10.5,10.5,10.5,12,10.5,10.5,10.5,10.5];
column_member = zeros(10,4);

lm = 12*[12,10.5,10.5,10.5,10.5,12,10.5,10.5,10.5,10.5];


for j = 1:10
    kc = (Ec*0.5*Ic)/(lm(j)^3)*[12, -12, 6*lm(j), 6*lm(j);...
    -12, 12, -6*lm(j), -6*lm(j);...
    6*lm(j), -6*lm(j), 4*(lm(j))^2, 2*(lm(j))^2;...
    6*lm(j), -6*lm(j), 2*(lm(j))^2, 4*(lm(j))^2];
    if j == 1 || j == 6
        Lj = L(j,:); % extract member locator vector
        column_member(j,:) = (kc * [0; dx1(Lj(2)); 0; dx1(Lj(4))]);  
    else 
        column_member(j,:) = (kc * dx1(L(j,:)));

    end
end


beam_member = zeros(5,2);

%beam demands/moment
Lb = L(:,4);
beam_coef = stiff_coef*(Eb*Ib)/(Li^3);
kbb = [4*(Li^2), 2*(Li^2);2*(Li^2), 4*(Li^2)];
kb = beam_coef*kbb;
for i = 1:5
    beam_member(i,:) = kb*dx1([Lb(i), Lb(i+5)]);
 
end

shear_demand = zeros(5,1);
for k = 1:size(beam_member,1)
    shear_demand(k,:) = (beam_member(k,1) + beam_member(k,2))/Li;

end
%% Response Spectrum Analysis (RSA)
% mode shapes and participation factors
phi = zeros(size(V,2),size(V,2));
% reverse mode order of phi
for i = 1:size(phi,2)
    phi(:,i) = V(:,6-i);
end
% mass normalize phi
alpha = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    alpha(i) = sqrt(inv(phi(:,i)'*Mtr*(phi(:,i))));
end 
phi_norm = phi;
for i = 1:size(phi,2)
    phi_norm(:,i) = alpha(i)*phi(:,i);
end 

% solve for gamma
r = zeros(n_story,1) + 1;
gamma_n = zeros(n_story,1);
for i = 1:size(phi,2)
    gamma_n(i) = (phi_norm(:,i)' * Mtr * r) / (phi_norm(:,i)'* Mtr *phi_norm(:,i));
end

% solve for Modal Mass Participation (MMP)
MMP = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    MMP(i) = gamma_n(i)^2 /(5*m);
end

% solve for static story forces
s_n = zeros(size(phi,2),size(phi,2));
for i = 1:size(phi,2)
    s_n(:,i) = gamma_n(i) * Mtr * phi_norm(:,i);
end

% find pseudo-accelerations
Ts = Sd1/Sds;
T0 = 0.2*Sd1/Sds;
% two modes required 90% mass participation per ASCE 7-16, but use all
PSa_n = zeros(size(phi,2),1); % in/s^2
for i = 1:size(phi,2)
    if Period(i) < T0
        PSa_n(i) = 0.6 * (Sds/T0) * Period(i) + 0.4*Sds;
    elseif Period(i) < Ts
        PSa_n(i) = Sds;
    else
        PSa_n(i) = Sd1 /Period(i);
    end
end

% calculate dynamic story forces
fx_rsa = zeros(size(phi,2),size(phi,2));
for i = 1:size(phi,2)
    fx_rsa(:,i) = PSa_n(i) * g * s_n(:,i)/R;
end

% calculate base shear for each mode
V_rsa_modal = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    V_rsa_modal(i) = sum(fx_rsa(:,i));
end
% compute SRSS combination for base shear and compute scaled V_rsa
V_rsa = sqrt(sum(V_rsa_modal.^2));
scale_factor = V_elf/V_rsa;
V_rsa_scaled = V_rsa * scale_factor;
fx_rsa_scaled = fx_rsa * scale_factor;


% compute SRSS combination for dynamic story forces
Fx_rsa = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    Fx_rsa(i) = sqrt(sum(fx_rsa(i,:).^2));
end

% calculate modal story displacements and SRSS combination
x_rsa_modal = zeros(size(phi,2),size(phi,2));
for i = 1:size(phi,2)
    x_rsa_modal(:,i) = inv(ktr)*fx_rsa(:,i) * Cd;
end
x_rsa = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    x_rsa(i) = sqrt(sum(x_rsa_modal(i,:).^2));
end

%% BEAM AND COLUMN DEMANDS FOR RSA


dx2 = zeros(size(phi,2)+10, size(phi,2));
for k = 1:size(phi,2)
    fx1_rsa = [fx_rsa(:,k);0;0;0;0;0;0;0;0;0;0];
    dx2(:,k) = K\fx1_rsa;
end

%using just the first two modes
%first mode
dx3 = dx2(:,1);
for j = 1:10

    kc = (Ec*0.5*Ic)/(lm(j)^3)*[12, -12, 6*lm(j), 6*lm(j);...
        -12, 12, -6*lm(j), -6*lm(j);...
        6*lm(j), -6*lm(j), 4*(lm(j))^2, 2*(lm(j))^2;...
        6*lm(j), -6*lm(j), 2*(lm(j))^2, 4*(lm(j))^2];

    if j == 1 || j == 6
        Lj = L(j,:); % extract member locator vector
            column_member_rsa1(j,:) = (kc * [0; dx3(Lj(2)); 0; dx3(Lj(4))]);  
    else 
            column_member_rsa1(j,:) = (kc * dx3(L(j,:)));
    
    end
end


dx4 = dx2(:,2);
for j = 1:10

    kc = (Ec*0.5*Ic)/(lm(j)^3)*[12, -12, 6*lm(j), 6*lm(j);...
        -12, 12, -6*lm(j), -6*lm(j);...
        6*lm(j), -6*lm(j), 4*(lm(j))^2, 2*(lm(j))^2;...
        6*lm(j), -6*lm(j), 2*(lm(j))^2, 4*(lm(j))^2];

    if j == 1 || j == 6
        Lj = L(j,:); % extract member locator vector
            column_member_rsa2(j,:) = (kc * [0; dx4(Lj(2)); 0; dx4(Lj(4))]);  
    else 
            column_member_rsa2(j,:) = (kc * dx4(L(j,:)));
    
    end
end



%combining the two modes with SRSS
for i = 1:4
    mode2_member = column_member_rsa2(:,i);
    mode1_member = column_member_rsa1(:,i);
    for j = 1:10
        combined_rsa_member(j,i) = -scale_factor * sqrt(mode1_member(j)^2 + mode2_member(j)^2);

    end
end


%beam demands/moment for the first mode
beam_member1 = zeros(5,2);
beam_member2 = zeros(5,2);
Lb = L(:,4);
beam_coef = stiff_coef*(Eb*Ib)/(Li^3);
kbb = [4*(Li^2), 2*(Li^2);2*(Li^2), 4*(Li^2)];
kb = beam_coef*kbb;
for i = 1:5
    beam_member1(i,:) =  kb*dx3([Lb(i), Lb(i+5)]);
    beam_member2(i,:) =  kb*dx4([Lb(i), Lb(i+5)]);

end

%combining the two modes
for j = 1:2
    mode1_beam = beam_member1(:,j);
    mode2_beam = beam_member2(:,j);
    for i = 1:5
        beam_member_combined(i,j) = scale_factor * sqrt(mode1_beam(i)^2 + mode2_beam(i)^2);
    end
    
end
% shear demand for RSA
shear_demand1 = zeros(5,1);
for k = 1:size(beam_member,1)
    shear_demand1(k,:) = (beam_member_combined(k,1) + beam_member_combined(k,2))/Li;

end

%% RESPONSE HISTORY ANALYSIS
p = dir("*.txt");
%number of data points
l = zeros(1, numel(p));
for j = 1:numel(p)
    l(1,j) = (size((table2array(readtable(p(j).name))),1)-1)*5;
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
    for j = 1:size((pt1),1)
    data(i,n:n+4) = pt1(j,:);
    n = n+5;
    end

end

%number of data points
npt = zeros(1, numel(p));
for i = 1:size(data,1)
    np = data(i,:);
    np1 = np(~isnan(np));
    npt(1,i) = size(np1(np1 ~= 0),2);
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

% acceleration and displacement for the first gm
w2 = sort(w1);
ngamma = 1/2; beta=1/4; zeta = 0.05;
D_kocaeli = zeros(2*n_story, size(kocaeli_00,2));
A_kocaeli = zeros(2*n_story, size(kocaeli_90,2));
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
D_niigata = zeros(2*n_story, size(niigata_ew,2));
A_niigata = zeros(2*n_story, size(niigata_ns,2));
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
D_northr = zeros(2*n_story, size(northr_00,2));
A_northr = zeros(2*n_story, size(northr_90,2));
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
A1_kocaeli = zeros(2*n_story, size(kocaeli_90,2));
D1_kocaeli = zeros(2*n_story, size(kocaeli_90,2));
d = 1;
for i = 1:2
    A1_kocaeli(d:d+4, :) = phi*diag(gamma_n)*A_kocaeli(d:d+4, 1:size(kocaeli_90,2));
    D1_kocaeli(d:d+4, :) = phi*diag(gamma_n)*D_kocaeli(d:d+4, 1:size(kocaeli_90,2));
    d = d+5;
end

%total response for acceleration and displacement in both directions for
%GM2
A1_niigata = zeros(2*n_story, size(niigata_ew,2));
D1_niigata = zeros(2*n_story, size(niigata_ns,2));
d = 1;
for i = 1:2
    A1_niigata(d:d+4, :) = phi*diag(gamma_n)*A_niigata(d:d+4, 1:size(niigata_ew,2));
    D1_niigata(d:d+4, :) = phi*diag(gamma_n)*D_niigata(d:d+4, 1:size(niigata_ns,2));
    d = d+5;
end

%total response for acceleration and displacement in both directions for
%GM3
A1_northr = zeros(2*n_story, size(northr_00,2));
D1_northr = zeros(2*n_story, size(northr_00,2));
d = 1;
for i = 1:2
    A1_northr(d:d+4, :) = phi*diag(gamma_n)*A_northr(d:d+4, 1:size(northr_00,2));
    D1_northr(d:d+4, :) = phi*diag(gamma_n)*D_northr(d:d+4, 1:size(northr_00,2));
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

% calculate dyanmic story forces rha calculation
fx_rha = zeros(size(phi,2),size(phi,2));
for i = 1:size(phi,2)
    fx_rha(:,i) = pseudo_acc_niigata(i) * s_n(:,i)/R;
end

% calculate base shear for each mode rha
V_rha_modal = zeros(size(phi,2));
for i = 1:size(phi,2)
    V_rha_modal(i) = sum(fx_rha(:,i));
end

%scale factor
V_rha = sqrt(sum(V_rha_modal.^2));
scale_factor1 = V_elf/sum(V_rha);


% calculate modal story displacements and SRSS combination rha
x_rha_modal = zeros(size(phi,2),size(phi,2));
for i = 1:size(phi,2)
    x_rha_modal(:,i) = inv(ktr)*fx_rha(:,i) * Cd;
end

x_rha = zeros(size(phi,2),1);
for i = 1:size(phi,2)
    x_rha(i) = sqrt(sum(x_rha_modal(i,:).^2));
end
%% BEAM AND COLUMN DEMAND FOR RHA

dx5 = zeros(size(phi,2)+10, size(phi,2));
for k = 1:size(phi,2)
    fx1_rha = [fx_rha(:,k);0;0;0;0;0;0;0;0;0;0];
    dx5(:,k) = K\fx1_rha;
end

%using just the first two modes
%first mode
dx6 = dx5(:,1);
for j = 1:10

    kc = (Ec*0.5*Ic)/(lm(j)^3)*[12, -12, 6*lm(j), 6*lm(j);...
        -12, 12, -6*lm(j), -6*lm(j);...
        6*lm(j), -6*lm(j), 4*(lm(j))^2, 2*(lm(j))^2;...
        6*lm(j), -6*lm(j), 2*(lm(j))^2, 4*(lm(j))^2];

    if j == 1 || j == 6
        Lj = L(j,:); % extract member locator vector
            column_member_rha1(j,:) = (kc * [0; dx6(Lj(2)); 0; dx6(Lj(4))]);  
    else 
            column_member_rha1(j,:) = (kc * dx6(L(j,:)));
    
    end
end


dx7 = dx5(:,2);
for j = 1:10

    kc = (Ec*0.5*Ic)/(lm(j)^3)*[12, -12, 6*lm(j), 6*lm(j);...
        -12, 12, -6*lm(j), -6*lm(j);...
        6*lm(j), -6*lm(j), 4*(lm(j))^2, 2*(lm(j))^2;...
        6*lm(j), -6*lm(j), 2*(lm(j))^2, 4*(lm(j))^2];

    if j == 1 || j == 6
        Lj = L(j,:); % extract member locator vector
            column_member_rha2(j,:) = (kc * [0; dx7(Lj(2)); 0; dx7(Lj(4))]);  
    else 
            column_member_rha2(j,:) = (kc * dx7(L(j,:)));
    
    end
end


%combining the two modes with SRSS
for i = 1:4
    mode2_member_rha = column_member_rha2(:,i);
    mode1_member_rha = column_member_rha1(:,i);
    for j = 1:10
        combined_rha_member(j,i) = -scale_factor1 * sqrt(mode1_member_rha(j)^2 + mode2_member_rha(j)^2);

    end
end


%beam demands/moment for the first mode
beam_member1_rha = zeros(5,2);
beam_member2_rha = zeros(5,2);
Lb = L(:,4);
beam_coef = stiff_coef*(Eb*Ib)/(Li^3);
kbb = [4*(Li^2), 2*(Li^2);2*(Li^2), 4*(Li^2)];
kb = beam_coef*kbb;
for i = 1:5
    beam_member1_rha(i,:) =  kb*dx6([Lb(i), Lb(i+5)]);
    beam_member2_rha(i,:) =  kb*dx7([Lb(i), Lb(i+5)]);

end

%combining the two modes
for j = 1:2
    mode1_beam_rha = beam_member1(:,j);
    mode2_beam_rha = beam_member2(:,j);
    for i = 1:5
        beam_member_combined_rha(i,j) = scale_factor1 * sqrt(mode1_beam_rha(i)^2 + mode2_beam_rha(i)^2);
    end
    
end
% shear demand for RSA
shear_demand_rha = zeros(5,1);
for k = 1:size(beam_member,1)
    shear_demand_rha(k,:) = (beam_member_combined_rha(k,1) + beam_member_combined_rha(k,2))/Li;

end

%% compare results
disp('lateral displacement for ELF vs RSA vs RHA')
displacement = [x_elf, x_rsa, x_rha];
disp(displacement)
disp('ELF vs RSA vs RHA for the column demand (shear demand is the first two columns and the moment demand is the last two columns)')
disp('ELF values')
disp(column_member)
disp('RSA values')
disp(combined_rsa_member)
disp('RHA values')
disp(combined_rha_member)

disp('ELF vs RSA vs RHA for the beam demand ')
disp('ELF values beam moment demand')
disp(beam_member)

disp('ELF values beam shear demand')
disp(shear_demand)

disp('RSA values beam moment demand')
disp(beam_member_combined)

disp('RSA values beam shear demand')
disp(shear_demand1)

disp('RHA values beam moment demand')
disp(beam_member_combined_rha)

disp('RHA values beam shear demand')
disp(shear_demand_rha)