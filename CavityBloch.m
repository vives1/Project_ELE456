% Code for solving Cavity-Bloch equations + generating numerical figures
% John McEnany, Harsh Babla, Marti Vives
% Figure numbers refer to figures in Bianchetti et al., not our report

%% FIGURE 2

% This is a pretty standard run of the Cavity-Bloch equations, with a pi
% pulse and a measurement tone. This code is repeated in large part for the
% other experiments, so it's commented in detail here, while the other
% sections focus on the differences from previous figures.

%Define arrays to store quadrature values
Qvals = zeros(601,226);
Ivals = zeros(601,226);
k=0;
for Drm = (3:-0.01:-3)*2*pi % Iterate across Delta_rm
k=k+1;
% Time measured in microseconds
chi = 0.69*2*pi; % Cavity pull
Das = -3*chi; % Drive-qubit off resonance, Delta_as -- will be updated once we know <n>
K = 2*pi*1.69; % Photon decay rate kappa
Em = sqrt(K/2); % Measurement tone amplitude
Om = 314.6; % Drive tone amplitude (found by doing an amplitude sweep to produce a pi pulse
G1 = 2*pi*0.19; % Qubit decay rate gamma_1
Gp = 1; % Qubit dephasing (arbitrary value)

% Variables are a, sigma_z, sigma_x, sigma_y, a*sigma_z, a*sigma_x, a*sigma_y, n
% Drive tone amplitude is set to be nonzero only when -10 ns < t < 0
% Other than that, just transcribing Cavity-Bloch equations
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*imag(Y(1))-K*Y(8)];

% Initial conditions, start a in arbitrary value and qubit in GS
% This run is only to determine steady-state value of a and n
inits = [1 -1 0 0 -1 0 0 1];

% First solve region before pi pulse, then during, then after
% Otherwise solver's adaptive time steps will skip over some/all of pulse
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);

% Adjust Delta_as based on steady-state <n>
Das = -2*chi*(deval(sol2,10,8)+1/2);

% Run solver again with steady-state initial conditions
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*imag(Y(1))-K*Y(8)];

inits = [deval(sol2,10,1) -1 0 0 -1*(deval(sol2,10,1)) 0 0 deval(sol2,10,8)];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);

% Define theta to rotate our quadrature axes so that I starts at 0 and Q is
% positive
theta = atan2(real(deval(sol2,10,1)),imag(deval(sol2,10,1)));

% Rotate quadrature axes
Qvals(k,1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
Qvals(k,26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
Ivals(k,1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
Ivals(k,26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));
end

% Plot!
figure()
subplot(2,2,1)
plot(-0.25:0.01:2,Qvals(230,:),'k','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('Q''','Interpreter','latex','FontSize',14)
xlim([-0.25 2])
subplot(2,2,2)
plot(-0.25:0.01:2,Ivals(230,:),'k','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('I''','Interpreter','latex','FontSize',14)
xlim([-0.25 2])
subplot(2,2,3)
imagesc(Qvals,[0 max(max(Qvals))])
colormap([ones(101,1) (1:-0.01:0)' (1:-0.01:0)'])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})

subplot(2,2,4)
imagesc(Ivals,[0 max(max(Ivals))])
colormap([ones(101,1) (1:-0.01:0)' (1:-0.01:0)'])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})
%% FIGURE 3

% Here we basically do the same thing, but plot only a slice of values at a
% given time for a range of Delta_rm's. We also don't perform a special
% rotation (which wasn't made super clear in the paper), and we do a run
% where the qubit starts in the excited state.

Qvals = zeros(601,226);
Ivals = zeros(601,226);
Qexc = zeros(601,1);
Iexc = zeros(601,1);

k=0;
for Drm = (3:-0.01:-3)*2*pi
k=k+1;
chi = 0.71*2*pi;
Das = -3*chi;
K = 2*pi*1.69;
Em = sqrt(K/2);
Om = 314.6;
G1 = 2*pi*0.19;
Gp = 1;

CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*imag(Y(1))-K*Y(8)];

inits = [1 -1 0 0 -1 0 0 1];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);

Das = -2*chi*(deval(sol2,10,8)+1/2);

CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*imag(Y(1))-K*Y(8)];

inits = [deval(sol2,10,1) -1 0 0 -1*(deval(sol2,10,1)) 0 0 1];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = pi; % We need to do this rotation to match the paper

Qvals(k,1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
Qvals(k,26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
Ivals(k,1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
Ivals(k,26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));

% Initial conditions now in excited state
inits = [deval(sol2,10,1) 1 0 0 1*(deval(sol2,10,1)) 0 0 deval(sol2,10,8)];
G1=0;
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*imag(Y(1))-K*Y(8)];

sol = ode45(CBEqs,[0 10],inits);
theta = pi;

Qexc(k) = cos(theta)*imag(deval(sol,10,1))+sin(theta)*real(deval(sol,10,1));
Iexc(k) = -sin(theta)*imag(deval(sol,10,1))+cos(theta)*real(deval(sol,10,1));
end

subplot(2,1,1)
DrmVals = 3:-0.01:-3;
plot(DrmVals,Qvals(:,1),'r','LineWidth',1.25)
hold on
plot(DrmVals,Qvals(:,44),'b','LineWidth',1.25)
plot(DrmVals,Qvals(:,100),'g','LineWidth',1.25)
plot(DrmVals,Qexc,'k--','LineWidth',1.25)
xlabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
ylabel('Q','Interpreter','latex','FontSize',14)
legend({'Ground','180 ns','740 ns','Excited'},'location','northwest','FontSize',14)

subplot(2,1,2)
DrmVals = 3:-0.01:-3;
plot(DrmVals,Ivals(:,1),'r','LineWidth',1.25)
hold on
plot(DrmVals,Ivals(:,44),'b','LineWidth',1.25)
plot(DrmVals,Ivals(:,100),'g','LineWidth',1.25)
plot(DrmVals,Iexc,'k--','LineWidth',1.25)
xlabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
ylabel('I','Interpreter','latex','FontSize',14)
legend({'Ground','180 ns','740 ns','Excited'},'location','northwest','FontSize',14)

%% FIGURE 4

% Now, we do a pulsed measurement, so the measurement tone becomes
% time-dependent. We initialize the system with zero photon occupation
% rather than at steady-state and let it become populated with photons, and
% whether the system is in ground or excited state is determined by whether
% or not we do a pi pulse before turning on the measurement tone. Still no
% fancy rotations; we just take theta = pi to match the paper. It turns out
% we don't have to change the amplitude of the pi pulse at all; the old one
% is still correct. Note that the AC-Stark shift in Delta_as is gone, so
% it's just -chi without a term related to <n>.

QvalsG = zeros(1,226);
IvalsG = zeros(1,226);
QvalsE = zeros(1,226);
IvalsE = zeros(1,226);
%a, sigma_z, sigma_x, sigma_y, a*sigma_z, a*sigma_x, a*sigma_y, n
% Time measured in microseconds
chi = 0.71*2*pi;
Drm = chi;
Das = -chi;
K = 2*pi*1.69;
Em = sqrt(K/2);
Om = 0; % For ground state, there is no pi pulse
G1 = 2*pi*0.19;
Gp = 1;

% Now we make measurement tone Em nonzero only for t > 0

CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0]; % Zero photon occupation, ground state
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = pi;

QvalsG(1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QvalsG(26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
IvalsG(1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
IvalsG(26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));

figure()
subplot(2,1,1)
plot(-0.25:0.01:2,QvalsG,'b','LineWidth',1.25)
xlim([-0.25 2])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('Q','Interpreter','latex','FontSize',14)


subplot(2,1,2)
plot(-0.25:0.01:2,IvalsG,'b','LineWidth',1.25)
xlim([-0.25 2])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('I','Interpreter','latex','FontSize',14)

Om = 314.6; % Now for excited state, we add a pi pulse as usual
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0]; % Start in ground state, then add a pi pulse
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = pi;

QvalsE(1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QvalsE(26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
IvalsE(1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
IvalsE(26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));

subplot(2,1,1)
hold on
plot(-0.25:0.01:2,QvalsE,'r','LineWidth',1.25)
xlim([-0.25 2])

subplot(2,1,2)
hold on
plot(-0.25:0.01:2,IvalsE,'r','LineWidth',1.25)
xlim([-0.25 2])
legend({'Ground','Excited'},'FontSize',14,'location','northeast')

figure()
plot(QvalsG,IvalsG,'b','LineWidth',1.25)
hold on
plot(QvalsE,IvalsE,'r','LineWidth',1.25)
xlabel('Q','Interpreter','latex','FontSize',14)
ylabel('I','Interpreter','latex','FontSize',14)
legend({'Ground','Excited'},'FontSize',14,'location','northeast')

%% FIGURE 5

% Now we're doing the pulsed measurement for a range of value of Delta_rm,
% and bringing back the rotation so that I is 0 at t -> \infty. Nothing
% changes from Figure 4 besides that.

QvalsG = zeros(601,226);
IvalsG = zeros(601,226);
QvalsE = zeros(601,226);
IvalsE = zeros(601,226);

k= 0;
for Drm=(3:-0.01:-3)*2*pi
k=k+1;
chi = 0.71*2*pi;
Das = -chi;
K = 2*pi*1.69;
Em = sqrt(K/2);
Om = 0; % Ground state, no pi pulse
G1 = 2*pi*0.19;
Gp = 1;

CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = atan2(real(deval(sol2,10,1)),imag(deval(sol2,10,1))); % Rotate axes

QvalsG(k,1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QvalsG(k,26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
IvalsG(k,1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
IvalsG(k,26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));

Om = 314.6; % Add a pi pulse to get into excited state
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta=atan2(real(deval(sol2,10,1)),imag(deval(sol2,10,1)));

QvalsE(k,1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QvalsE(k,26:end) = cos(theta)*imag(deval(sol2,0:0.01:2,1))+sin(theta)*real(deval(sol2,0:0.01:2,1));
IvalsE(k,1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
IvalsE(k,26:end) = -sin(theta)*imag(deval(sol2,0:0.01:2,1))+cos(theta)*real(deval(sol2,0:0.01:2,1));
end

figure()
subplot(2,2,3)
plot(-0.25:0.01:2,QvalsG(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,QvalsG(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,QvalsG(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('Q''','Interpreter','latex','FontSize',14)
subplot(2,2,4)
plot(-0.25:0.01:2,IvalsG(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,IvalsG(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,IvalsG(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('I''','Interpreter','latex','FontSize',14)
subplot(2,2,1)
imagesc(QvalsG,[-max(max(abs(QvalsG))) max(max(abs(QvalsG)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})
subplot(2,2,2)
imagesc(IvalsG,[-max(max(abs(IvalsG))) max(max(abs(IvalsG)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})

figure()
subplot(2,2,3)
plot(-0.25:0.01:2,QvalsE(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,QvalsE(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,QvalsE(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('Q''','Interpreter','latex','FontSize',14)
subplot(2,2,4)
plot(-0.25:0.01:2,IvalsE(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,IvalsE(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,IvalsE(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('I''','Interpreter','latex','FontSize',14)
subplot(2,2,1)
imagesc(QvalsE,[-max(max(abs(QvalsE))) max(max(abs(QvalsE)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})
subplot(2,2,2)
imagesc(IvalsE,[-max(max(abs(IvalsE))) max(max(abs(IvalsE)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})

% Difference between excited and ground states
QvalsD = QvalsE - QvalsG;
IvalsD = IvalsE - IvalsG;

figure()
subplot(2,2,3)
plot(-0.25:0.01:2,QvalsD(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,QvalsD(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,QvalsD(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('Q''','Interpreter','latex','FontSize',14)
subplot(2,2,4)
plot(-0.25:0.01:2,IvalsD(371,:),'r','LineWidth',1.25)
hold on
plot(-0.25:0.01:2,IvalsD(271,:),'k','LineWidth',1.25)
plot(-0.25:0.01:2,IvalsD(161,:),'b','LineWidth',1.25)
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('I''','Interpreter','latex','FontSize',14)
subplot(2,2,1)
imagesc(QvalsD,[-max(max(abs(QvalsD))) max(max(abs(QvalsD)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})
subplot(2,2,2)
imagesc(IvalsD,[-max(max(abs(IvalsD))) max(max(abs(IvalsD)))])
colormap([[(0:0.01:1)' (0:0.01:1)' ones(101,1)]; [ones(101,1) (1:-0.01:0)' (1:-0.01:0)']])
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',14)
ylabel('$\Delta_\mathrm{rm}$ (MHz)','Interpreter','latex','FontSize',14)
xticks([26 76 126 176 226])
xticklabels({'0','0.5','1','1.5','2'})
yticks([1 101 201 301 401 501 601])
yticklabels({'3','2','1','0','-1','-2','-3'})
%%

% FIGURE 6

% We're now at the main result of the paper -- using the quadrature
% amplitudes to predict the probability that the qubit is in the excited
% state. We first calculate the reference response of a qubit prepared in
% the ground or excited states, as in Fig. 4. Next, we perform experiments
% with varying drive pulse lengths (which will put the qubit in a different
% superposition of ground and excited states), and compare the actual
% probability that it's in the excited state to the probability calculated
% from the quadrature amplitudes, using different discrete timesteps.

% Note that our parameter values are now different
chi = -1.02*2*pi;
Das = -chi;
K = 2*pi*1.69;
Em = sqrt(K/2); 
G1 = 1/0.86;
Gp = 1;
Drm = chi;

% Ground state:
Om = 0; % No pi pulse
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
mqG = (imag(deval(sol2,0:0.01:2,1))); % Quadrature amplitudes Q and I
miG = (real(deval(sol2,0:0.01:2,1)));

% Excited state:
Om = 314.6; % Pi pulse
CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
    Om*(t >= -0.01 & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
    (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -0.01 & t <= 0)*Y(2); ...
    -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -0.01 & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
    -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
    -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -0.01 & t <= 0)*Y(5); ...
    -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

inits = [0 -1 0 0 0 0 0 0];
sol0 = ode45(CBEqs,[-0.25 -0.01],inits);
sol1 = ode45(CBEqs,[-0.01 0],[deval(sol0,-0.01,1) deval(sol0,-0.01,2) deval(sol0,-0.01,3) deval(sol0,-0.01,4) deval(sol0,-0.01,5) deval(sol0,-0.01,6) deval(sol0,-0.01,7) deval(sol0,-0.01,8)]);
sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
mqE = (imag(deval(sol2,0:0.01:2,1)));
miE = (real(deval(sol2,0:0.01:2,1)));

% Now we run the drive pulse for a variety of times ending at t=0, and
% query the probability that the qubit is in the excited state at t=0
Om = 50*2*pi; % New value of omega as reported in paper

pe = zeros(51,1); % Actual excited-state probability
peEst = zeros(51,1); % Probability estimated from quadratures
peEst2 = zeros(51,1); % Estimated probability with bigger timesteps
pe(1)=0;
peEst(1)=0;
peEst2(1)=0;
k=1;
for tPulse = 0.0001:0.0001:0.05
    k=k+1;
    CBEqs = @(t,Y)[-1i*Drm*Y(1)-1i*chi*Y(5)-1i*Em*(t>0)-(K/2)*Y(1); ...
        Om*(t >= -tPulse & t <= 0)*Y(4)-G1*(1+Y(2)); -(Das+2*chi*(Y(8)+1/2))*Y(4)-(G1/2+Gp)*Y(3); ...
        (Das+2*chi*(Y(8)+1/2))*Y(3)-(G1/2+Gp)*Y(4)-Om*(t >= -tPulse & t <= 0)*Y(2); ...
        -1i*Drm*Y(5)-1i*chi*Y(1)+Om*(t >= -tPulse & t <= 0)*Y(7)-1i*Em*(t>0)*Y(2)-G1*Y(1)-(G1+K/2)*Y(5); ...
        -1i*Drm*Y(6)-(Das+2*chi*(Y(8)+1))*Y(7)-1i*Em*(t>0)*Y(3)-(G1/2+Gp+K/2)*Y(6); ...
        -1i*Drm*Y(7)+(Das+2*chi*(Y(8)+1))*Y(6)-1i*Em*(t>0)*Y(4)-(G1/2+Gp+K/2)*Y(7)-Om*(t >= -tPulse & t <= 0)*Y(5); ...
        -2*Em*(t>0)*imag(Y(1))-K*Y(8)];

    inits = [0 -1 0 0 0 0 0 0];
    sol0 = ode45(CBEqs,[-0.25 -tPulse],inits);
    sol1 = ode45(CBEqs,[-tPulse 0],[deval(sol0,-tPulse,1) deval(sol0,-tPulse,2) deval(sol0,-tPulse,3) deval(sol0,-tPulse,4) deval(sol0,-tPulse,5) deval(sol0,-tPulse,6) deval(sol0,-tPulse,7) deval(sol0,-tPulse,8)]);
    sol2 = ode45(CBEqs,[0 10],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);

    % Quadrature amplitudes to compare to ground/excited states
    mq = (imag(deval(sol2,0:0.01:2,1)));
    mi = (real(deval(sol2,0:0.01:2,1)));
    pe(k) = (deval(sol1,0,2)+1)/2; % Query sigma-z at t=0
    
    % Estimated pe, based on estimate from I quadrature
    peEst(k) = (0.01/1.99)*sum((mi(2:end)-miG(2:end))./(miE(2:end)-miG(2:end)));
    peEst2(k) = 0;
    for j = 2:10:192
        peEst2(k) = peEst2(k) + (0.1/1.91)*((mi(j)-miG(j))/(miE(j)-miG(j)));
    end
end


figure()
plot(1000*(0:0.0001:0.05),pe,'k','LineWidth',1.25)
hold on
plot(1000*(0:0.0001:0.05),peEst,'r','LineWidth',1.25)
plot(1000*(0:0.0001:0.05),peEst2,'b','LineWidth',1.25)
xlabel('Pulse Length (ns)','FontSize',14,'Interpreter','latex')
ylabel('$p_e$','FontSize',20,'Interpreter','latex')
legend({'Actual','Estimated ($\Delta t$ = 0.01 $\mu$s)','Estimated ($\Delta t$ = 0.1 $\mu$s)'},'FontSize',14,'Interpreter','latex')
