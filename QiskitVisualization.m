%% Comparison of Experimental Data to Cavity-Bloch Equations

% Load data into quadratures
data = csvread('measAmp0.50_driveFreq4.97GHz_measFreq6.99GHz.csv',1);
times = data(:,1);
Ignd = data(:,3);
Qgnd = data(:,2);
Iexc = data(:,5);
Qexc = data(:,4);

% Truncated quadratures for the first 6.8 microseconds, where the data is
% well-behaved
timesTrunc = times(1:69);
IgndTrunc = Ignd(1:69);
IexcTrunc = Iexc(1:69);
QgndTrunc = Qgnd(1:69);
QexcTrunc = Qexc(1:69);

% Fit the ground-state response to a line in the QI-plane to determine the
% rotation angle, which lets us make sure the I response is 0
m = polyfit(QgndTrunc, IgndTrunc,1);
theta = pi-atan2(-m(1),-1);

% Define rotated and rotated truncated quadratures
QgndRot = cos(theta)*Qgnd - sin(theta)*Ignd;
IgndRot = sin(theta)*Qgnd + cos(theta)*Ignd;
QexcRot = cos(theta)*Qexc - sin(theta)*Iexc;
IexcRot = sin(theta)*Qexc + cos(theta)*Iexc;

IgndTruncRot = IgndRot(1:69);
IexcTruncRot = IexcRot(1:69);
QgndTruncRot = QgndRot(1:69);
QexcTruncRot = QexcRot(1:69);

% Now, calculate solution of Cavity-Bloch equations with our new parameter
% values for chi, kappa, and Gamma_1. Code is otherwise the same as that
% used to replicate Bianchetti et al.'s Fig. 4
QCBEgnd = zeros(1,1526);
ICBEgnd = zeros(1,1526);
QCBEexc = zeros(1,1526);
ICBEexc = zeros(1,1526);

%Ground state response
chi = 0.421;
Drm = chi;
Das = -chi;
K = 0.49;
Em = sqrt(K/2);
Om = 0; % For ground state, there is no pi pulse
G1 = 0.00840357;
Gp = 1;

% We make measurement tone Em nonzero only for t > 0

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
sol2 = ode45(CBEqs,[0 15],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = pi;

QCBEgnd(1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QCBEgnd(26:end) = cos(theta)*imag(deval(sol2,0:0.01:15,1))+sin(theta)*real(deval(sol2,0:0.01:15,1));
ICBEgnd(1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
ICBEgnd(26:end) = -sin(theta)*imag(deval(sol2,0:0.01:15,1))+cos(theta)*real(deval(sol2,0:0.01:15,1));

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
sol2 = ode45(CBEqs,[0 15],[deval(sol1,0,1) deval(sol1,0,2) deval(sol1,0,3) deval(sol1,0,4) deval(sol1,0,5) deval(sol1,0,6) deval(sol1,0,7) deval(sol1,0,8)]);
theta = pi;

QCBEexc(1:25) = cos(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+sin(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
QCBEexc(26:end) = cos(theta)*imag(deval(sol2,0:0.01:15,1))+sin(theta)*real(deval(sol2,0:0.01:15,1));
ICBEexc(1:25) = -sin(theta)*imag(deval(sol0,-0.25:0.01:-0.01,1))+cos(theta)*real(deval(sol0,-0.25:0.01:-0.01,1));
ICBEexc(26:end) = -sin(theta)*imag(deval(sol2,0:0.01:15,1))+cos(theta)*real(deval(sol2,0:0.01:15,1));

% Plot our data!
figure()
subplot(2,1,1)
p1=plot(times,QgndRot,'b.','MarkerSize',10);
hold on
p2=plot(-0.25:0.01:15,QCBEgnd/mean(QCBEgnd(26:706))*mean(QgndTruncRot),'k--','LineWidth',1.25);
p3=plot(times,QexcRot,'r.','MarkerSize',10);
plot(-0.25:0.01:15,QCBEexc/mean(QCBEexc(26:706))*mean(QexcTruncRot),'k--','LineWidth',1.25)
xlim([0 15])
xlabel('Time ($\mu$s)','FontSize',14,'LineWidth',1.25,'Interpreter','latex')
ylabel('$Q$','FontSize',14,'LineWidth',1.25,'Interpreter','latex')
legend([p1 p3 p2],{'Ground','Excited','Numerical'},'location','northwest')

subplot(2,1,2)
plot(times,IgndRot,'b.','MarkerSize',10)
hold on
plot(-0.25:0.01:15,ICBEgnd/mean(ICBEgnd(26:706))*mean(IgndTruncRot),'k--','LineWidth',1.25)
plot(times,IexcRot,'r.','MarkerSize',10)
plot(-0.25:0.01:15,ICBEexc/mean(ICBEexc(26:706))*mean(IexcTruncRot),'k--','LineWidth',1.25)
xlim([0 15])
xlabel('Time ($\mu$s)','FontSize',14,'LineWidth',1.25,'Interpreter','latex')
ylabel('$I$','FontSize',14,'LineWidth',1.25,'Interpreter','latex')

figure()
plot(QgndTruncRot,IgndTruncRot,'b.','MarkerSize',10)
hold on
plot(QexcTruncRot,IexcTruncRot,'r.','MarkerSize',10)
plot(QCBEgnd/mean(QCBEgnd(26:706))*mean(QgndTruncRot),ICBEgnd/mean(ICBEgnd(26:706))*mean(IgndTruncRot),'k--','LineWidth',1.25)
plot(QCBEexc/mean(QCBEexc(26:706))*mean(QexcTruncRot),ICBEexc/mean(ICBEexc(26:706))*mean(IexcTruncRot),'k--','LineWidth',1.25)
xlabel('Q','FontSize',14,'Interpreter','latex')
ylabel('I','FontSize',14,'Interpreter','latex')

% Define reference value for I excited state response for next section
IexcRef = zeros(1,20);

% Time values don't quite match up, where they don't we average adjacent
% values
for i = 1:20
    if mod(i,2) == 0
        IexcRef(i) = abs((IexcRot(3+5*(i-2)/2) + IexcRot(4+5*(i-2)/2))/2);
    else
        IexcRef(i) = abs(IexcRot(1+5*(i-1)/2));
    end
end

%% Experimental Reconstruction of Qubit State
% Relies on "IexcRef" variable from previous section

% Load data into quadratures
data = csvread('rabiDrives_driveFreq4.97GHz_measFreq6.99GHz.csv',1);
times = 0.1:0.25:4.85;
amps = data(:,1);
Imat = real(data(:,2:end));
Qmat = imag(data(:,2:end));

% Fit the ground-state response to a line in the QI-plane to determine the
% rotation angle, which lets us make sure the I response is 0
m = polyfit(Qmat(1,:), Imat(1,:),1);
theta = pi-atan2(-m(1),-1);

% Define rotated quadrature (we only need I)
IRot = abs(sin(theta)*Qmat + cos(theta)*Imat);

% Calculate probability we're in the excited state based on ref values
pe = zeros(1,length(amps));
for i = 1:length(pe)
pe(i)=(0.25/4.75)*sum(IRot(i,:)./IexcRef);
end

% Plot!
plot(amps,pe,'k.','MarkerSize',20)
xVals = 0:0.01:0.8;
yVals = (1-cos(11.78*xVals))/2;
hold on
plot(xVals,yVals,'r--','LineWidth',1.25)
legend({'Experimental','Best-Fit Cosine'})
xlabel('Drive Pulse Amplitude','FontSize',14,'Interpreter','latex')
ylabel('$p_e$','FontSize',14,'Interpreter','latex')