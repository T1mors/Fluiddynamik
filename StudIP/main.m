clc
clear all
close all

pathData='k6mf_processed_FD2.txt';
Fs=20000;

%Opening the file
fileID=fopen(pathData,'r');
data=textscan(fileID, '%f','CollectOutput',1);
fclose(fileID);
data=data{1};
dt=1/Fs;
time=[dt:dt:length(data)/Fs]';

%% Exercise 5 (18.05.22)
% mean value
mean_u=sum(data)/length(data);
%mean(data)

% standard deviation
std_u=sqrt(sum(((data-mean_u).^2)/length(data)));
%std(data)

% spatial field u(x) -> x=<u>t
xr=time.*mean_u;
max_xr=max(xr);

figure()
maxplot=10000;
plot(xr(1:maxplot),data(1:maxplot))
xlabel('$x$ [m]', 'interpreter','latex')
ylabel('$u(x)$  [m/s]','interpreter','latex')

%fluctuations
u_prime=data-mean_u;


%% Correlation
max_r=10000;
tVect_reduced=[1:10:max_r];     %For reducing computational time
rVect_reduced=tVect_reduced.*mean_u./Fs;
corr=[];
count=1;

for r=tVect_reduced
    corr(count)=sum((u_prime(1:size(u_prime,1)-r)).*(u_prime(1+r:end)))/(length(u_prime)-r);
    count=count+1;
end

figure()
plot(rVect_reduced,corr)
hold on
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$C_r$  [m/s]','interpreter','latex')
ylim([0 0.03])

%
% [autocorr, lags]=xcorr(u_prime);
% autocorr=autocorr(lags>0);
% figure()
% plot(lags(lags>0), autocorr(lags>0))
% xlabel('$r$ [m]', 'interpreter','latex')
% ylabel('$C_r$  [m/s]','interpreter','latex')


% Normalizing autocorrelation Cr(r)/Cr(0)
corrNorm=corr./max(corr);

figure()
plot(rVect_reduced,corrNorm)
hold on
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$C_r$  [m/s]','interpreter','latex')
ylim([0 1])

% Integral Length
lim_L=find(corrNorm<0,1,'first');
int_L1     = trapz(rVect_reduced(1:lim_L),corrNorm(1:lim_L))

%Exponential fit
fitExp=fit(tVect_reduced',corrNorm','exp1');
corr_exp=fitExp.a*exp(fitExp.b*tVect_reduced);

hold on
plot(rVect_reduced,corr_exp)
int_L_Exp=trapz(rVect_reduced,corr_exp)

% Taylor Length 
% I will calculate the correlation for all points(not reduced) but only for
% 100 data points.
corrTL=[];
count=1;
for r=1:100
    corrTL(count)=sum((u_prime(1:size(u_prime,1)-r)).*(u_prime(1+r:end)))/(length(u_prime)-r);
    count=count+1;
end

% Normalizing autocorrelation Cr(r)/Cr(0)
corrTLNorm=corrTL./max(corrTL);

figure()
plot(xr(1:100),corrTLNorm(1:100))
hold on

% Parabolic fit
% Fit will be done for the first 8 values...
p1=polyfit(xr(1:8), corrTLNorm(1:8), 2);
zG=polyval(p1,xr(1:30));

plot(xr(1:30),zG)
hold on

TaylorLength=sqrt(-corrNorm(1)/p1(1))
%% Structure functions
Sn2=[];
Sn3=[];
Sn4=[];
Sn3c=[];

count=1;
for r=tVect_reduced
    %Sn1(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^1)/(length(data)-r);
    Sn2(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^2)/(length(data)-r);
    Sn3(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^3)/(length(data)-r);
    Sn4(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^4)/(length(data)-r);
    Sn5(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^5)/(length(data)-r);
    Sn6(count)=sum(((data(1:size(data,1)-r))-(data(1+r:end))).^6)/(length(data)-r);
    count=count+1;
end

%% EVEN -structure functions
figure()
plot(rVect_reduced, Sn2, 'DisplayName', '$S_2$')
hold on
plot(rVect_reduced, Sn4, 'DisplayName', '$S_4$')
hold on
plot(rVect_reduced, Sn6, 'DisplayName', '$S_6$')
set(gca,'xscale','log','Yticklabel',[])
set(gca,'yscale','log')
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$S_n $ [a.u]','interpreter','latex')
legend('show','interpreter', 'latex')

rlow=find(rVect_reduced>TaylorLength,1,'first');
rhigh=find(rVect_reduced>int_L1,1,'first');

p1=polyfit(log(rVect_reduced(rlow:rhigh)), log(Sn2(rlow:rhigh)),1);
z1=polyval(p1,log(rVect_reduced(rlow:rhigh)));
a2=p1(1);
hold on
plot(rVect_reduced(rlow:rhigh), exp(z1),'--','DisplayName',strcat('$\propto r^{',num2str(a2,'%.2f'),'}$'))

p1=polyfit(log(rVect_reduced(rlow:rhigh)), log(Sn4(rlow:rhigh)),1);
z1=polyval(p1,log(rVect_reduced(rlow:rhigh)));
a4=p1(1);
hold on
plot(rVect_reduced(rlow:rhigh), exp(z1),'--','DisplayName',strcat('$\propto r^{',num2str(a4,'%.2f'),'}$'))

p1=polyfit(log(rVect_reduced(rlow:rhigh)), log(Sn6(rlow:rhigh)),1);
z1=polyval(p1,log(rVect_reduced(rlow:rhigh)));
a6=p1(1);
hold on
plot(rVect_reduced(rlow:rhigh), exp(z1),'--','DisplayName',strcat('$\propto r^{',num2str(a6,'%.2f'),'}$'))

hold on
xline(rVect_reduced(rlow), 'DisplayName', 'Taylor Length')
hold on
xline(rVect_reduced(rhigh),'DisplayName','Integral Length')

% n vs exponent 
figure()
plot([2 4 6], [a2 a4 a6],'-o')
xlabel('$n, S^n$', 'interpreter','latex')
ylabel('Exponent','interpreter','latex')

% ODD - Structure functions
figure()
plot(rVect_reduced, abs(Sn3), 'DisplayName', '$S_3$')
hold on
plot(rVect_reduced, abs(Sn5), 'DisplayName', '$S_5$')
set(gca,'xscale','log','Yticklabel',[])
set(gca,'yscale','log')
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$S_n $ [a.u]','interpreter','latex')
legend('show','interpreter', 'latex')

%% Exercise 6 08.06.2022

% Bins are defined ror smoothing
increment_bin = ceil((max(range(data))/mean(nanstd(data))*10));
if mod(increment_bin,2) == 0
    increment_bin = increment_bin +1;
end

L=size(data,1);
if mod(L,2)==1
    data=data(1:end-1,:); % This is to make length of data array even
    L=size(data,1);
end

f                 = Fs/2*linspace(0,1,L/2+1).';        % Nyquist frequency of the signal==Fs/2
f                 = f(2:end);                                   % Remove zero Hz component

spek=abs(fft(data,L)).^2/(L*Fs);                      % Energy spectral density(ESD) using a fft
spek              = 2.*spek(2:L/2+1);
%spek              = 2.*spek(2:L/2+1,:);               % FFt will yield half number of unique points (!!!Here different for a matrix!!!)


% !!If moving average is required (as here I will avg over several time series, I am not using moving avg)
% moving average with equally spaced frequency interval in log-space
plot_length         = increment_bin*10;
if mod(plot_length,2) == 0
    plot_length = plot_length +1;
end

% moving average with equally spaced frequency interval in lin-space
% intervall=unique(round(linspace(1,L/2,plot_length)),'stable');
intervall=unique(round(logspace(0,log10(L/2),plot_length)),'stable');
plot_length=length(intervall);      % Number of points for plotting averaged spectrum

% Initializing the arrays/preallocation
spek_smooth=zeros(plot_length-2,size(data,2));
x_achse_f_spektrum=zeros(plot_length-2,1);

% Averaging of spectrum and hence smoothing
for i=2:(plot_length-1)
    x_achse_f_spektrum(i-1,1) = mean(f(intervall(i-1):intervall(i+1)));
    spek_smooth(i-1,:)        = mean(spek(intervall(i-1):intervall(i+1),:));
end

figure()
loglog(f,spek,'-','LineWidth',1.5,'DisplayName','Spectrum')
hold on
loglog(x_achse_f_spektrum,spek_smooth,'-','LineWidth',1.5,'DisplayName','Avg. Spectrum')
hold on
plot(f,f.^(-5/3),'DisplayName','$f^{-5/3}$');
set(gca,'FontSize',16)
xlabel('f [1/s]','interpreter','latex');
ylabel('E(f) [$\frac{m^2}{s}$]', 'interpreter','latex')
legend('show','interpreter','latex')

%Increments
%Values of tau for whose increments will be calculated
% Vector logarithmically distrbuted from 10^0 to 10^2.5

tauVec=[1*dt 10*dt 100*dt  1000*dt]; %[s]
dataIncrement=cell(1,length(tauVec));

for i=1:length(tauVec)
    tau=round(tauVec(i)*Fs); %[steps]
    incr_range  =   (data((tau+1):size(data,1),:) - data(1:(size(data,1)-tau),:));
    std_increments= nanstd(incr_range); % Standard deviation of all increments
    
    dataIncrement{i}= incr_range./std_increments; % Normalized increments (by standard deviation)
end

rangeBin=[-10 10];
nBins=31;
minEvents=2;
binEdges=linspace(rangeBin(1),rangeBin(2),nBins+1);
dx=diff(binEdges);
binCenters=binEdges(1:length(binEdges)-1)+dx;
pX_scale=[];

for i=1:length(tauVec)
    data=dataIncrement{i};
    figHist=figure();
    h1=histogram(data,binEdges,'Normalization','pdf');
    countsX=h1.BinCounts;
    countsX(countsX<minEvents)=nan;
    N_eventsX=nansum(countsX);
    statErrorX=sqrt(countsX)./N_eventsX;
    close(figHist)
    
    pXi= countsX./(N_eventsX*dx);
    pX_scale(:,i)=pXi;
end

figure()
for i=1:length(tauVec)
    scatter(binCenters,1/20^i.*pX_scale(:,i),...
        'DisplayName',strcat('$$r = ',num2str(tauVec(i)*mean_u),'\,\mathrm {m}$$'))
    hold on
end
set(gca,'yscale','log')
yticklabels({'','',''})
xlabel('$u_{r} / \sigma_{u_{r}}$', 'interpreter','latex')
ylabel('$p(u_{r})$ (log scale)','interpreter','latex')
lgd = legend('show','Interpreter','Latex','NumColumns',1,'Location','eastoutside');
xlim(rangeBin)
set(gca,'FontSize',16)

%% Exercise 7  22.06.22
% Local energy disipation rate
visc=1.5e-5;

% Derivative of velocity field
dr=mean_u/Fs;
dudx  =   ((data(2:size(data,1),:) - data(1:(size(data,1)-1),:))./dr).^2;

e_deriv=dudx*2*visc;
e_mean_deriv=mean(e_deriv)

eta=(visc^3/e_mean_deriv)^0.25
TaylorLength2=sqrt(15*visc/e_mean_deriv)*std(data);

rlow=find(rVect_reduced>eta,1,'first');
rhigh=find(rVect_reduced>TaylorLength,1,'first');

% From S3
figure()
hold on
plot(rVect_reduced, abs(Sn3), 'DisplayName', '$S_3$')
hold on
plot(rVect_reduced,rVect_reduced.^(4/5),'DisplayName', '$r^{4/5}$')
hold on
set(gca,'xscale','log','Yticklabel',[])
set(gca,'yscale','log')
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$S_n $ [a.u]','interpreter','latex')
legend('show','interpreter', 'latex')

er_s3=(-5/4).*((Sn3))./rVect_reduced;
idx=find(er_s3==max(er_s3));
e_mean_s3=mean(er_s3(idx-4:idx+4)) % averaging over 8 steps

figure()
plot(rVect_reduced, abs(er_s3), 'DisplayName', '$S_3$')
hold on
yline(max(er_s3))
set(gca,'xscale','log')
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('$\epsilon [m^2/s^3]$','interpreter','latex')

%% Scaled resolved dissipation rate
a=1;
rangeBin=[-10 10];
nBins=81;
minEvents=2;
binEdges=linspace(rangeBin(1),rangeBin(2),nBins+1);
dx=diff(binEdges);
binCenters=binEdges(1:length(binEdges)-1)+dx;
pX_scale=[];
std_er=[];

for dr_inc=dr*[2 10 100 1000 10000]

    step_dx=round(dr_inc*Fs);
    
    up_vec=cumsum(e_deriv(step_dx+1:size(e_deriv,1)));
    low_vec=cumsum(e_deriv(1:size(e_deriv,1)-step_dx));

    c=up_vec-low_vec;
    std_er(a)=std(c);
    c=c./std(c);
    
    figHist=figure();
    h1=histogram(c,binEdges,'Normalization','pdf');
    countsX=h1.BinCounts;
    countsX(countsX<minEvents)=nan;
    N_eventsX=nansum(countsX);
    statErrorX=sqrt(countsX)./N_eventsX;
    close(figHist)
    
    pXi= countsX./(N_eventsX*dx);
    pX_scale(:,a)=pXi;
    a=a+1;
end


dr_inc=dr*[2 10 100 1000 10000];
figure()
for i=1:size(pX_scale,2)
    scatter(binCenters,pX_scale(:,i),...
        'DisplayName',strcat('$$r = ',num2str(dr_inc(i)),'\,\mathrm {m}$$'))
    hold on
end
set(gca,'yscale','log')
yticklabels({'','',''})
xlabel('$\varepsilon_{r} / \sigma_{\varepsilon_{r}}$', 'interpreter','latex')
ylabel('$p(\varepsilon_{r})$ (log scale)','interpreter','latex')
lgd = legend('show','Interpreter','Latex','NumColumns',1,'Location','eastoutside');
set(gca,'FontSize',16)

%Standard deviation vs r
figure()
plot(dr_inc, std_er.^2,'-o')
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel('$\sigma(\varepsilon_{r})^2$', 'interpreter','latex')
xlabel('$r$ [m]','interpreter','latex')
set(gca,'FontSize',16)

p1=polyfit(log(dr_inc), log(std_er.^2),1);
zG=polyval(p1,log(dr_inc));
plot(dr_inc,dr_inc.^(p1(1)))

    
%% Exercise 13.07.22
Flatness_r=Sn4./(3.*(Sn2.^2));

l1=3;
l2=8;
%p1=polyfit(Flatness_r(l1:l2), log(rVect(l1:l2).*mean_u./Fs),1);
p1=polyfit(log(rVect_reduced(l1:l2)), log(Flatness_r(l1:l2)),1);
%zG=0.55.*(rVect(l1:l2).*mean_u./Fs).^-0.185;
zG=(rVect_reduced(l1:l2)).^(p1(1))+p1(2);

figure()
plot(rVect_reduced, Flatness_r,'DisplayName','F(r)','LineWidth',2)
hold on
yline(1,'DisplayName','F=1')
hold on
plot(rVect_reduced(l1:l2),zG,'LineWidth',2)
%hold on
%plot(rVect(l1:l2).*mean_u./Fs,(rVect(l1:l2).*mean_u./Fs).^((-4/9)*miu_intermittent))
xlabel('$r$ [m]', 'interpreter','latex')
ylabel('Flatness [-]','interpreter','latex')
set(gca,'FontSize',14)
legend('show','interpreter', 'latex')
set(gca,'xscale','log')
%set(gca,'yscale','log')
%ylim([0 3])
miu_intermittent=(-9/4)*p1(1)


