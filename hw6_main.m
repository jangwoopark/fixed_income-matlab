 
[data, text] = xlsread('homework5data.xlsx');
global data;

 opt = lsqnonlin(@fitTwoFactorVas,[0.0039 0.1882 0.1983 0.0038 0.0125 0.0327],[0 0 0 0 0 0]);

CMT25 = data(:,11)/100;
CMT2 = data(:,12)/100
CMT3 =data(:,13)/100
CMT5 = data(:,14)/100
CMT7 = data(:,15)/100
CMT10 = data(:,16)/100;

% Golden Numbers 0.0039 0.1882 0.1983 0.0038 0.0125 0.0327


%clear text data;
% pick inital values
ax = opt(1);%0.0089;
ay = opt(2);%-0.0022;
bx = opt(3);%0.1059;
by = opt(4);%0.0447;
sx = opt(5);%0.0206;
sy = opt(6)%0.0109;

testpoints = [ax ay bx by sx sy]

% caculate the implied coefficents for CMT 0.25 and CMT 10
Ax1 = A(ax,bx,sx,0.25);
Bx1 = B(bx,0.25);
Ay1 = A(ay,by,sy,0.25);
By1 = B(by,0.25);
Ax2 = A(ax,bx,sx,10);
Bx2 = B(bx,10);
Ay2 = A(ay,by,sy,10);
By2 = B(by,10);

K = [Bx1/0.25 By1/0.25; Bx2/10 By2/10];
% back out X and Y at each time step using CMT yields
for i=1:length(CMT25)
Vars = K^(-1)*([CMT25(i) + log(Ax1)/0.25 + log(Ay1)/0.25 CMT10(i) + log(Ax2)/10 + log(Ay2)/10]');
X(i,1) = Vars(1);
Y(i,1) = Vars(2);
end

D = [];
% create D(T) function at each point in time from 0.5 to 10
% and obtain the yield amounts
ti = [0.5:0.5:7]; count = 1;
while(count <= length(X))
        for t=1:length(ti)
       D(count,t) = A(ax,bx,sx,ti(t))*A(ay,by,sy,ti(t))*exp(-B(bx,ti(t))*X(count)-B(by,ti(t))*Y(count));
       Yield(count,t) = (D(count,t))^(-1/ti(t))-1;
        end
    count = count+1;
end
% Generate t = 0.25 and t = 10
count = 1;
while(count <= length(X))
    D25 = A(ax,bx,sx,0.25)*A(ay,by,sy,0.25)*exp(-B(bx,0.25)*X(count)-B(by,0.25)*Y(count));
       D10 = A(ax,bx,sx,10)*A(ay,by,sy,10)*exp(-B(bx,10)*X(count)-B(by,10)*Y(count));
       Yield10(count) = (D10)^(-1/10)-1;
       Yield25(count) = (D25)^(-1/0.25)-1;
    count = count+1;
end

% compare the yields to the values given and calculate RMSE
points = [4 6 10 14];
E25 = (Yield25' - CMT25);
E2 = (Yield(:,points(1)) - CMT2);
E3 = (Yield(:,points(2)) - CMT3);
E5 = (Yield(:,points(3)) - CMT5);
E7 = (Yield(:,points(4)) - CMT7);
E10 = (Yield10' - CMT10);

RMSE = sqrt(sum(E2.^2 + E3.^2 + E5.^2 + E7.^2))/4;

% Question 5 - fit to implied moments
figure(1);
subplot(1,2,1); plot(1:length(X), X); title('Implied X Series') ;
subplot(1,2,2); plot(1:length(X), Y); title('Implied Y Series');

smeanX  = mean(X);     smeanY = mean(Y);
sSDX = std(X);         sSDY = std(Y);
imeanX = ax/bx;       imeanY = ay/by;
iSDX = sqrt( sx^2/(2*bx));   iSDY = sqrt( sy^2/(2*by));

display('Question 5: Moment matching')
display(' ');
display(['X: sm = ' num2str(smeanX) '     im = ' num2str(imeanX)]);
display(['   sSD = ' num2str(sSDX) '   iSD = ' num2str(iSDX)]);
display(' ');
display(['Y: sm = ' num2str(smeanY) '     im = ' num2str(imeanY)]);
display(['   sSD = ' num2str(sSDY) '   iSD = ' num2str(iSDY)]);

% Question 6

% Generate t = 0.25 and t = 10
while(count <= length(X))
    D25 = A(ax,bx,sx,0.25)*A(ay,by,sy,0.25)*exp(-B(bx,0.25)*X(count)-B(by,0.25)*Y(count));
       D10 = A(ax,bx,sx,10)*A(ay,by,sy,10)*exp(-B(bx,10)*X(count)-B(by,10)*Y(count));
       Yield10(count) = (D10)^(-1/10)-1;
       Yield25(count) = (D25)^(-1/0.25)-1;
    count = count+1;
end


figure(2)
subplot(2,3,1); plot(1:length(E25), E25); title('3 Month Maturity'); 
subplot(2,3,2); plot(1:length(E2), E2); title('2 Year Maturity');
subplot(2,3,3); plot(1:length(E3), E3); title('3 Year Maturity');
subplot(2,3,4); plot(1:length(E5), E5); title('5 Year Maturity');
subplot(2,3,5); plot(1:length(E7), E7); title('7 Year Maturity');
subplot(2,3,6); plot(1:length(E10), E10); title('10 Year Maturity');

% fitting AR models to the deviations

%{
figure(3)
subplot(2,3,1); parcorr(E25,400); 
subplot(2,3,2); parcorr(E2,400)
subplot(2,3,3); parcorr(E3,400)
subplot(2,3,4); parcorr(E5,400)
subplot(2,3,5); parcorr(E7,400)
subplot(2,3,6); parcorr(E10,400)

figure(4)
subplot(2,3,1); autocorr(E25,400); 
subplot(2,3,2); autocorr(E2,400);
subplot(2,3,3); autocorr(E3,400);
subplot(2,3,4); autocorr(E5,400);
subplot(2,3,5); autocorr(E7,400);
subplot(2,3,6); autocorr(E10,400);

%}
