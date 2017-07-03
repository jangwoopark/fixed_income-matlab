function values = fitTwoFactorVas(params)
ax = params(1); ay = params(2); bx = params(3);by = params(4); sx = params(5);sy = params(6); 
global data;
CMT25 = data(:,11)/100;
CMT2 = data(:,12)/100;
CMT3 =data(:,13)/100;
CMT5 = data(:,14)/100;
CMT7 = data(:,15)/100;
CMT10 = data(:,16)/100;

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
ti = [0.25:0.25:10]; count = 1;
while(count <= length(X))
        for t=1:length(ti)
       D(count,t) = A(ax,bx,sx,ti(t))*A(ay,by,sy,ti(t))*exp(-B(bx,ti(t))*X(count)-B(by,ti(t))*Y(count));
       Yield(count,t) = (D(count,t))^(-1/ti(t))-1;
        end
    count = count+1;
end

% compare the yields to the values given and calculate RMSE
points = [1 8 12 20 28 40];
%{
E2 = sum( (Yield(:,points(1)) - CMT2).^(2) );
E3 = sum( (Yield(:,points(2)) - CMT3).^(2) );
E5 = sum( (Yield(:,points(3)) - CMT5).^(2) );
E7 = sum( (Yield(:,points(4)) - CMT7).^(2) );

RMSE = sqrt(E2 + E3 + E5 + E7)/4;

E2 = (Yield(:,points(1)) - CMT2);
E3 = (Yield(:,points(2)) - CMT3);
E5 = (Yield(:,points(3)) - CMT5);
E7 = (Yield(:,points(4)) - CMT7);
%}
E25 = (Yield(:,points(1)) - CMT25);
E2 = (Yield(:,points(2)) - CMT2);
E3 = (Yield(:,points(3)) - CMT3);
E5 = (Yield(:,points(4)) - CMT5);
E7 = (Yield(:,points(4)) - CMT7);
E10 = (Yield(:,points(4)) - CMT10);

values = [E25; E2; E3; E5; E7; E10];

end