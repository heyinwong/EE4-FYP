clear all
%extract info from dataset
%data=readtable('data.xlsx');
time_horizon=50;
P_size=500;
t=zeros(P_size,1);
data=-1+2*rand(time_horizon,P_size);
for Width=1:P_size
%r=data{6:5+time_horizon,2:Width+1};
r=data(:,1:Width);
mu=mean(r);
%Initilise variables
epsilon=0.0005;
%Set up the return target
if mean(mu)<=epsilon
R_target=epsilon;
else 
R_target=mean(mu);
end
n=time_horizon;
%Obtain weighting through two cvx optimisations
t_start=tic;
cvx_begin
  variable t_1(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_1 )
    subject to
        y_1>=mu*w-r*w-epsilon;
        y_1>=-mu*w+r*w+epsilon;
        y_2>=mu*w-r*w+epsilon;
        y_2>=-mu*w+r*w-epsilon;
        w>=0;
        w<=1;
        t_1 >= (1/n)*sum(y_1)+epsilon;
        t_1 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon >= R_target;
        sum(w) == 1;
cvx_end
w_1=w;
cvx_begin
  variable t_2(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_2 )
    subject to
        y_1>=r*w-R_target;
        y_1>=-r*w+R_target;
        y_2>=mu*w-r*w+epsilon;
        y_2>=-mu*w+r*w-epsilon;
        t_2 >= (1/n)*sum(y_1)+epsilon;
        t_2 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon <= R_target;
        mu*w+epsilon >= R_target;
        w>=0;
        w<=1;
        sum(w) == 1;
cvx_end
w_2=w;
t_end=toc(t_start);
t(Width)=t_end;
end
%% 
a=sum(t);
x=1:1:P_size;
plot(x,t,'-b');
xlabel('Portfolio Size') ;
ylabel('Times/s');
title('Analysis of Portfolio Size on Computation Time')
grid on