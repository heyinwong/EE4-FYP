clear all
%extract info from dataset
data=readtable('test_data_2.xlsx','PreserveVariableNames',true);
Width=3;
r=data{:,2:4};
mu=mean(r);
%Initilise variables
weight=zeros(100,Width);
epsilon=0;
Epsilon=zeros(25,1);
n=length(r);
for j=1:200
%Set up the return target
R_target=0;
Epsilon(j)=epsilon;
n=length(r);
%Obtain weighting through two cvx optimisations
cvx_begin
  variable t_1(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_1 )
    subject to
        for i=1:n
        y_1(i) >= mu*w-r(i,:)*w-epsilon;
        y_1(i) >= -mu*w+r(i,:)*w+epsilon;
        y_2(i)>= mu*w-r(i,:)*w+epsilon ;    
        y_2(i) >= -mu*w+r(i,:)*w-epsilon;
        end
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
        for i=1:n
        y_1(i) >= r(i,:)*w-R_target;
        y_1(i) >= -r(i,:)*w+R_target;
        y_2(i)>= mu*w-r(i,:)*w+epsilon ;    
        y_2(i) >= -mu*w+r(i)*w-epsilon;
        w>=0;
        w<=1;
        end
        t_2 >= (1/n)*sum(y_1)+epsilon;
        t_2 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon <= R_target;
        mu*w+epsilon >= R_target;
        sum(w) == 1;
cvx_end
w_2=w;
if t_1<=t_2
    weight(j,:)=w_1;
else
    weight(j,:)=w_2;
end
epsilon=epsilon+2/10000;
end
%% 
weight(isnan(weight))=0;
cvx_begin
  variable t_1(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_1 )
    subject to
        for i=1:n
        y_1(i) >= mu*w-r(i,:)*w;
        y_1(i) >= -mu*w+r(i,:)*w;
        y_2(i)>= mu*w-r(i,:)*w ;    
        y_2(i) >= -mu*w+r(i,:)*w;
        end
        w>=0;
        w<=1;
        t_1 >= (1/n)*sum(y_1);
        t_1 >= (1/n)*sum(y_2);
        mu*w >= 0;
        sum(w) == 1;
cvx_end
bench_mad=w;
%% 
r=zeros(25,1)
for i=1:100
    r(i)=corr(weight(i,:).',w)
end
plot(Epsilon,r);
title('Correlation between DR-MAD and MAD')
xlabel('Epsilon') 
ylabel('Correlation coefficient') 