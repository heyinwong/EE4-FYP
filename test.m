clear all
%extract info from dataset
data=readtable('test_data_2.xlsx','PreserveVariableNames',true);
%Initilise variables
Width=width(data)-1;
r=data{:,2:end};
%range should be 0.05-0.20
epsilon=0.05;
mu=mean(r);
%Set up the return target
R_target=0;
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