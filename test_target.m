clear all
data=readtable('data.xlsx');
Width=width(data)-1;
test_width=50;
weight_1=zeros(test_width,Width);
weight_15=zeros(test_width,Width);
weight_05=zeros(test_width,Width);
w_1=zeros(test_width,Width);
w_2=zeros(test_width,Width);
mu=zeros(1,Width);
epsilon=0.0004;
for k=1:test_width
r=data{5+k:end+k-test_width-1,2:end};
n=length(r);
mu=mean(r);
R_target=mean(mu);
%setting up the return target
R=epsilon+R_target;
R_15=1.5*R;
R_05=0.5*R;
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
        t_1 >= (1/n)*sum(y_1)+epsilon;
        t_1 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon >= R;
        sum(w) == 1;
        w>=0;
        w<=1;
cvx_end
w_1(k,:)=w;
cvx_begin
  variable t_2(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_2 )
    subject to
        for i=1:n
        y_1(i) >= r(i,:)*w-R;
        y_1(i) >= -r(i,:)*w+R;
        y_2(i)>= mu*w-r(i,:)*w+epsilon ;    
        y_2(i) >= -mu*w+r(i,:)*w-epsilon;
        end
        w>=0;
        w<=1;
        t_2 >= (1/n)*sum(y_1)+epsilon;
        t_2 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon <= R;
        mu*w+epsilon >= R;
        sum(w) == 1;
cvx_end
w_2(k,:)=w;
if t_1<=t_2
    weight_1(k,:)=w_1(k,:);
else
    weight_1(k,:)=w_2(k,:);
end
%model with 1.5 target
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
        t_1 >= (1/n)*sum(y_1)+epsilon;
        t_1 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon >= R_15;
        sum(w) == 1;
        w>=0;
        w<=1;
cvx_end
w_1(k,:)=w;
cvx_begin
  variable t_2(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_2 )
    subject to
        for i=1:n
        y_1(i) >= r(i,:)*w-R_15;
        y_1(i) >= -r(i,:)*w+R_15;
        y_2(i)>= mu*w-r(i,:)*w+epsilon ;    
        y_2(i) >= -mu*w+r(i,:)*w-epsilon;
        end
        w>=0;
        w<=1;
        t_2 >= (1/n)*sum(y_1)+epsilon;
        t_2 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon <= R_15;
        mu*w+epsilon >= R_15;
        sum(w) == 1;
cvx_end
w_2(k,:)=w;
if t_1<=t_2
    weight_15(k,:)=w_1(k,:);
else
    weight_15(k,:)=w_2(k,:);
end
%model with 0.5 target
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
        t_1 >= (1/n)*sum(y_1)+epsilon;
        t_1 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon >= R_05;
        sum(w) == 1;
        w>=0;
        w<=1;
cvx_end
if t_1<=t_2
    weight_05(k,:)=w_1(k,:);
else
    weight_05(k,:)=w_2(k,:);
end
cvx_begin
  variable t_2(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_2 )
    subject to
        for i=1:n
        y_1(i) >= r(i,:)*w-R_05;
        y_1(i) >= -r(i,:)*w+R_05;
        y_2(i)>= mu*w-r(i,:)*w+epsilon ;    
        y_2(i) >= -mu*w+r(i,:)*w-epsilon;
        end
        w>=0;
        w<=1;
        t_2 >= (1/n)*sum(y_1)+epsilon;
        t_2 >= (1/n)*sum(y_2)+epsilon;
        mu*w-epsilon <= R_05;
        mu*w+epsilon >= R_05;
        sum(w) == 1;
cvx_end
w_2(k,:)=w;
if t_1<=t_2
    weight_05(k,:)=w_1(k,:);
else
    weight_05(k,:)=w_2(k,:);
end
end
%%
weight_1(isnan(weight_1))=0;
weight_15(isnan(weight_15))=0;
weight_05(isnan(weight_1))=0;
return_1=cal_return(Width,test_width,weight_1);
return_15=cal_return(Width,test_width,weight_15);
return_05=cal_return(Width,test_width,weight_05);
%calculate cumulative return
c_1=cal_creturn(return_1);
c_15=cal_creturn(return_15);
c_05=cal_creturn(return_05);
%plot graph
x=1:1:test_width+1;
plot(x,c_1,'-*b',x,c_15,'magenta',x,c_05,'black'); %line,color and marker
xlabel('days') ;
ylabel('cummulative return');
legend('1 R','1.5 R','0.5 R');