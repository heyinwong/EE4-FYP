clear all
data=readtable('data.xlsx');
Width=width(data)-1;
test_width=50;
weight=zeros(test_width,Width);
w_1=zeros(test_width,Width);
w_2=zeros(test_width,Width);
mu=zeros(1,Width);
steps=100;
epsilon=0.0075;
c=zeros(steps,1);
R=0;
for z=1:steps
for k=1:test_width
r=data{5+k:end+k-test_width-1,2:end};
n=length(r);
mu=mean(r);
%setting up the return target
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
    weight(k,:)=w_1(k,:);
else
    weight(k,:)=w_2(k,:);
end
end
%%
weight(isnan(weight))=0;
return_1=cal_return(Width,test_width,weight);
%calculate cumulative return
c_overall=cal_creturn(return_1);
c(z)=c_overall(end);
R=R+0.0005;
end
%plot graph
%% 
x=0:0.0005:0.0495;
plot(x,c,'-*b'); 
xline(0.0075,'r');
xline(0.014,'black')
xlabel('target return') ;
ylabel('Cummulative Return');
legend('Return target','target = eplisode','10 X average stock return'); 
a=mean(mu);