clear all
Width=20;
horizon=100;
test_width=10;
data=-1+2*rand(horizon,Width);
weight_final=zeros(test_width,Width);
bench_mad=zeros(test_width,Width);
bench_mv=zeros(test_width,Width);
w_1=zeros(test_width,Width);
w_2=zeros(test_width,Width);
mu=zeros(1,Width);
epsilon=0.01;
for k=1:test_width
r=data(k:end+k-test_width-1,1:end);
mu=mean(r);
R_target=mean(mu);
R=R_target;
n=length(r);
%calculate sigma
s=cov(r);
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
    weight_final(k,:)=w_1(k,:);
else
    weight_final(k,:)=w_2(k,:);
end
%mad_bench without eplison
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
        y_2(i)>= mu*w-r(i,:)*w;    
        y_2(i) >= -mu*w+r(i,:)*w;
        w>=0;
        w<=1;
        end
        t_1 >= (1/n)*sum(y_1);
        t_1 >= (1/n)*sum(y_2);
        mu*w>= R_target;
        sum(w) == 1;
cvx_end
bench_mad(k,:)=w;
%create classical mean variance model
cvx_begin
  variable t(1);
  variable w(Width);
  minimize( t )
    subject to
        t>=w.'*s*w;
        w>=0;
        w<=1;
        sum(w) == 1;
        mu*w>=R_target;
cvx_end
bench_mv(k,:)=w;
end
%%
weight_final(isnan(weight_final))=0;
bench_mad(isnan(bench_mad))=0;
bench_mv(isnan(bench_mad))=0;
final_return=cal_sreturn(Width,test_width,weight_final,data);
%create bench for compare: equally invest each stock
bench_weight=ones(test_width,Width)/Width;
bench_return=cal_sreturn(Width,test_width,bench_weight,data);
mad_return=cal_sreturn(Width,test_width,bench_mad,data);
mv_return=cal_sreturn(Width,test_width,bench_mv,data);
%calculate cumulative return
c=cal_creturn(final_return);
c_bench=cal_creturn(bench_return);
c_mad=cal_creturn(mad_return);
c_mv=cal_creturn(mv_return);
%plot graph
x=1:1:test_width+1;
plot(x,c,'-*b',x,c_bench,'r',x,c_mad,'magenta',x,c_mv,'black'); %line,color and marker
xlabel('days') ;
ylabel('cummulative return/%');
legend('DR-MAD','1/N','MAD','Mean Variance');
grid on
title('Test on synthetic data')
average=mean(final_return);
a_mv=mean(mv_return);
a_stupid=mean(bench_return);
a_mad=mean(mad_return);