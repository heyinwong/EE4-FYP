clear all
data=readtable('data.xlsx');
Width=width(data)-1;
test_width=20;
miss_width=35;
weight_final=zeros(test_width,Width);
bench_mad=zeros(test_width,Width);
bench_mv=zeros(test_width,Width);
w_1=zeros(test_width,Width);
w_2=zeros(test_width,Width);
mu=zeros(1,Width);
epsilon=0.0075;
for k=1:test_width
r=data{5+k:end+k-test_width-1-miss_width,2:end};
mu=mean(r);
if mean(mu) < epsilon
    R_target=epsilon;
else
    R_target=mean(mu);
end
R=mean(mu);
n=length(r);
s=cov(r);
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
w_1(k,:)=w;
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
w_2(k,:)=w;
if t_1<=t_2
    weight_final(k,:)=w_1(k,:);
else
    weight_final(k,:)=w_2(k,:);
end
%mad_bench without epsilon
cvx_begin
  variable t_1(1);
  variable y_2(n);
  variable y_1(n);
  variable w(Width);
  minimize( t_1 )
    subject to
        y_1>=mu*w-r*w;
        y_1>=-mu*w+r*w;
        y_2>=mu*w-r*w;
        y_2>=-mu*w+r*w;
        w>=0;
        w<=1;
        t_1 >= (1/n)*sum(y_1);
        t_1 >= (1/n)*sum(y_2);
        mu*w >= R;
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
        mu*w>=R;
cvx_end
bench_mv(k,:)=w;
end
%%
weight_final(isnan(weight_final))=0;
bench_mad(isnan(bench_mad))=0;
bench_mv(isnan(bench_mad))=0;
final_return=cal_sreturn(Width,test_width,weight_final,miss_width);
%create bench for compare: equally invest each stock
bench_weight=ones(test_width,Width)/Width;
bench_return=cal_sreturn(Width,test_width,bench_weight,miss_width);
mad_return=cal_sreturn(Width,test_width,bench_mad,miss_width);
mv_return=cal_sreturn(Width,test_width,bench_mv,miss_width);
%calculate cumulative return
c=cal_creturn(final_return);
c_bench=cal_creturn(bench_return);
c_mad=cal_creturn(mad_return);
c_mv=cal_creturn(mv_return);
%plot graph
x=1:1:test_width+1;
plot(x,c,'-*b',x,c_bench,'r',x,c_mad,'magenta',x,c_mv,'black'); %line,color and marker
xlabel('days') ;
ylabel('cummulative return');
legend('DR-MAD','1/N','MAD','Mean Variance');
grid on
title('Out-of-sample performance')
average=mean(final_return);
a_mv=mean(mv_return);
a_stupid=mean(bench_return);
a_mad=mean(mad_return);