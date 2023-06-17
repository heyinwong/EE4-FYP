clear all
data=readtable('test_data_2.xlsx','PreserveVariableNames',true);
Width=width(data)-1;
test_width=30;
weight_final=zeros(test_width,Width);
w_1=zeros(test_width,Width);
w_2=zeros(test_width,Width);
mu=zeros(1,Width);
epsilon=0;
steps=100;
p=zeros(steps,1);
Eplison=zeros(steps,1);
for z=1:steps
Eplison(z)=epsilon;
for k=1:test_width
r=data{k:end+k-test_width-1,2:end};
r=log(r+1);
mu=mean(r);
%setting up the return target
if mean(mu)<=epsilon
R_target=epsilon;
else 
R_target=mean(mu);
end
n=length(r);
%Obtain weighting through two cvx optimisations
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
end
weight_final(isnan(weight_final))=0;
final_return=cal_return(Width,test_width,weight_final);
p(z)=mean(final_return);
epsilon=epsilon+0.0001;
end
%% 
%plot graph
plot(Eplison,p,'-*b'); %line,color and marker
xlabel('Eplison') ;
ylabel('mean return');