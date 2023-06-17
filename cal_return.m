
function output = cal_return(Width,test_width,weight_final)
%open up the actual return and compare
r=readtable('data.xlsx');
overall_return_one=0;
overall_return=zeros(test_width,1);
actual=r{end-test_width+1:end,2:end};
for j=1:test_width
for i=1:Width
    overall_return_one=overall_return_one+weight_final(j,i)*actual(j,i);
end
overall_return(j)=overall_return_one;
overall_return_one=0;
end
output=overall_return;
end