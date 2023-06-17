function output = cal_creturn(final_return)
final_return=final_return+1;
final_return=[1;final_return];
final_return=cumprod(final_return);
output=final_return;
