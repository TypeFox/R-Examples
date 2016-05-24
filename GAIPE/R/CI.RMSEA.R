CI.RMSEA <-
function(rmsea,df,N,clevel=.95){
if(length(rmsea)!=1|rmsea<0){stop('rmsea has to be correctly specified.')}
if(length(df)!=1|df<=0){stop('df has to be correctly specified.')}
if(length(N)!=1|N<=0){stop('N has to be correctly specified.')}
if(length(clevel)!=1|clevel<=0|clevel>=1){stop('clevel has to be correctly specified.')}
tol=.0000001;T=rmsea^2*df*(N-1)+df;C=qchisq(1-(1-clevel)/2,df)
if((T<C|rmsea==0)&N>1){;start=rmsea;O_U_P=(1-clevel)/2;P_LU=1
A=c(.01,.001,.0001,.00001,.000001,.0000001)
for(i in c(1,3,5)){;while(P_LU>O_U_P){
P_LU=pchisq(T,df,start^2*df*(N-1));if(P_LU>O_U_P){start=start+A[i]}}
while(P_LU<O_U_P){;P_LU=pchisq(T,df,start^2*df*(N-1))
if(P_LU<O_U_P){start=start-A[i+1]};if(start<0){start=0}}}
k=0;while(abs(P_LU-O_U_P)>tol&k<20){;P_LU=pchisq(T,df,start^2*df*(N-1))
if(P_LU>O_U_P){start=start+.00000001};if(P_LU<O_U_P){start=start-.00000001}
k=k+1;if(start<0){start=0}};if(k==20){;start1=start;start2=start+.00000001
s1=pchisq(T,df,start1^2*df*(N-1));s2=pchisq(T,df,start2^2*df*(N-1))
start=c(start1,start2)[which((abs(c(s1,s2)-O_U_P))==min(abs(c(s1,s2)-O_U_P)))]}
L.rmsea=0;U.rmsea=start};if(T>C&N>1){;O_U_P=(1-clevel)/2;P_LU=1;start=rmsea
A=c(.01,.001,.0001,.00001,.000001,.0000001)
for(i in c(1,3,5)){;while(P_LU>O_U_P){
P_LU=pchisq(T,df,start^2*df*(N-1));if(P_LU>O_U_P){start=start+A[i]}}
while(P_LU<O_U_P){;P_LU=pchisq(T,df,start^2*df*(N-1))
if(P_LU<O_U_P){start=start-A[i+1]};if(start<0){start=0}}}
k=0;while(abs(P_LU-O_U_P)>tol&k<20){;P_LU=pchisq(T,df,start^2*df*(N-1))
if(P_LU>O_U_P){start=start+.00000001};if(P_LU<O_U_P){start=start-.00000001}
k=k+1;if(start<0){start=0}};if(k==20){;start1=start;start2=start+.00000001
s1=pchisq(T,df,start1^2*df*(N-1));s2=pchisq(T,df,start2^2*df*(N-1))
start=c(start1,start2)[which((abs(c(s1,s2)-O_U_P))==min(abs(c(s1,s2)-O_U_P)))]}
U.rmsea=start;start=0;O_L_P=1-(1-clevel)/2;P_LL=1
A=c(.01,.001,.0001,.00001,.000001,.0000001)
for(i in c(1,3,5)){;while(P_LL>O_L_P){
P_LL=pchisq(T,df,start^2*df*(N-1));if(P_LL>O_L_P){start=start+A[i]}}
while(P_LL<O_L_P){;P_LL=pchisq(T,df,start^2*df*(N-1))
if(P_LL<O_L_P){start=start-A[i+1]};if(start<0){start=0}}}
k=0;while(abs(P_LL-O_L_P)>tol&k<20){;P_LL=pchisq(T,df,start^2*df*(N-1))
if(P_LL>O_L_P){start=start+.00000001};if(P_LL<O_L_P){start=start-.00000001}
k=k+1;if(start<0){start=0}};if(k==20){;start1=start;start2=start+.00000001
s1=pchisq(T,df,start1^2*df*(N-1));s2=pchisq(T,df,start2^2*df*(N-1))
start=c(start1,start2)[which((abs(c(s1,s2)-O_U_P))==min(abs(c(s1,s2)-O_U_P)))]}
L.rmsea=start};if(N==1){L.rmsea=0;U.rmsea=Inf}
list(Lower.CI=L.rmsea,RMSEA=rmsea,Upper.CI=U.rmsea)}
