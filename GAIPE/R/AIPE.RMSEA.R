AIPE.RMSEA <-
function(rmsea,df,width,clevel=.95){
if(length(width)!=1|width<=0){stop('width has to be correctly specified.')}
A=c(100,50,10,5,2,1)
if(width<.03){A=c(500,200,100,50,10,5,2,1)}
N=50;W=1;for(i in seq(1,length(A),2)){
while(W>width){
CI=CI.RMSEA(rmsea,df,N,clevel);W=CI$U-CI$L;N=N+A[i]}
while(W<width){
CI=CI.RMSEA(rmsea,df,N,clevel);W=CI$U-CI$L;N=N-A[i+1];if(N<=0){N=1}
}};N=N+2;N}
