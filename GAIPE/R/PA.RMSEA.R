PA.RMSEA <-
function(df,method=c("exact.fit","close.fit","not.close.fit"),
H0rmsea,HArmsea,power=.8,alpha=.05){
if(length(df)!=1|df<=0){stop('df has to be correctly specified.')}
if(length(method)!=1|!is.element(method,c("exact.fit","close.fit","not.close.fit")))
stop('method must be one of "exact.fit","close.fit", or "not.close.fit".')
if(length(H0rmsea)!=1|H0rmsea<0){stop('H0rmsea has to be correctly specified.')}
if(length(HArmsea)!=1|HArmsea<0){stop('HArmsea has to be correctly specified.')}
if(length(power)!=1|power<=0|power>1)
{stop('power has to be correctly specified.')}
if(length(alpha)!=1|alpha<=0|alpha>1)
{stop('alpha has to be correctly specified.')}
p=0;N=50;if(method=="not.close.fit"){;if(HArmsea>=H0rmsea)
stop("For the test of not-close fit, H0rmsea has to be larger than HArmsea.")
while(p<power){
p=pchisq(qchisq(alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N+10};while(p>power){
p=pchisq(qchisq(alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N-1};N=N+2};if(method=="close.fit"){;if(HArmsea<=H0rmsea)
stop("For the test of close fit, H0rmsea has to be smaller than HArmsea.")
while(p<power){
p=1-pchisq(qchisq(1-alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N+10};while(p>power){
p=1-pchisq(qchisq(1-alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N-1};N=N+2};if(method=="exact.fit"){;if(H0rmsea!=0)
stop("For the test of exact fit, H0rmsea has to be 0.");while(p<power){
p=1-pchisq(qchisq(1-alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N+10};while(p>power){
p=1-pchisq(qchisq(1-alpha,df,H0rmsea^2*df*(N-1)),df,HArmsea^2*df*(N-1))
N=N-1};N=N+2};N}
