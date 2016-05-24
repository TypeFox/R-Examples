clopper.pearson.ci <-
function(k,n,alpha=0.1,CI="upper"){
l<-round(k)
if ( is.na(k) || k < 0 || max(abs(k - l)) > 1e-07) 
        stop("'k' must be nonnegative and integer")
m<-round(n)
if ( is.na(n) || n < k || max(abs(n - m)) > 1e-07) 
        stop("'n' must be nonnegative and integer >= k")
if(alpha<0 || alpha>1) {
stop("'alpha' must be a number between 0 and 1")}
if(k==0||CI=="upper")
{ll<-0
 ul<-qbeta(1-alpha,k+1,n-k)}
else if (k==n||CI=="lower") 
{ul<-1
 ll<-qbeta(alpha,k,n-k+1)}
else if (CI=="two.sided")
{ll<-qbeta(alpha/2,k,n-k+1)
 ul<-qbeta(1-alpha/2,k+1,n-k)}
else 
stop("undefined CI detected")
data.frame(Confidence.Interval=CI,Lower.limit=ll,Upper.limit=ul,alpha=alpha,row.names="")}
