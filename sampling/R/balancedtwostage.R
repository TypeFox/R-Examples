"balancedtwostage" <-
function(X,selection,m,n,PU,comment=TRUE,method=1)
{
N=dim(X)[1]
p=dim(X)[2]
str=cleanstrata(PU)
M=max(PU)
res1=balancedcluster(X,m,PU,method,comment)
if(selection==2) 
        {
         pik2=rep(n/N*M/m,times=N);
         if(n/N*M/m>1) stop("at the second stage, inclusion probabilities larger than 1");
         }
if(selection==1) 
        {
        pik2=inclusionprobastrata(str,rep(n/m ,times=max(str)));
        if(max(pik2)>1) stop("at the second stage, inclusion probabilities larger than 1");
        }
liste=(res1[,1]==1)
sf=rep(0,times=N)
sf[liste]=balancedstratification(array(X[liste,]/res1[,2][liste],c(sum(as.integer(liste)),p)),cleanstrata(str[liste]),pik2[liste],comment,method)
x=cbind(sf,res1[,2]*pik2,res1[,1],res1[,2],pik2)
colnames(x)=c("second_stage","final_pik", "primary","pik_first_stage", "pik_second_stage")
x 
}


