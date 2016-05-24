`sj.test` <-
function(x, crit.values=c("t.approximation", "empirical"), N=0) 
{
crit.values=match.arg(crit.values)


if ((crit.values=="empirical")&(N==0)) 
{stop("number of Monte Carlo simulations N should be provided for the empirical critical values")}

### SJ Test - New Directional Test

DNAME = deparse(substitute(x))

n<-length(x)
J<-sqrt(pi/2)*mean(abs(x-median(x))) 
x<-sort(x)
cw1<-sd(x)/J
statistic = sqrt(n)*(cw1-1)/sqrt((pi-3)/2)

if(crit.values=="empirical")
  {

#### computes empirical critical values for the SJ statistic####

   sj<-double(N)

   for (k in 1:N)
      {
       e <- rnorm(length(x), mean=0, sd = sqrt(1))  
       J<-sqrt(pi/2)*mean(abs(e-median(e))) 
       sj[k]<-sd(e)/J
      }
   y<-sort(sj)
   if (cw1>=max(y)) {p.value=0}
   else if (cw1<=min(y)) {p.value=1}
   else
      {
   bn<-which(y==min(y[I(y>=cw1)]))
   an<-which(y==max(y[I(y<cw1)]))
   a<-max(y[I(y<cw1)])
   b<-min(y[I(y>=cw1)])
   pa<-(an - 1) / (N - 1)
   pb<-(bn - 1) / (N - 1)
   alpha<-(cw1-a)/(b-a)  
   p.value=1-alpha*pb-(1-alpha)*pa
      }


 }

else if (crit.values=="t.approximation") {p.value=1-pt(statistic, df=(sqrt(n)+3)/2)}


METHOD = "Test of Normality - SJ Test"


### Display Output ###

STATISTIC=statistic
names(STATISTIC)="Standardized SJ Statistic"
PARAMETER = cw1
names(PARAMETER) = "ratio of S to J"

structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value=p.value, method = METHOD, data.name = DNAME), class = "htest")

}

