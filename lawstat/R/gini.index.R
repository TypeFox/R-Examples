`gini.index` <-
function(x){

DNAME = deparse(substitute(x))

##Strip NAs
x<-na.omit(x)

x=sort(x)
n=length(x)

### Calculate the delta and Gini Index ###

a<-0
for (i in 1:(n-1))
for (j in (i+1):n)
          {
           a<-a+abs(x[j]-x[i])
          }

delta<-2*a/n/(n-1)
GI<-delta/(2*mean(x))

METHOD = "Measures of Relative Variability - Gini Index"


### Display Output ###

STATISTIC=GI
names(STATISTIC)="Gini Index"
PARAMETER = delta
names(PARAMETER) = "delta"

structure(list(statistic = STATISTIC, parameter = PARAMETER, method = METHOD, data.name = DNAME), class = "htest")

}

