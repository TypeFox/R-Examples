`cd` <-
function(x){

DNAME = deparse(substitute(x))

##Strip NAs
x<-na.omit(x)


### Calculate the Coefficient of Dispersion

x=sort(x)
n=length(x)
M=median(x)
coef=mean(abs(x-M))/M

METHOD = "Measures of Variability - Coefficient of Dispersion"


### Display Output

STATISTIC=coef
names(STATISTIC)="Coefficient of Dispersion"

structure(list(statistic = STATISTIC, method = METHOD, data.name = DNAME), class = "htest")

}

