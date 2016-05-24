`bartels.test` <-
function(y, alternative=c("two.sided", "positive.correlated", "negative.correlated"))
{
alternative<-match.arg(alternative)
DNAME = deparse(substitute(y))

##Strip NAs
y<-na.omit(y)


### Calculate the rank of the input data ###
R<-rank(y)
sum<-0
T<-length(y)
for (i in 1:(T-1)) {sum<-sum+(R[i]-R[i+1])^2}
num<-sum
sum<-0
for (i in 1:T){sum<-sum+(R[i]-(T+1)/2)^2}
den<-sum

#### Ratio of Von Neumann ####
RVN<-num/den
#statistic<-((RVN-2)*sqrt(5*T*(T+1)*(T-1)^2))/sqrt(4*(T-2)*(5*T^2-2*T-9))
statistic<-(RVN-2)*sqrt(T)/2


### One-sided or Two-sided Test ###
### Users will select the test alternative. Two-sided test is default ###

if(alternative=="positive.correlated"){
p.value = pnorm(statistic);
METHOD = "Bartels Test - Positive Correlated"}

else if(alternative=="negative.correlated"){
p.value = 1-pnorm(statistic);
METHOD = "Bartels Test - Negative Correlated"}

else {p.value=2*min(pnorm(statistic), 1-pnorm(statistic)); 
alternative="two.sided";
METHOD = "Bartels Test - Two sided"} 

### Display Output ###


PARAMETER=RVN
names(PARAMETER)="RVN Ratio"
STATISTIC=statistic
names(STATISTIC)="Standardized Bartels Statistic"

structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value=p.value, method = METHOD, data.name = DNAME), class = "htest")

}

