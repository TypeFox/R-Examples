itemanal <- function(object, SE.par=1.96, boots=1000){

if (is.null(row.names(t(object)))) {print("Dataset must be coerced into matrix from prior data frame.")}

else {

if (!is.numeric(object)) {print("Your dataset is not a numeric object.")} 

else {

object<-(na.omit(object))
yy<-sort(abs(as.numeric(var(object))))
yy<-yy[1]
yy<-ifelse(yy==0, NA, yy)

if (is.na(yy)) {print("One of your variables is a constant. Constants are disallowed as part of a scale.")}

else {

if ( dim(cbind(object))[2]< 2 ) {print("Your set contains only one variable. At least two are required.")}

else {

if ( dim(cbind(object))[1]< 2 ) {print("Your set contains only one case. You need at least two cases.")}

else {


same <- function(x) {
for (i in 1:ncol(x)) {
	for (j in 1:ncol(x)) {
		if (i!=j) {
			test <- x[,i]-x[,j]
			if (sum(abs(test))==0) {print(paste("WARNING: Items ",i," (",colnames(x)[i],") ","and ",j," (",colnames(x)[j],") ","are duplicates."))} 
			}
		}
	}
	test
}

notneed <- same(object)


object<-(na.omit(object))

options(digits=5)

Type<-mode(object)
Cases<-apply(object,2,length)
Minimum<-apply(object, 2, min)
Maximum<-apply(object, 2, max)
Mean<-apply(object, 2, mean)
Median<-apply(object, 2, median)

sum1<-sum((object-mean(object))^3)
N1<-(length(object))-1
SDF<- function(object) {
sqrt((sum((object-mean(object))^2))/N1)}
SD <- apply(object, 2, SDF)

Sum<-apply(object,2,sum)
SE.mean<-(SD)/(sqrt(length(object)))
Lower<-(Mean)-(SE.par*SD)
Upper<-(Mean)+(SE.par*SD)

Skew<- function(object){
(sum1/N1)*((SDF(object))^3)}

Skewness<-apply(object, 2, Skew)

er.skew <- sqrt(6/(dim(object)[1]))
lower.skew<- SE.par*er.skew
upper.skew<- -SE.par*er.skew

sum2<-sum((object-mean(object))^4)

Kurt<- function(object) {
(sum2/N1)*((SDF(object))^4)}

Kurtosis<-apply(object, 2, Kurt)

er.kurt <- sqrt(24/(dim(object)[1]))
lower.kurt<- SE.par*er.kurt
upper.kurt<- -SE.par*er.kurt

eivar<- apply(object, 2, var)

summary.statement1<-data.frame(cbind("Variable.type"=Type, "No.of.cases"=Cases, Minimum, Maximum, Sum))
summary.statement2<-cbind(Mean, Median, SD, SE.mean, Lower, Upper, "Variance"=eivar)
summary.statement3a<-cbind(Skewness,"SE.skew"=er.skew, "Lower"=lower.skew, "Upper"=upper.skew)
summary.statement3b<-cbind(Kurtosis,"SE.kurtosis"=er.kurt, "Lower"=lower.kurt, "Upper"=upper.kurt)

Covariances<-cov(object,object)

Correlations<-cor(object,object)

# Cronbach's Alpha

Columns<-ncol(object)
Rows<-nrow(object)
Alpha<-(Columns/(Columns-1))*(1-sum(apply(object,2,var))/(var(apply(object,1,sum))))

# Duhachek SE

numbitem = Columns
numbsubj = Rows
itemcov = cov(object)

one=matrix(1,numbitem,1)
jtphij = t(one)
jtphij=jtphij%*%itemcov
jtphij=jtphij%*%one 
trmy=sum(diag(itemcov)) 
trmy=trmy/jtphij 
myalpha=1-trmy 
nn1=numbitem-1 
nn1=numbitem/nn1 
myalpha<-nn1*myalpha

trphisq=itemcov%*%itemcov 
trphisq=sum(diag(trphisq))
trsqphi=sum(diag(itemcov))
trsqphi=trsqphi^2

ttp=itemcov%*%itemcov 
jtphisqj=t(one)
jtphisqj=jtphisqj%*%ttp
jtphisqj=jtphisqj%*%one 

omega=trphisq+trsqphi 
omega=jtphij*omega 
omegab=sum(diag(itemcov))
omegab=omegab*jtphisqj 
omega=omega-(2*omegab) 
omega=(2/(jtphij^3))*omega 

s2=(numbitem^2)/((numbitem-1)^2) 
s2=s2*omega

se<-sqrt(s2/numbsubj)

cimin95<-myalpha-(SE.par*se) 
cimax95<-myalpha+(SE.par*se) 

summary.statement4<-data.frame("No.of.items"=Columns, "Cronbachs.alpha"=Alpha)
summary.statement5<-data.frame("SE.for.alpha"=se, "Lower"=cimin95, "Upper"=cimax95)

# Duhachek SE standardized alpha

itemcov = cor(object)

one=matrix(1,numbitem,1)
jtphij = t(one)
jtphij=jtphij%*%itemcov
jtphij=jtphij%*%one 
trmy=sum(diag(itemcov)) 
trmy=trmy/jtphij 
myalpha=1-trmy 
nn1=numbitem-1 
nn1=numbitem/nn1 
myalpha<-nn1*myalpha

trphisq=itemcov%*%itemcov 
trphisq=sum(diag(trphisq))
trsqphi=sum(diag(itemcov))
trsqphi=trsqphi^2

ttp=itemcov%*%itemcov 
jtphisqj=t(one)
jtphisqj=jtphisqj%*%ttp
jtphisqj=jtphisqj%*%one 

omega=trphisq+trsqphi 
omega=jtphij*omega 
omegab=sum(diag(itemcov))
omegab=omegab*jtphisqj 
omega=omega-(2*omegab) 
omega=(2/(jtphij^3))*omega 

s2=(numbitem^2)/((numbitem-1)^2) 
s2=s2*omega

se<-sqrt(s2/numbsubj)

cimin95<-myalpha-(SE.par*se) 
cimax95<-myalpha+(SE.par*se) 

summary.statement6<-data.frame("No.of.items"=Columns, "Standardized.alpha"=myalpha)
summary.statement7<-data.frame("SE.for.standardized.alpha"=se, "Lower"=cimin95, "Upper"=cimax95)


if.deleted <- function(object) {

item.sum <- apply(object, 2, sum)
scale.sum <- sum(item.sum)
scale.mean <- (scale.sum)/(length(apply(object,1,length)))

scale.mean.var <- data.frame("Variable"=row.names(t(object)), scale.mean)

n.0 <- length(apply(object,2,length))
n.10 <- length(apply(object,1,length))

ave <- function(x) {
(scale.sum-x)/(length(apply(object,1,length)))}

mean.if.deleted <- c()
for (i in 1:length(item.sum)) {
	mean.if.deleted<-rbind(mean.if.deleted,ave(item.sum[i]))
}

deleted <- data.frame("Scale.mean"=scale.mean.var[,-1],"Mean.if.deleted"=mean.if.deleted[,1])
}

summary.statement8a <- if.deleted(object)


# Variance If Item Deleted


if.deleted.variance <- function(object) {

item.sum.variance <- apply(object, 1, sum)
item.variance <- var(item.sum.variance)

if (dim(object)[2]<3) {del.variance <- c(var((object)[,1]),var((object)[,2]))

item.variance <- c(item.variance, item.variance)
del.variance <- data.frame("Scale.variance"=item.variance, "Variance.if.deleted"=del.variance)
}

else {

del.variance <- c(1)
for (i in 1:dim(object)[2]) {del.variance<-(c(del.variance,(var(apply(object[,-1*i],1,sum)))))}
del.variance <- del.variance[-1]

del.variance <- data.frame("Scale.variance"=item.variance, "Variance.if.deleted"=del.variance)
}

}

summary.statement8b <- if.deleted.variance(object)

# Corrected total item correlation if item deleted

item.cor<- function(object) {

if (dim(object)[2]<3) {

cor1<- cor(object[,1], object[,2])
cor2<- cor(object[,2], object[,1])

del.cor<- c(cor1, cor2)
del.cor<-data.frame("Corrected.item.total.cor"=del.cor)
}

else

{
del.cor <- c(1)
for (i in 1:dim(object)[2]){ del.cor <- c(del.cor, cor((sapply(object[,i],sum)), (apply(object[,-1*i],1,sum))))}
del.cor <- del.cor[-1]
del.cor<- data.frame("Corrected.item.total.cor"=del.cor)
}
}

summary.statement8 <- data.frame(summary.statement8a, summary.statement8b)
summary.statement9a <- data.frame("Variable"=row.names(t(object)), item.cor(object))


# Square Multiple correlation If Item Deleted

r.if.del<- function(object) {

if (dim(object)[2]<3) {

del.lm<-lm(object[,1]~object[,2])
del.lm.sum1<-summary(del.lm)

del.lm<-lm(object[,2]~object[,1])
del.lm.sum2<-summary(del.lm)

r.if.item.del<-c(del.lm.sum1$r.squared, del.lm.sum2$r.squared)

}

else {

r.if.item.del<- c(1)
for (i in 1:dim(object)[2]) { 

r.if.item.del<- c(r.if.item.del, (summary(lm ( object[,i] ~ object[,-1*i])))$r.squared ) }

r.if.item.del <- r.if.item.del[-1]
}
}

summary.statement10<-data.frame("Squared.multiple.cor"=r.if.del(object))



# Alpha if item deleted

del.alpha <- function(object) {

if (dim(object)[2]<3) {make.alpha<-c("N/A","N/A")
}

else {

make.alpha <- c(1)

n.items <- dim(object)[2]

for (i in 1:n.items) { 

make.alpha <- c(make.alpha ,(((length(1:n.items)-1) / ((length(1:n.items))-2)) * (1-sum(apply(object[,-1*i],2,var)) / (var(apply(object[,-1*i],1,sum))))))
}

make.alpha <- make.alpha[-1]

}
}

summary.statement9b <- del.alpha(object)
summary.statement9 <- data.frame(summary.statement9a,summary.statement10, "Alpha.if.deleted"=summary.statement9b)


# Bootstrap based alpha standard error and confidence interval
boot.one <-
function(object, rep=boots, ci.par=SE.par, nd=4) {
new <- c(1)

func <- function(object) {
Columns<-ncol(object)
Rows<-nrow(object)
numbitem = Columns
numbsubj = Rows
itemcov = cov(object)
one=matrix(1,numbitem,1)
jtphij = t(one)
jtphij=jtphij%*%itemcov
jtphij=jtphij%*%one 
trmy=sum(diag(itemcov)) 
trmy=trmy/jtphij 
myalpha=1-trmy 
nn1=numbitem-1 
nn1=numbitem/nn1 
myalpha<-nn1*myalpha
myalpha }

{ repeat
{ord <- sample(1:length(object[,1]), length(object[,1]), replace=TRUE )
XX <- object[ord,]
XXz <- round(func(XX), digits=nd)
new <- c(new, XXz)
if(length(new)+.1>rep+1){break} }}

new <- new[-1]
orig <- round(func(object), digits=nd)
SE <- round(sd(new), digits=nd)
ci.lower <- round(orig-(ci.par*SE), digits=nd)
ci.upper <- round(orig+(ci.par*SE), digits=nd)
MM <- round(mean(new), digits=nd)
final <- list("new"=new,MM,SE, ci.lower, ci.upper)
final}

b.o.o <- boot.one(object)
b.o.o.sim <- as.numeric(b.o.o$new)

summary.statement10a <- data.frame(cbind(b.o.o[2],b.o.o[3],b.o.o[4],b.o.o[5]))
colnames(summary.statement10a) <- c("Mean", "SE", "Lower", "Upper")

# Bootstrap based standardized alpha standard error and confidence interval
boot.two <-
function(object, rep=boots, ci.par=SE.par, nd=4) {
new <- c(1)

func <- function(object) {
Columns<-ncol(object)
Rows<-nrow(object)
numbitem = Columns
numbsubj = Rows
itemcov = cor(object)
one=matrix(1,numbitem,1)
jtphij = t(one)
jtphij=jtphij%*%itemcov
jtphij=jtphij%*%one 
trmy=sum(diag(itemcov)) 
trmy=trmy/jtphij 
myalpha=1-trmy 
nn1=numbitem-1 
nn1=numbitem/nn1 
myalpha<-nn1*myalpha
myalpha }

{ repeat
{ord <- sample(1:length(object[,1]), length(object[,1]), replace=TRUE )
XX <- object[ord,]
XXz <- round(func(XX), digits=nd)
new <- c(new, XXz)
if(length(new)+.1>rep+1){break} }}

new <- new[-1]
orig <- round(func(object), digits=nd)
SE <- round(sd(new), digits=nd)
ci.lower <- round(orig-(ci.par*SE), digits=nd)
ci.upper <- round(orig+(ci.par*SE), digits=nd)
MM <- round(mean(new), digits=nd)
final <- list("new"=new,MM,SE, ci.lower, ci.upper)
final}

b.o.o2 <- boot.two(object)
b.o.o2.sim <- as.numeric(b.o.o2$new)

summary.statement11 <- data.frame(cbind(b.o.o2[2],b.o.o2[3],b.o.o2[4],b.o.o2[5]))
colnames(summary.statement11) <- c("Mean", "SE", "Lower", "Upper")

item.analysis.output<-list(Variables=summary.statement1, Tendency=summary.statement2,
Skewness=summary.statement3a, Kurtosis=summary.statement3b, Covariance=Covariances, Correlation=Correlations,
Alpha=summary.statement4, Conf.Alpha=summary.statement5,Bootstrap.Simmulations=b.o.o.sim,
Alpha.Bootstrap=summary.statement10a, Std.Alpha=summary.statement6, Conf.Std.Alpha=summary.statement7,
Bootstrap.Std.Simmulations=b.o.o2.sim,Std.Alpha.Bootstrap=summary.statement11,
Scale.Stats=summary.statement8, Alpha.Stats= summary.statement9)

ob <- as.character(match.call()[2])
cl <- call("itemanal","object"=ob,"SE.par"=SE.par, "boots"=boots)

item.analysis.output$call <- cl
item.analysis.output$items <- names(item.analysis.output)

class(item.analysis.output) <- c("itemanal", class(item.analysis.output))

item.analysis.output

}}}}}}