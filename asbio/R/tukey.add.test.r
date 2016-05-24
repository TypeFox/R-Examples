tukey.add.test<-function(y,A,B){
dname<-paste(deparse(substitute(A)), "and",deparse(substitute(B)), "on ",deparse(substitute(y)))
A<-factor(A);B<-factor(B)
ybar..<-mean(y)
ybari.<-tapply(y,A,mean)
ybar.j<-tapply(y,B,mean)
len.means<-c(length(levels(A)),length(levels(B)))
SSAB<-sum(rep(ybari.-ybar..,len.means[2])*rep(ybar.j-ybar.., rep(len.means[1],len.means[2]))*tapply(y,interaction(A,B),mean))^2/(sum((ybari.-ybar..)^2)*sum((ybar.j-ybar..)^2))
aovm<-anova(lm(y~A+B))
SSrem<-aovm[3,2]-SSAB
dfdenom<-aovm[3,1]-1
STATISTIC<-SSAB/SSrem*dfdenom
names(STATISTIC)<- "F"
DEGREES.OF.FREEDOM<-c(1,dfdenom)
names(DEGREES.OF.FREEDOM)<-c("num.df","denom.df")
D<-sqrt(SSAB/(sum((ybari.-ybar..)^2)*sum((ybar.j-ybar..)^2)))
names(D)="D"
RVAL<-list(statistic=STATISTIC,df=DEGREES.OF.FREEDOM,p.value=1-pf(STATISTIC,1,dfdenom),head="Tukey's one df test for additivity",data.name=dname)
attr(RVAL, "class")<-"addtest"
RVAL
}

print.addtest <- function(x, digits= max(3, getOption("digits")),...){
cat("\n")
cat(x$head,"\n")
p <- paste("F =", round(x$statistic[1], digits), " ", "Denom df =", x$df[2], "  ", "p-value =", round(x$p.value[1], digits))
cat(p)
cat("\n")
invisible(x)
}

