pr1 <-function(data, plateau=FALSE, x.plateau=NULL){
plat=ifelse(class(x.plateau)=="NULL", 1,2)
d1 = as.data.frame(data[,1])
    d2 = as.data.frame(data[,-1])
    f = function(h) {
        data.frame(d1, d2[h])
    }
    h = length(d2)
    h = 1:h
    l = lapply(h, f)

R2 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        return(r1)
    }
R3 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        p1 <- (gl/((gl + 1) - (length(coef(m) + 1))) * (1 - 
            r1))
        r2 <- 1 - p1; return(r2)
    }
reg=
function(data){
names(data)=c("x","y")
m1=lm(data[,2]~data[,1])
m2=lm(data[,2]~data[,1]+I(data[,1]^2))
c3=coef(m2)[[1]]
c4=coef(m2)[[2]]
c5=coef(m2)[[3]]
s1=summary(m1)
s2=summary(m2)
c11=s1[[4]][,1]
c22=s2[[4]][,1]
p1=s1[[4]][,4]
p2=s2[[4]][,4]
A1=AIC(m1)
A2=AIC(m2)
B1=BIC(m1)
B2=BIC(m2)
r1=R2(m1)
r2=R2(m2)
ar1=R3(m1)
ar2=R3(m2)
pm=c3-(c4^2)/(4*c5)
pc=-0.5*c4/c5
linear=c(c11,NA,p1,NA,r1,ar1,A1,B1,NA,NA); linear=round(linear,4)
quadratic=c(c22,p2,r2,ar2,A2,B2, pm, pc); quadratic=round(quadratic,4)
resp=data.frame(linear, quadratic)
rownames(resp)=c("coefficient a","coefficient b","coefficient c","p-value t.test for a","p-value t.test for b","p-value t.test for c", "r-squared","adjusted r-squared", "AIC", "BIC", "maximum or minimum value for y","critical point in x")
return(resp)
        }
funreg=
function(data){
names(data)=c("x","y")
ml=lm(data[,2]~data[,1])
mq=lm(data[,2]~data[,1]+I(data[,1]^2))
c1=coef(ml)[[1]]
c2=coef(ml)[[2]]
c3=coef(mq)[[1]]
c4=coef(mq)[[2]]
c5=coef(mq)[[3]]
m1<-nls(y~a+b*x, start=list(a=c1, b=c2), data=data)
m2<-nls(y~a+b*x+c*x^2, start=list(a=c3, b=c4, c=c5), data=data)
s1=summary(m1)
s2=summary(m2)
c11=s1[[11]][,1]
c22=s2[[11]][,1]
p1=s1[[11]][,4]
p2=s2[[11]][,4]
A1=AIC(m1)
A2=AIC(m2)
B1=BIC(m1)
B2=BIC(m2)
r1=R2(m1)
r2=R2(m2)
ar1=R3(m1)
ar2=R3(m2)
pm=c3-(c4^2)/(4*c5)
pc=-0.5*c4/c5
pc1=-0.5*c4/c5
pc=list(pc, x.plateau)
pc=pc[[plat]]
m3 <- nls(y~a+b*(x-pc)*(x<=pc), start=list(a=c1, b=c2, pc=pc), data=data, control = nls.control(maxiter = 6000))
m4=nls(y~(a+b*x+c*I(x^2))*(x<=-0.5*b/c)+(a+I(-b^2/(4*c)))*(x>-0.5*b/c),start=list(a=c3, b=c4, c=c5),data=data, control = nls.control(maxiter = 6000))
ss1=summary(m3)
ss2=summary(m4)
cc11=ss1[[11]][,1]
cc22=ss2[[11]][,1]
pp1=ss1[[11]][,4]
pp2=ss2[[11]][,4]
Aa1=AIC(m3)
Aa2=AIC(m4)
Bb1=BIC(m3)
Bb2=BIC(m4)
ra1=R2(m3)
ra2=R2(m4)
ara1=R3(m3)
ara2=R3(m4)
ccc3=round(cc22[1],4)
ccc4=round(cc22[2],4)
ccc5=round(cc22[3],4)
pmm=(ccc3+I(-ccc4^2/(4*ccc5)))
pcc=-0.5*ccc4/ccc5
linear=c(c11,NA,p1,NA,r1,ar1,A1,B1,NA,NA); linear=round(linear,4)
quadratic=c(c22,p2,r2,ar2,A2,B2, pm, pc1); quadratic=round(quadratic,4)
linear.plateau=c(cc11[-3],NA,pp1[-3],NA,ra1,ara1,Aa1,Bb1,cc11[1],cc11[3]); linear.plateau=round(linear.plateau,4)
quadratic.plateau=c(cc22,pp2,ra2,ara2,Aa2,Bb2, pmm, pcc); quadratic.plateau=round(quadratic.plateau,4)
resp=data.frame(linear, quadratic, linear.plateau, quadratic.plateau)
rownames(resp)=c("coefficient a","coefficient b","coefficient c","p-value t.test for a","p-value t.test for b","p-value t.test for c", "r-squared","adjusted r-squared", "AIC", "BIC", "maximum or minimum value for y","critical point in x")
return(resp)
        }
i=ifelse(plateau==FALSE,1,2)
lf=list(reg,funreg)
FUN=lf[[i]]
rep=lapply(l, FUN)
names(rep)= names(data)[-1]
return(rep)
}

