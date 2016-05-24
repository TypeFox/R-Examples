regplot <-
function(data, xlab=NULL, ylab=NULL, poly=1, position=6, colors=TRUE, mean=TRUE, variable=1, x.plateau=NULL){
variable=variable+1
plat=ifelse(class(x.plateau)=="NULL", 1,2)
data=data[,c(1,variable)]
dreg=data   
means=function(data){
t=as.factor(data[,1])
d=data.frame(t,data[,-1])
s=split(data.frame(d[,-1]), d$t)
r=lapply(s, colMeans, na.rm=TRUE)
r=lapply(r, round,2)
rr=t(data.frame(r)); rr=data.frame(rr);rownames(rr)=NULL
treat=levels(t)
rr=data.frame(treat, rr); colnames(rr)=colnames(data)
return(rr)
}
n=names(data)
	y=n[2]
        x=n[1]
x=ifelse(class(xlab)=="NULL", x,xlab)
y=ifelse(class(ylab)=="NULL", y,ylab)

     rrr=means(data)
oi=ifelse(mean==TRUE,2,1)
l=list(data, rrr)
data=l[[oi]]

f1=function(dreg){
ml=lm(dreg[,2]~dreg[,1])
rsl=summary(ml)[[8]]
fl=function(qi){coef(ml)[[1]]+coef(ml)[[2]]*qi}
c1=round(coef(ml)[[1]],4)
c2=round(coef(ml)[[2]],4)
r2=rsl
r2=round(r2,2)
sin1=ifelse(c2>0,"+","")
e1=substitute(y==c1*sin1*c2*x*"  "* R^2*" = "*r2,  list(c1 = c1, c2 = c2,r2=r2, sin1=sin1))
c2=ifelse(colors==TRUE, 1,2)
        cor=list(2,1)
        cor=cor[[c2]]
        t=list("top", "bottomright", "bottom", "bottomleft", "left", "topleft", "topright", "right", "center")
        position=position
        p=t[[position]]
minx=min(dreg[,1])-sd(dreg[,1])/2
maxx=max(dreg[,1])+sd(dreg[,1])/2
miny=min(dreg[,2])-sd(dreg[,2])/2
maxy=max(dreg[,2])+sd(dreg[,2])/2
        plot(as.numeric(as.character(data[,1])),data[,2], xlab = x, 
              ylab = y, col = 1, pch=20, xlim=c(minx,maxx), ylim=c(miny,maxy), main="linear")
	curve(fl, min(dreg[,1]), max(dreg[,1]), add=TRUE, col = cor, lty = 1)
        legend(p, legend=e1, bty="n")
}

f2=function(dreg){
mq=lm(dreg[,2]~dreg[,1]+I(dreg[,1]^2))
rsq=summary(mq)[[8]]
fq=function(qi){coef(mq)[[1]]+coef(mq)[[2]]*qi+coef(mq)[[3]]*qi^2}
c3=round(coef(mq)[[1]],4)
c4=round(coef(mq)[[2]],4)
c5=round(coef(mq)[[3]],4)
r2=rsq
r2=round(r2,2)
sin1=ifelse(c4>0,"+","")
sin2=ifelse(c5>0,"+","")
e2=substitute(y==c3*sin1*c4*x*sin2*c5*x^2*"  "* R^2*" = "*r2,  list(c3 = c3, c4 = c4,c5=c5,r2=r2, sin1=sin1, sin2=sin2))
c2=ifelse(colors==TRUE, 1,2)
        cor=list(2,1)
        cor=cor[[c2]]
        t=list("top", "bottomright", "bottom", "bottomleft", "left", "topleft", "topright", "right", "center")
        position=position
        p=t[[position]]
minx=min(dreg[,1])-sd(dreg[,1])/2
maxx=max(dreg[,1])+sd(dreg[,1])/2
miny=min(dreg[,2])-sd(dreg[,2])/2
maxy=max(dreg[,2])+sd(dreg[,2])/2
        plot(as.numeric(as.character(data[,1])),data[,2], xlab = x, 
              ylab = y, col = 1, pch=20,xlim=c(minx,maxx), ylim=c(miny,maxy), main="quadratic")
	curve(fq, min(dreg[,1]), max(dreg[,1]), add=TRUE, col = cor, lty = 1)
        legend(p, legend=e2, bty="n")
}

R2 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        return(r1)
    }

f3=function(dreg){
ml=lm(dreg[,2]~dreg[,1])
c1=round(coef(ml)[[1]],4)
c2=round(coef(ml)[[2]],4)
mq=lm(dreg[,2]~dreg[,1]+I(dreg[,1]^2))
fq=function(qi){coef(mq)[[1]]+coef(mq)[[2]]*qi+coef(mq)[[3]]*qi^2}
c3=round(coef(mq)[[1]],4)
c4=round(coef(mq)[[2]],4)
c5=round(coef(mq)[[3]],4)
pm=c3-(c4^2)/(4*c5)
pc=-0.5*c4/c5
pc=list(pc, x.plateau)
pc=pc[[plat]]
names(dreg)=c("x","y")
m3=nls(y~a+b*(x-pc)*(x<=pc), start=list(a=c1, b=c2, pc=pc), data=dreg, control = nls.control(maxiter = 6000))
ss1=summary(m3)
cc11=ss1[[11]][,1]
c3=round(cc11[1],4)
c4=round(cc11[2],4)
cp=round(cc11[3],4)
fp1=function(x){c3+c4*(x-cp)*(x<=cp)}
r2=R2(m3)
r2=round(r2,2)
pc=round(pc,2)
pm=round(pm,2)
cp=round(cp,2)
sin1=ifelse(c4>0,"+","")
sin2=ifelse(cp>0,"-","")
e2=substitute(y==c3*sin1*c4*(x*sin2*cp)*(x<=cp) * "  x.plateau = "*cp*"  "* R^2*" = "*r2,  list(c3 = c3, c4 = c4,r2=r2, cp=cp, sin1=sin1, sin2=sin2))
c2=ifelse(colors==TRUE, 1,2)
        cor=list(2,1)
        cor=cor[[c2]]
        t=list("top", "bottomright", "bottom", "bottomleft", "left", "topleft", "topright", "right", "center")
        position=position
        p=t[[position]]
minx=min(dreg[,1])-sd(dreg[,1])/2
maxx=max(dreg[,1])+sd(dreg[,1])/2
miny=min(dreg[,2])-sd(dreg[,2])/2
maxy=max(dreg[,2])+sd(dreg[,2])/2
        plot(as.numeric(as.character(data[,1])),data[,2], xlab = x, 
              ylab = y, col = 1, pch=20, xlim=c(minx,maxx), ylim=c(miny,maxy), main="linear.plateau")
	curve(fp1, min(dreg[,1]), max(dreg[,1]), add=TRUE, col = cor, lty = 1)
        legend(p, legend=e2, bty="n")
}

f4=function(dreg){
ml=lm(dreg[,2]~dreg[,1])
c1=round(coef(ml)[[1]],4)
c2=round(coef(ml)[[2]],4)
mq=lm(dreg[,2]~dreg[,1]+I(dreg[,1]^2))
fq=function(qi){coef(mq)[[1]]+coef(mq)[[2]]*qi+coef(mq)[[3]]*qi^2}
c3=round(coef(mq)[[1]],4)
c4=round(coef(mq)[[2]],4)
c5=round(coef(mq)[[3]],4)
pm=(c3+I(-c4^2/(4*c5)))
pc=-0.5*c4/c5
names(dreg)=c("x","y")
m4=nls(y~(a+b*x+c*I(x^2))*(x<=-0.5*b/c)+(a+I(-b^2/(4*c)))*(x>-0.5*b/c),start=list(a=c3, b=c4, c=c5),data=dreg, control = nls.control(maxiter = 6000))
ss1=summary(m4)
cc11=ss1[[11]][,1]
c3=round(cc11[1],4)
c4=round(cc11[2],4)
c5=round(cc11[3],4)
pm=(c3+I(-c4^2/(4*c5)))
pc=-0.5*c4/c5
fp1=function(x){(c3+c4*x+c5*I(x^2))*(x<=-0.5*c4/c5)+(c3+I(-c4^2/(4*c5)))*(x>-0.5*c4/c5)}
r2=R2(m4)
r2=round(r2,2)
pc=round(pc,2)
pm=round(pm,2)
sin1=ifelse(c4>0,"+","")
sin2=ifelse(c5>0,"+","")
e2=substitute(y==c3*sin1*c4*x*sin2*c5*x^2*(x<=pc)* "  x.plateau = "*pc*"  "* R^2*" = "*r2,  list(c3 = c3, c4 = c4, c5=c5,r2=r2, pc=pc, pm=pm, sin1=sin1, sin2=sin2))
c2=ifelse(colors==TRUE, 1,2)
        cor=list(2,1)
        cor=cor[[c2]]
        t=list("top", "bottomright", "bottom", "bottomleft", "left", "topleft", "topright", "right", "center")
        position=position
        p=t[[position]]
minx=min(dreg[,1])-sd(dreg[,1])/2
maxx=max(dreg[,1])+sd(dreg[,1])/2
miny=min(dreg[,2])-sd(dreg[,2])/2
maxy=max(dreg[,2])+sd(dreg[,2])/2
        plot(as.numeric(as.character(data[,1])),data[,2], xlab = x, 
              ylab = y, col = 1, pch=20, xlim=c(minx,maxx), ylim=c(miny,maxy), main="quadratic.plateau")
	curve(fp1, min(dreg[,1]), max(dreg[,1]), add=TRUE, col = cor, lty = 1)
        legend(p, legend=e2, bty="n")
}
FUN=list(f1,f2, f3, f4)
FUN[[poly]](dreg)
   
    }
