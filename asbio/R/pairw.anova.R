pairw.anova<-function(y,x,conf.level=0.95,method="tukey",MSE=NULL,df.err=NULL,control=NULL){
indices <- c("lsd","bonf","tukey","scheffe","dunnett")
method <- match.arg(method, indices)
if(method=="lsd"){
res <- lsdCI(y,x,conf.level=conf.level,MSE=MSE,df.err=df.err)}
if(method=="bonf"){
res <- bonfCI(y,x,conf.level=conf.level,MSE=MSE,df.err=df.err)}
if(method=="tukey"){
res <- tukeyCI(y,x,conf.level=conf.level,MSE=MSE,df.err=df.err)}
if(method=="scheffe"){
res <- scheffeCI(y,x,conf.level=conf.level,MSE=MSE,df.err=df.err)}
if(method=="dunnett"){
res <- dunnettCI(y,x,conf.level=conf.level,control=control)}
res
}


lsdCI<-function(y, x, conf.level = 0.95, MSE = NULL, df.err = NULL){
fitted<-tapply(y,factor(x),mean)
nis<-tapply(y,factor(x),length)
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)

dif.mat<-outer(fitted,fitted,"-")
diffs<-dif.mat[upper.tri(dif.mat)]

SE.diff.mat<-sqrt(MSE*outer(1/nis,1/nis,"+"))
SE.diff<-SE.diff.mat[upper.tri(SE.diff.mat)]

p.val<- 2*pt(abs(diffs/SE.diff),df.error,lower.tail=FALSE)
t<-qt(1-(1-conf.level)/2,df.error)
hwidths<-t*SE.diff
LSD<-round(qt(1-((1-conf.level)/2),df.error)*SE.diff,5)
Decision<- ifelse(abs(diffs)>LSD,"Reject H0","FTR H0")
val<-as.data.frame(cbind(LSD,round(diffs, 5),round(diffs-hwidths,5),round(diffs+hwidths,5),Decision,round(p.val,5)))
lvl<-outer(levels(x),levels(x),function(x1,x2){paste(paste("mu",x1,sep=""),paste("mu",x2,sep=""),sep="-")})
dimnames(val)<-list(lvl[upper.tri(lvl)],
c("LSD","Diff","Lower","Upper","Decision","Adj. p-value"))
head<-paste(paste(as.character(conf.level*100),"%",sep=""),c("LSD confidence intervals"))
###
res <- list()
res$head <- head
res$conf <- conf.level
comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
res$comp <- comp[upper.tri(comp)]
res$summary <- val
res$band <- cbind(diffs-hwidths, diffs+hwidths)
res$fitted <- fitted
res$x <- x
res$y <- y
res$method <- "Fisher LSD"
class(res)<-"pairw"
res
}
 

bonfCI<-function(y,x, conf.level=0.95,MSE=NULL,df.err=NULL){
fitted<-tapply(y,factor(x),mean)
nis<-tapply(y,factor(x),length)
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)
r<-length(fitted)
dif.mat<-outer(fitted,fitted,"-")
diffs<-dif.mat[upper.tri(dif.mat)]
SE.diff.mat<-sqrt(MSE*outer(1/nis,1/nis,"+"))
L.hat.SE<-SE.diff.mat[upper.tri(SE.diff.mat)]
p.val<- 2*pt(abs(diffs/L.hat.SE),df=df.error,lower.tail=FALSE)
p.adj<-ifelse(p.val*((r^2-r)/2)>=1,1,round(p.val*((r^2-r)/2),6))
t<-qt(1-(1-conf.level)/(r*(r-1)),df.error)
hwidths<-t*L.hat.SE
val<-round(cbind(diffs,diffs-hwidths,diffs+hwidths),5)
Decision<-ifelse((val[,2]>0&val[,3]>0)|val[,2]<0&val[,3]<0,"Reject H0","FTR H0")
val<-as.data.frame(cbind(val,Decision,p.adj))
lvl<-outer(levels(x),levels(x),function(x1,x2){paste(paste("mu",x1,sep=""),paste("mu",x2,sep=""),sep="-")})
dimnames(val)<-list(lvl[upper.tri(lvl)],
c("Diff","Lower","Upper","Decision","Adj. p-value"))
head<-paste(paste(as.character(conf.level*100),"%",sep=""),c("Bonferroni confidence intervals"))
###
res <- list()
res$head <- head
res$conf <- conf.level
comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
res$comp <- comp[upper.tri(comp)]
res$summary <- val
res$band <- cbind(diffs-hwidths, diffs+hwidths)
res$fitted <- fitted
res$x <- x
res$y <- y
res$method <- "Bonferroni"
class(res)<-"pairw"
res
}


tukeyCI<-function(y,x, conf.level=0.95,MSE=NULL,df.err=NULL){
fitted<-tapply(y,factor(x),mean)
nis<-tapply(y,factor(x),length)
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)
r<-length(fitted)

dif.mat<-outer(fitted,fitted,"-")
diffs<-dif.mat[upper.tri(dif.mat)]

SE.diff.mat<-sqrt(MSE*outer(1/nis,1/nis,"+"))
SE.diff<-SE.diff.mat[upper.tri(SE.diff.mat)]

Q.star<-abs((sqrt(2)*diffs)/ SE.diff) 
p.val<-round(ptukey(Q.star,r,df=df.error,lower.tail=FALSE),6)
T<-qtukey(conf.level,r,df.error)/sqrt(2)
hwidths<-T*SE.diff
val<-round(cbind(diffs,diffs-hwidths,diffs+hwidths),5)
Decision<-ifelse((val[,2]>0&val[,3]>0)|val[,2]<0&val[,3]<0,"Reject H0","FTR H0")
val<-as.data.frame(cbind(val,Decision,p.val))
lvl<-outer(levels(x),levels(x),function(x1,x2){paste(paste("mu",x1,sep=""),paste("mu",x2,sep=""),sep="-")})
dimnames(val)<-list(lvl[upper.tri(lvl)],
c("Diff","Lower","Upper","Decision","Adj. p-value"))
head<-paste(paste(as.character(conf.level*100),"%",sep=""),c("Tukey-Kramer confidence intervals"))
###
res <- list()
res$head <- head
res$conf <- conf.level
comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
res$comp <- comp[upper.tri(comp)]
res$summary <- val
res$band <- cbind(diffs-hwidths, diffs+hwidths)
res$fitted <- fitted
res$x <- x
res$y <- y
res$method <- "Tukey HSD" 
class(res)<-"pairw"
res
}

 
scheffeCI<-function(y,x, conf.level=0.95,MSE=NULL,df.err=NULL){
fitted<-tapply(y,factor(x),mean)
nis<-tapply(y,factor(x),length)
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)
r<-length(fitted)
dif.mat<-outer(fitted,fitted,"-")
diffs<-dif.mat[upper.tri(dif.mat)]
SE.diff.mat<-sqrt(MSE*outer(1/nis,1/nis,"+"))
SE.diff<-SE.diff.mat[upper.tri(SE.diff.mat)]
F.star<-diffs^2/((r-1)*SE.diff^2)
p.val<-round(pf(F.star,r-1,df.error,lower.tail=FALSE),6)
S<-sqrt((r-1)*qf(conf.level,r-1,df.error))
hwidths<-S*SE.diff
val<-round(cbind(diffs,diffs-hwidths,diffs+hwidths),5)
Decision<-ifelse((val[,2]>0&val[,3]>0)|val[,2]<0&val[,3]<0,"Reject H0","FTR H0")
val<-as.data.frame(cbind(val,Decision,p.val))
lvl<-outer(levels(x),levels(x),function(x1,x2){paste(paste("mu",x1,sep=""),paste("mu",x2,sep=""),sep="-")})
dimnames(val)<-list(lvl[upper.tri(lvl)],
c("Diff","Lower","Upper","Decision","Adj. P-value"))
head<-paste(paste(as.character(conf.level*100),"%",sep=""),c("Scheffe confidence intervals"))
###
res <- list()
res$head <- head
res$conf <- conf.level
comp <- outer(levels(x),levels(x),function(x1,x2){paste(x1, x2, sep="-")})
res$comp <- comp[upper.tri(comp)]
res$summary <- val
res$band <- cbind(diffs-hwidths, diffs+hwidths)
res$fitted <- fitted
res$x <- x
res$y <- y
res$method <- "Scheffe"
class(res)<-"pairw"
res
}

dunnettCI<-function(y, x, conf.level = 0.95, control = NULL){
x <- factor(x)
if(is.null(control)) stop(call. = FALSE, "Please specify the control")
nis <- tapply(y, x, length)
nu <- sum(nis) - nlevels(x)
    
    p.sd <- function(y, x){
        rs <- seq(1, nlevels(x))
        for(i in levels(x)) rs[i] <- sum((y[x == i] - mean(y[x == i]))^2)
        rs <- rs[-c(1:nlevels(x))]
        sqrt(sum(rs)/nu)
    }
                                                  
s <- p.sd(y, x)
controlm <- mean(y[x == control])
fittedm <- tapply(y[x != control], factor(x)[x != control], mean)
fittedm <-fittedm[complete.cases(fittedm)]
controln <- length(y[x == control])
fittedn <- tapply(y[x != control], factor(x)[x != control], length)
fittedn <-fittedn[complete.cases(fittedn)]
diffs <- fittedm - controlm

    Dj <- seq(1, (nlevels(x) - 1))
    for(i in 1 : (nlevels(x)-1)) Dj[i] <- diffs[i]/(s*sqrt((1/fittedn[i]) + (1/controln)))
    
    Rij <- seq(1, (nlevels(x) - 1))
    for(i in 1 : (nlevels(x)-1)) Rij[i] <- sqrt(fittedn[i]/(fittedn[i] + controln))
    R <- outer(Rij, Rij, "*")
    diag(R) <- rep(1, (nlevels(x) - 1))
    
    upper <- seq(1, (nlevels(x) - 1))
    for(i in 1 : (nlevels(x)-1)) upper[i] <-  diffs[i] + s*sqrt((1/fittedn[i])+ (1/controln))*qmvt((1 - (1 - conf.level)/2), df = nu, sigma = R, tail="lower.tail")$quantile
    lower <- seq(1, (nlevels(x) - 1))
    for(i in 1 : (nlevels(x)-1)) lower[i] <- diffs[i] - s*sqrt((1/fittedn[i])+ (1/controln))*qmvt((1 - (1 - conf.level)/2), df = nu, sigma = R, tail="lower.tail")$quantile

fl <- levels(x)[levels(x)!=control]
con <- levels(x)[levels(x)==control]
decision<-ifelse((lower>0&upper>0)|lower<0&upper<0,"Reject H0","FTR H0")
lvl <- paste(paste("mu", fl, sep=""),paste("mu", con, sep=""),sep="-")
val <- as.data.frame(cbind(round(apply(cbind(lower, upper), 1, mean), 6),round(lower, 6),round(upper, 6),decision))
dimnames(val)<-list(lvl,c("Diff","Lower","Upper","Decision"))
head<-paste(paste(as.character(conf.level*100),"%",sep=""),c("Dunnett confidence intervals"))
###
res <- list()
res$head <- head
res$conf <- conf.level
comp <- paste(fl, con, sep="-")
res$comp <- comp
res$summary <- val
res$band <- cbind(lower, upper)
res$fitted <- fittedm
res$x <- x
res$y <- y
method <- "Dunnett"
class(res)<-"pairw"
res
}


scheffe.cont <- function(y, x, lvl = c("x1", "x2"), conf.level = 0.95, MSE = NULL, df.err = NULL){
alpha <- 1 - conf.level
n <- length(y); x <- factor(x)
r <- nlevels(x)
Y.new <- y[x == lvl[1]|x ==lvl[2]]
X.new <- x[x == lvl[1]|x ==lvl[2]]
means <- tapply(Y.new, X.new, mean)
means <- means[complete.cases(means)]
ni <- tapply(Y.new, X.new, length)
ni <- ni[complete.cases(ni)]
D <- means[1] - means[2]
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)
S.D <- sqrt(MSE * ((1/ni[1])+(1/ni[2])))
S <- sqrt((r - 1)*qf(1 - alpha, r - 1, n - r))
CI <- list(lower = D - S.D * S, upper = D + S.D * S)
F.star <- D^2/((r - 1)*(S.D^2))
P.val <- pf(F.star, r - 1, n - r, lower.tail = FALSE)
Dec <- ifelse(P.val <= alpha, "Reject H0", "FTR H0")
res <- data.frame(Diff = D, Lower = CI$lower, Upper = CI$upper, Decision = Dec, P = P.val)
names(res) <- c("Diff", "Lower", "Upper", "Decision", "Adj. P-value")
row.names(res) <- paste("mu",lvl[1],"-","mu",lvl[2], sep = "")
res
}


bonf.cont <- function(y, x, lvl = c("x1", "x2"), conf.level = 0.95, MSE = NULL, df.err = NULL, comps = 1){
alpha <- 1 - conf.level
n <- length(y); x <- factor(x)
Y.new <- y[x == lvl[1]|x ==lvl[2]]
X.new <- x[x == lvl[1]|x ==lvl[2]]
means <- tapply(Y.new, X.new, mean)
means <- means[complete.cases(means)]
ni <- tapply(Y.new, X.new, length)
ni <- ni[complete.cases(ni)]
D <- means[1] - means[2]
a<-anova(lm(y~factor(x)))
df.error<-ifelse(is.null(df.err),a$Df[length(a$Df)],df.err)
MSE<-ifelse(is.null(MSE),a$"Mean Sq"[length(a$"Mean Sq")],MSE)
S.D <- sqrt(MSE * ((1/ni[1])+(1/ni[2])))
t.crit <- qt(1-(alpha/(2*comps)), df.error)
CI <- list(lower = D - S.D * t.crit, upper = D + S.D * t.crit)
P.val <- 2*pt(abs(D/S.D), df.error, lower.tail = FALSE)*comps
Dec <- ifelse(P.val <= alpha, "Reject H0", "FTR H0")
res <- data.frame(Diff = D, Lower = CI$lower, Upper = CI$upper, Decision = Dec, P = ifelse(P.val>1,1,P.val))
names(res) <- c("Diff", "Lower", "Upper", "Decision", "Adj. P-value")
row.names(res) <- paste("mu",lvl[1],"-","mu",lvl[2], sep = "")
res
}

print.pairw <- function(x, ...){
cat("\n")
cat(x$head,"\n")
cat("\n")
rq<-structure(x$summary)
print(rq)
invisible(x)
}

plot.pairw <- function(x, type = 1, lcol = 1, lty = NULL, lwd = NULL, cap.length = 0.1, xlab = "", main = NULL, ...){
if(class(x)!="pairw") stop("Requires object of class pairw")
if(type == 1){
levels <- factor(names(x$fitted))
cont <- outer(levels, levels, function(x1,x2)paste(x1,x2,sep="-"))
cont1 <- cont[upper.tri(cont)]
dec <- x$summary$Decision
dec1 <- ifelse(dec=="FTR H0", FALSE, TRUE)
xlvl <- x$x
y <- x$y
names(dec1) <- cont1
lett <- multcompLetters(dec1)$Letters #Requires multcompView
bplot(y, xlvl, simlett = TRUE, lett = lett, xlab = xlab, ...)
if(x$method == "kBonferroni"){
cat(paste("The true effects of factor levels with the same letter are not", "\n", "significantly different at alpha = ", 1- x$conf, " using the ", strsplit(x$method, "k")[[1]][2], " method.", "\n\n", sep = ""))
} else if (x$method == "mBonferroni"){
cat(paste("The true effects of factor levels in blocks with the same letter are not", "\n", "significantly different at alpha = ", 1- x$conf, " using the ", strsplit(x$method, "m")[[1]][2], " method.", "\n\n", sep = ""))
} else{
cat(paste("The population means of factor levels with the same letter are not", "\n", "significantly different at alpha = ", 1- x$conf, " using the ", x$method, " method.", "\n\n", sep = ""))
}
}
if(type == 2){
ci <- x$band
w <- diff(range(min(ci[,1]),max(ci[,2])))*.05
r <- nrow(x$summary)
b <- barplot(rep(0, r), names.arg = x$comp, main = ifelse(is.null(main), x$head, main), xlab = xlab, horiz = TRUE, xlim = c(min(ci[,1]) - w , max(ci[,2])), border=NA, axis.lty = 1, ...)
abline(h = b, col = "gray")
arrows(ci[,1], b, ci[,2], b, angle = 90, col = lcol, lty = lty, lwd = lwd, length = cap.length, code = 3)
abline(v=0, col="gray", lty=2)
points(apply(x$band,1,mean), b, pch = 19)
}
}