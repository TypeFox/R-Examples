`correlation` <-
function (x, y = NULL, method = c("pearson","kendall","spearman"), alternative = "two.sided")
{
x1<-x
y1<-y
    method <- match.arg(method)
    if (is.data.frame(y)) y1 <- as.matrix(y)
    if (is.data.frame(x)) x1 <- as.matrix(x)

if (!is.null(y1)) {
   if (!is.matrix(x1) &  !is.matrix(y1)) {
#-------------
# cor.vector
name.xy<-paste(deparse(substitute(x)), "and", deparse(substitute(y)))
xx<-cbind(x,y)
yy<-na.omit(xx)
nn<-length(yy[,1])
x<-yy[,1]
y<-yy[,2]
corr<-correl(x,y,method=method,alternative=alternative)
stat<-corr$stat
coef<-corr$rho
pvalue<-corr$pvalue
if(method=="pearson") {
gl<-nn-2
cat("\nPearson's product-moment correlation\n\n")
cat("data:",name.xy,"\n")
cat("t =",stat,", df =",gl,", p-value =",pvalue,"\n")
cat("alternative hypothesis: true rho is ")
if(alternative == "two.sided" ) cat("not equal to 0")
if(alternative == "less" ) cat("less than 0")
if(alternative == "greater") cat("greater than 0")
cat("\nsample estimates:\ncor\n",coef,"\n")
#list(t=t,df=gl,p.value=pvalue,rho=coef)
}
if(method=="spearman"){
cat("\nSpearman's rank correlation rho\n\n")
cat("data:",name.xy,"\n")
cat("p-value =",pvalue,"\n")
cat("alternative hypothesis: true rho is ")
if(alternative == "two.sided" ) cat("not equal to 0")
if(alternative == "less" ) cat("less than 0")
if(alternative == "greater") cat("greater than 0")
cat("\nsample estimates:\nrho\n",coef,"\n")
#list(S=t,p.value=pvalue,rho=coef)
}
if(method=="kendall"){
cat("\nKendall's rank correlation tau\n\n")
cat("data:",name.xy,"\n")
cat("z-norm = ",stat,"p-value =",pvalue,"\n")
cat("alternative hypothesis: true rho is ")
if(alternative == "two.sided" ) cat("not equal to 0")
if(alternative == "less" ) cat("less than 0")
if(alternative == "greater") cat("greater than 0")
cat("\nsample estimates:\ntau\n",coef,"\n")
#list(z.norm=stat,p.value=pvalue,tau=coef)
}
#-------------
}
   else {
# realizar cor.mv
#-------------
     if (is.data.frame(x) | is.matrix(x)) {
        nvarx <- ncol(x)
        if (is.matrix(x)) x<-data.frame(x)
        nombrex <- names(x)
        if (is.matrix(x)) nombrex<-colnames(x)
        x <- as.matrix(x)
    }
    else {
        nvarx <- 1
        nombrex<- deparse(substitute(x))
        x <- as.matrix(x)
    }
    if (is.data.frame(y) | is.matrix(y)) {
        nvary <- ncol(y)
        nombrey <- names(y)
        if (is.matrix(y)) nombrey<-colnames(y)
        y <- as.matrix(y)
    }
    else {
        nvary <- 1
        nombrey<- deparse(substitute(y))
        y <- as.matrix(y)
    }
    estimate <- rep(0, nvarx*nvary)
    dim(estimate) <- c(nvarx, nvary)
    dimnames(estimate) <- list(nombrex, nombrey)
    pvalue <- estimate
    nn <- round(estimate, 0)
    for (i in 1:nvarx) {
        for (j in 1:nvary) {
            xx <- cbind(x[, i], y[, j])
            yy <- na.omit(xx)
            nn[i, j] <- length(yy[, 1])
            xa <- yy[, 1]
            yb <- yy[, 2]
            corr <- correl(xa, yb, method = method, alternative=alternative)
            estimate[i, j] <- corr$rho
            pvalue[i, j] <- corr$pvalue
         }
    }
    names(method) = ""
    estimate <- round(estimate, 2)
    pvalue <- round(pvalue, 4)
    n1<-unique(c(nn))
    if(length(n1)==1)nn<-n1
cat("\nCorrelation Analysis\n\nMethod     :",method)
cat("\nAlternative:",alternative,"\n\n")
    lista <- list(correlation = estimate, pvalue = pvalue,
        n.obs = nn)
    return(lista)
#-------------
}
}
else {
#-------------
# "cor.matrix"
nvar<-ncol(x)
estimate<-rep(0,nvar*nvar)
nombres<-names(x)
if (is.matrix(x)) nombres<-colnames(x)
dim(estimate)<-c(nvar,nvar)
dimnames(estimate)<-list(nombres,nombres)
pvalue<-estimate
nn <- round(estimate,0)
x <-as.matrix(x)
for(i in 1:nvar){
for(j in 1:nvar){
xx<-cbind(x[,i],x[,j])
yy<-na.omit(xx)
nn[i,j]<-length(yy[,1])
x0<-yy[,1]
y0<-yy[,2]
if (i==j) {pvalue[i,j]=0 ; estimate[i,j]=1}
else {
corr<-correl(x0,y0,method=method,alternative=alternative)
estimate[i,j]<-corr$rho
pvalue[i,j]<-corr$pvalue
}
}
}
names(method)=""
estimate<-round(estimate,2)
diag(pvalue)<-1
cat("\nCorrelation Analysis\n\nMethod     :",method)
cat("\nAlternative:",alternative,"\n\n")
n1<-unique(c(nn))
if(length(n1)==1)nn<-n1
lista<-list(correlation=estimate,pvalue=pvalue, n.obs=nn)
return(lista)
#-------------
}
}

