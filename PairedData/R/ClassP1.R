##################################################" nouvelle classe


validPairedObject<-function(object){
test1<-(dim(object)[2]==2)
if(test1) {test<-TRUE} else 
{
return("paired data must have two columns")
}
test2<-class(object[,1])==class(object[,2])
if(test2) {test<-TRUE} else 
{
return("paired data must have the same classes")
}
if(is.factor(object[,1])) {test3<-all(levels(object[,1])==levels(object[,2]))
if(test3) {test<-TRUE} else 
{
return("paired factor data must have the same levels")
}
}
return(test)
}


### Les constructeurs

setClass(Class="paired",contains="data.frame",validity=validPairedObject)

paired<-function(x,y){
object<-new(Class="paired",data.frame(x,y))
name.x<-deparse(substitute(x))
name.y<-deparse(substitute(y))
colnames(object)<-c(name.x,name.y)
object
}



setMethod("summary",
          signature(object = "paired"),
          function(object,tr=0.2){
if(is.numeric(object[,1])){
X<-paired.summary(object[,1],object[,2],tr=tr)
rownames(X)<-c(paste(colnames(object)[1]," (x)",sep=""),paste(colnames(object)[2]," (y)",sep=""),"x-y","(x+y)/2")
Y<-matrix(numeric(2),ncol=2,nrow=1)
rownames(Y)<-"(x,y)"
colnames(Y)<-c("cor","wcor")
Y[1,1]<-cor(object[,1],object[,2])
Y[1,2]<-wincor(object[,1],object[,2],tr=tr)$cor
list(stat=X,cor=Y)
} 
else{
summary.data.frame(object)}
          }
)



setMethod(f="effect.size", 
signature(object = "paired"),
  definition=function(object,tr=0.2){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}
paired.effect.size(object[,1],object[,2],tr=tr)
}
)

setMethod("plot", signature="paired",
  function(x,groups=NULL,subjects=NULL,facet=TRUE,type=c("correlation","BA","McNeil","profile"),...){
# plot type
type <- match.arg(type)

# paired object
df<-x
if(is.null(colnames(df))){colnames(df)<-c("C1","C2")}
conditions<-colnames(df)

# subjects
if(is.null(subjects)){
n<-dim(df)[1]
subjects<-factor(paste("S",1:n,sep=""))
df<-data.frame(df,subjects)
}
else{
df<-data.frame(df,subjects)
name.subjects<-deparse(substitute(subjects))
colnames(df)[3]<-name.subjects 
}

if(!is.null(groups)){
# with groups
df<-data.frame(df,groups)
colnames(df)[4]<-deparse(substitute(groups))
if(type=="correlation"){
return(paired.plotCor(df,conditions[1],conditions[2],groups=colnames(df)[4],facet=facet,...))}
if(type=="BA"){
return(paired.plotBA(df, conditions[1], conditions[2], groups = colnames(df)[4], facet = facet,...))
}
if(type=="McNeil"){
return(paired.plotMcNeil(df, conditions[1], conditions[2], groups = colnames(df)[4], subjects= colnames(df)[3],facet = facet,...))
}
if(type=="profile"){
return(paired.plotProfiles(df, conditions[1], conditions[2], groups = colnames(df)[4], subjects= colnames(df)[3],facet = facet,...))
}

}
else{
# without groups
if(type=="correlation"){
return(paired.plotCor(df,conditions[1],conditions[2],groups=NULL,
facet=facet,...))}
if(type=="BA"){
return(paired.plotBA(df, conditions[1], conditions[2], groups = NULL, facet = facet,...))
}
if(type=="McNeil"){
return(paired.plotMcNeil(df, conditions[1], conditions[2], groups = NULL, subjects= colnames(df)[3],facet = facet,...))
}
if(type=="profile"){
return(paired.plotProfiles(df, conditions[1], conditions[2], groups = NULL, subjects= colnames(df)[3],facet = facet,...))
}



}
}
)



setMethod(f="slidingchart", 
signature="paired",
  definition=function(object,...){
if(!is.numeric(object[,1])){return("only suitable to numeric paired data")}

df<-object
condition1<-colnames(df)[1]
condition2<-colnames(df)[2]
return(paired.slidingchart(df, condition1, condition2,...))

}
)


#### Enfin les fonctions de simulation :

rpaired.contaminated<-
function (n, d1 = c(0.1, 10, 1), d2 = c(0.1, 10, 1), r = 0.5) 
{
    require(mvtnorm)
    eps1 <- d1[1]
    k1 <- d1[2]
    Sigma1 <- d1[3]
    eps2 <- d2[1]
    k2 <- d2[2]
    Sigma2 <- d2[3]
    X <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, r, r, 
        1), ncol = 2))
    u1 <- pnorm(X[, 1])
    b1 <- rbinom(n, size = 1, prob = eps1)
    SD1 <- Sigma1 * (b1 * k1 + (1 - b1) * 1)
    x <- qnorm(u1, mean = 0, sd = SD1)
    u2 <- pnorm(X[, 2])
    b2 <- rbinom(n, size = 1, prob = eps2)
    SD2 <- Sigma2 * (b2 * k2 + (1 - b2) * 1)
    y <- qnorm(u2, mean = 0, sd = SD2)
    return(paired(x,y))
}

rpaired.gld<-
function (n, d1=c(0.000,0.1974,0.1349,0.1349), d2=c(0.000,0.1974,0.1349,0.1349), r) 
{
    require(mvtnorm)
    require(gld)
    X <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, r, r, 
        1), ncol = 2))
    u1 <- pnorm(X[, 1])
    u2 <- pnorm(X[, 2])
    x <- qgl(u1, lambda1 = d1[1], lambda2 = d1[2], lambda3 = d1[3], 
        lambda4 = d1[4], param = "rs")
    y <- qgl(u2, lambda1 = d2[1], lambda2 = d2[2], lambda3 = d2[3], 
        lambda4 = d2[4], param = "rs")
    return(paired(x,y))
}

