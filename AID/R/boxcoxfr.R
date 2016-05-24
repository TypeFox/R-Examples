
boxcoxfr <- function(y, x, option="both",lam = seq(-2,2,0.02), alpha = 0.05){


x=factor(x)

k=length(levels(x))


if (length(y) != length(x)) {stop("The lengths of x and y must be equal")}

if (is.na(min(y)) == TRUE) {stop("Data include NA")}

if (min(y) <= 0) {stop("Data must include positive values")}

if(((option=="both")|(option=="nor")|(option=="var"))==F) {stop("Correct option argument")}


####################################
if ((option=="both")|(option=="nor")){

stor_w=NULL
for (i in 1:k){

for (j in 1:length(lam)) {

if (lam[j]!=0){
y1=y[which(x==(levels(x)[i]))]
w=(shapiro.test((y1^(lam[j]) - 1)/(lam[j])))
stor_w=rbind(stor_w,c(lam[j],w$statistic,w$p))
}

if (lam[j]==0){
y1=y[which(x==(levels(x)[i]))]
w=shapiro.test(log(y1))
stor_w=rbind(stor_w,c(lam[j],w$statistic,w$p))
}

} #closing for loop

lam=stor_w[which(stor_w[,3]>=alpha),1]
if (length(lam)==0) {stop("Feasible region is null set. No solution.")}
stor_w=NULL

} #closing for loop
}
################################


##########

if ((option=="both")|(option=="var")){
stor_w=NULL
for (j in 1:length(lam)) {

if (lam[j]!=0){
lt=bartlett.test((y^(lam[j]) - 1)/(lam[j]),x)
stor_w=rbind(stor_w,c(lam[j],lt$statistic,lt$p.value))
}

if (lam[j]==0){
lt=bartlett.test(log(y),x)
stor_w=rbind(stor_w,c(lam[j],lt$statistic,lt$p.value))
}
}
lam=stor_w[which(stor_w[,3]>=alpha),1]
if (length(lam)==0) {stop("Feasible region is null set. No solution.")}
}

##########



####

van=boxcox(y~x, lam, plotit = FALSE)
lam=van$x[which.max(van$y)]

####


################################




store=NULL
for(i in 1:k){

if(lam!=0){

kk=shapiro.test((y[which(x==(levels(x)[i]))]^lam-1)/lam)

}else{
kk=shapiro.test(log(y[which(x==(levels(x)[i]))]))
}

store=rbind(store,c(kk$stat,kk$p))

} 

rownames(store)=invisible(paste("Group",levels(x)))
colnames(store)=c("W","p-value")




if(lam!=0){

kk2=bartlett.test((y^lam-1)/lam,x)

}else{

kk2=bartlett.test(log(y),x)
}

store2=matrix(c(kk2$statistic,kk2$parameter,kk2$p.value),1,3)
colnames(store2)=c("K-squared","df","p-value")
rownames(store2)=c("")


out <- list()
out$method <- "MLEFR"
out$date <- date()
out$lambda.hat <-lam
out$shapiro.test <- store
out$bartlett.test <- store2


out

}

