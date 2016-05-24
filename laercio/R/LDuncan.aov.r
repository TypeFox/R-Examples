`LDuncan.aov` <- function(anova,which="",conf.level=0.95){
prob <- function(alfa=conf.level,nmedias,dfresidual) {
alfa <- 1-alfa
a <- (1-alfa)^(nmedias-1)
res <- qtukey(a,nmedias,dfresidual)
return(res)}
nomes <- names(anova$model)
for (i in 1:length(nomes)) {
if (nomes[i]==which) {
variavel=TRUE
vari = i
break}
else	{variavel=FALSE} }
cat("\n","DUNCAN TEST TO COMPARE MEANS","\n","\n",
"Confidence Level: ",conf.level,"\n",
"Dependent Variable: ",nomes[1])
CV <- sqrt(sum(anova$residuals^2)/anova$df.residual)/mean(anova$model[[1]])
cat("\n","Variation Coefficient: ",CV*100,"%","\n","\n")
if (variavel==TRUE){
r <- tapply(anova$model[[1]],anova$model[[vari]],length)
sum.sq.res <- tapply(anova$residuals^2,anova$model[[vari]],sum)
sy <- sum(sum.sq.res/r)/anova$df.residual
medias <- tapply(anova$model[[1]],anova$model[[vari]],mean)
nmeans=length(medias)
n <- nmeans
nn <- nmeans
cat("\n","Independent Variable: ", nomes[vari],"\n")
l <- rep("",length(medias))
nd <- 0
medias <- sort(medias,decreasing=T) #coloca as medias em ordem crescente
start <- 1
stop <- 0 #o ultimo numero que teve diferenca sifnificativa
Di <- 0
repeat{
dif <- medias[start]-medias[n]
nn <- n-start+1
while (dif > Di){
if (nn==1) {stop}
Zi <- prob(conf.level,nn,anova$df.residual)
Di <- Zi*sqrt(sy)
n <- n - 1
nn <- n-start+1			
dif <- medias[start]-medias[n]}
if (start<n & stop<n){
l[start:n] <- paste(l[start:n],letters[nd+1],sep="")
if (n < length(l)){
l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="")}
nd <- nd + 1}
if (start == n && start != stop){
l[n] <- paste(l[n],letters[nd+1],sep="")
nd <- nd + 1
l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="")}
if (start==length(medias)-1) {
Zi <- prob(conf.level,2,anova$df.residual)
Di <- Zi*sqrt(sy)
if (medias[length(medias)-1] - medias[length(medias)] > Di){
l[length(medias)] <- paste(l[length(medias)],letters[nd+1],sep="")}}
stop <- n
n <- nmeans
start <- start + 1;
if (start==n) {break}}
resultado <- cbind("Factors"=names(medias),"Means"=as.numeric(medias)," "=l)
resultado <- as.table(resultado)
row.names(resultado) <- rep(" ",nmeans)
print(resultado)}
else{
for (i in 2:length(nomes)){	
r <- tapply(anova$model[[1]],anova$model[[i]],length)
sum.sq.res <- tapply(anova$residuals^2,anova$model[[i]],sum)
sy <- sum(sum.sq.res/r)/anova$df.residual
medias <- tapply(anova$model[[1]],anova$model[[i]],mean)
nmeans=length(medias)
n <- nmeans
q <- qtukey(conf.level,nmeans,anova$df.residual)
DMS <- q*sqrt(sy)
cat("\n","Independent Variable: ", nomes[i],"\n")
l <- rep("",length(medias))
nd <- 0
medias <- sort(medias,decreasing=T)
start <- 1
stop <- 0
Di <- 0
repeat {
dif <- medias[start]-medias[n]
nn <- n-start+1
while (dif > Di){
if (nn==1) {stop}
Zi <- prob(conf.level,nn,anova$df.residual)
Di <- Zi*sqrt(sy)
n <- n - 1
nn <- n-start+1			
dif <- medias[start]-medias[n]}
if (start<n & stop<n){
l[start:n] <- paste(l[start:n],letters[nd+1],sep="")
if (n < length(l)){
l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="")}
nd <- nd + 1}
if (start == n && start != stop){
l[n] <- paste(l[n],letters[nd+1],sep="")
nd <- nd + 1
l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="")}
if (start==length(medias)-1) {
Zi <- prob(conf.level,2,anova$df.residual)
Di <- Zi*sqrt(sy)
if (medias[length(medias)-1] - medias[length(medias)] > Di){
l[length(medias)] <- paste(l[length(medias)],letters[nd+1],sep="")}}
stop <- n
n <- nmeans
start <- start + 1;
if (start==n) {break}}
resultado <- cbind("Factors"=names(medias),"Means"=as.numeric(medias)," "=l)
resultado <- as.table(resultado)
row.names(resultado) <- rep(" ",nmeans)
print(resultado)
cat("\n")}}
}