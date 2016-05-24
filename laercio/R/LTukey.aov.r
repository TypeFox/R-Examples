`LTukey.aov` <- function(anova,which="",conf.level=0.95)
{
	nomes <- names(anova$model)
	for (i in 1:length(nomes)){
		if (nomes[i]==which){
			variavel=TRUE
			vari = i;
			break	}
			else	{variavel=FALSE}	}
cat("\n","TUKEY TEST TO COMPARE MEANS","\n","\n",
"Confidence level: ",conf.level,"\n",
"Dependent variable: ",nomes[1])
CV <- sqrt(sum(anova$residuals^2)/anova$df.residual)/mean(anova$model[[1]])
cat("\n","Variation Coefficient: ",CV*100,"%","\n","\n")
if (variavel==TRUE){
	r <- tapply(anova$model[[1]],anova$model[[vari]],length)
	sum.sq.res <- tapply(anova$residuals^2,anova$model[[vari]],sum)
	sy <- sum(sum.sq.res/r)/anova$df.residual
	medias <- tapply(anova$model[[1]],anova$model[[vari]],mean)
	nmeans=length(medias)
	n <- nmeans
	q <- qtukey(conf.level,nmeans,anova$df.residual)
	DMS <- q*sqrt(sy)
	cat("Independent variable: ", nomes[vari],"\n")
	l <- rep("",length(medias))
	nd <- 0
	medias <- sort(medias,decreasing=T)
	start <- 1
	stop <- 0
	repeat {
		dif <- medias[start]-medias[n]
		while (dif > DMS)	{
			n <- n - 1
			dif <- medias[start]-medias[n]}
		if (start<n & stop<n){
			l[start:n] <- paste(l[start:n],letters[nd+1],sep="")
			if (n < length(l)) {
				l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="")}
			nd <- nd + 1 }
		if (start == n && start != stop) {
			l[n] <- paste(l[n],letters[nd+1],sep="")
			nd <- nd + 1
			l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="") }
		if (medias[length(medias)-1] - medias[length(medias)] > DMS && start==length(medias)-1) {
			l[length(medias)] <- paste(l[length(medias)],letters[nd+1],sep="") }
	stop <- n
	n <- nmeans
	start <- start + 1;
	if (start==n) {break} }
resultado <- cbind("Factors"=names(medias),"Means"=as.numeric(medias)," "=l)
resultado <- as.table(resultado)
row.names(resultado) <- rep(" ",nmeans)
print(resultado)
cat("\n","\n")}
else {
	for (i in 2:length(nomes)) {	
	r <- tapply(anova$model[[1]],anova$model[[i]],length)
	sum.sq.res <- tapply(anova$residuals^2,anova$model[[i]],sum)
	sy <- sum(sum.sq.res/r)/anova$df.residual
	medias <- tapply(anova$model[[1]],anova$model[[i]],mean)
	nmeans=length(medias)
	n <- nmeans
	q <- qtukey(conf.level,nmeans,anova$df.residual)
	DMS <- q*sqrt(sy)
	cat("Independent variable: ", nomes[i],"\n")
	l <- rep("",length(medias))
	nd <- 0
	medias <- sort(medias,decreasing=T)
	start <- 1
	stop <- 0
	repeat {
		dif <- medias[start]-medias[n]
		while (dif > DMS) {
			n <- n - 1
			dif <- medias[start]-medias[n] }
		if (start<n & stop<n) {
			l[start:n] <- paste(l[start:n],letters[nd+1],sep="")
			if (n < length(l)) {
				l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="") }
			nd <- nd + 1 }
		if (start == n && start != stop) {
			l[n] <- paste(l[n],letters[nd+1],sep="")
			nd <- nd + 1
			l[(n+1):length(l)] <- paste(l[(n+1):length(l)]," ",sep="") }
		if (medias[length(medias)-1] - medias[length(medias)] > DMS && start==length(medias)-1) {
			l[length(medias)] <- paste(l[length(medias)],letters[nd+1],sep="") }
	stop <- n
	n <- nmeans
	start <- start + 1;
	if (start==n) {break} }
resultado <- cbind("Factors"=names(medias),"Means"=as.numeric(medias)," "=l)
resultado <- as.table(resultado)
row.names(resultado) <- rep(" ",nmeans)
print(resultado)
cat("\n","\n") } } }