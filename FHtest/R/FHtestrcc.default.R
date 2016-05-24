FHtestrcc.default <-
function(L, R, group, rho=0, lambda=0,alternative,...){
	if(class(group)[1]=="labelled") group <- as.vector(group)
	if((!is.factor(group))&(!(is.vector(group)&is.numeric(group)))&(!(is.vector(group)&is.character(group)))) stop("group should be a factor, character, or numeric vector")
	call <- match.call(expand.dots =TRUE)
	L.name <- as.character(c(call$L))
	R.name <- as.character(c(call$R))
	group.name <- as.character(c(call$group))
    if ((length(L.name) > 1) | (length(R.name) > 1) | (length(group.name) > 
        1)) {
        L.name <- "L"
        R.name <- "R"
        group.name <- "group"
    }
    if ((sum(nchar(L.name)) + sum(nchar(R.name)) + sum(nchar(group.name))) > 
        50) {
        L.name <- "L"
        R.name <- "R"
        group.name <- "group"
    }
	if(missing(alternative)||all(alternative!=c("different","increasing","decreasing"))) alternative="different"
	if(sum(L==R)+sum((L<R)&(R==max(R)))!=length(L)) stop("This method needs right-censored data")
	times <- L
	status <- 1*(L==R)
	out<-rcc(times,status,group,rho,lambda)
	out$call <- call
 	out$n = table(group)
	names(out$n)<-paste(group.name,"=",names(out$n),sep="")
	out$data.name <- paste("Data:",paste("{", L.name, ",", R.name, "}", " by ", group.name, sep = ""))
	dif <- t(out$obs-out$exp)
	ug <- sort(unique(group))
	k <- length(ug)
	if ((is.numeric(group))&(k>2)){
		out$information <- paste("Trend FH test for right-censored data",sep = "")
		out$var <- c(t(ug)%*%out$var%*%ug)
		out$statistic <- c((dif%*%ug)/sqrt(out$var))
		if (alternative=="increasing"){
			out$pvalue <- pnorm(out$statistic)
			out$alt.phrase <- paste("Alternative hypothesis: increasing survival functions (higher ", group.name ," implies later event times)", sep = "")
			}
		else if (alternative=="decreasing"){
			out$pvalue <- 1-pnorm(out$statistic)
			out$alt.phrase <- paste("Alternative hypothesis: decreasing survival functions (higher ", group.name, " implies earlier event times)", sep = "")
			}
		else {
			out$pvalue <- 2-2*pnorm(abs(out$statistic))
			out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal", sep = "")
			}
		names(out$statistic) <- "Z"
		}
	else if (!(is.numeric(group))&(k>2)){
		out$information <- paste("K-sample test for right-censored data",sep = "")
		out$statistic= c(dif%*%ginv(out$var)%*%t(dif))
		if (alternative!="different") warning("alternative ignored, group is factor with more than 2 groups")
		out$pvalue <- 1 - pchisq(out$statistic,k-1)
		out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal")
		names(out$statistic) <- "Chi Square"
		}
	else {
		out$var <- out$var[1,1]
		out$statistic <- dif[2]/sqrt(out$var)
		out$information <- paste("Two-sample test for right-censored data",sep = "")
		if (alternative=="increasing"){
			out$pvalue <- pnorm(out$statistic)
			out$alt.phrase <- paste("Alternative hypothesis: increasing survival functions (", names(out$n)[2] ," has later event times)", sep = "")
			}
		else if (alternative=="decreasing"){
			out$pvalue <- 1-pnorm(out$statistic)
			out$alt.phrase <- paste("Alternative hypothesis: decreasing survival functions (", names(out$n)[2] ," has earlier event times)", sep = "")
			}
		else{
			out$pvalue <- 2-2*pnorm(abs(out$statistic))
			out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal", sep = "")
			}
		names(out$statistic) <- "Z"
		}

	out$information <- paste("\t",out$information,"\n\nParameters: rho=",as.character(rho), ", lambda=",as.character(lambda),"\nDistribution: counting process approach",sep = "")
	class(out)<-"FHtestrcc"
out
}
