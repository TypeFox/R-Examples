FHtestrcp.default <-
function(L, R, group, rho=0, lambda=0,alternative,method=NULL,methodRule = methodRuleIC1,exact = NULL,permcontrol=permControl(),...){
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

	cc <- rcp(times,status,rho,lambda)
	names(cc)<-1:length(times)

	if (is.null(method)) method <- methodRule(cc, group, exact)
	method.OK <- (method == "pclt" | method == "exact.mc" | method == "exact.network" | method == "exact.ce")
	if (!method.OK) stop("method not one of: 'pclt', 'exact.mc'. 'exact.network', 'exact.ce'")

	ug <- sort(unique(group))
	k <- length(ug)
	if ((is.numeric(group))&(k>2)) {
		icalternative <- c("two.sided","less","greater")[alternative==c("different","increasing","decreasing")]
		icp <- do.call("permTREND", list(x = cc, y = group, alternative=icalternative, method = method,control=permcontrol))
		}
	else if (!(is.numeric(group))&(k>2)) {
		icp <- do.call("permKS", list(x = cc, g = group, alternative="two.sided", method = method,control=permcontrol))
		}
	else {
		icalternative <- c("two.sided","greater","less")[alternative==c("different","increasing","decreasing")]
		icp <- do.call("permTS", list(x = cc[group==ug[1]], y = cc[group==ug[2]], alternative=icalternative, method = method,control=permcontrol))
		}

	out <- list(cc,icp$statistic,icp$p.value,icp$p.conf.int)
	names(out)<-c("scores","statistic","pvalue","p.conf.int")
	out$diff <- aggregate(out$scores,list(group),sum)$x
	out$call <- call
 	out$n = table(group)
	names(out$n)<-paste(group.name,"=",names(out$n),sep="")
	out$data.name <- paste("Data:",paste("{", L.name, ",", R.name, "}", " by ", group.name, sep = ""))

	n <- length(group)
	z <- 1*(matrix(rep(group,k)==rep(ug,each=n),n))
	out$var <- var(out$scores)*((t(z)%*%z) -n*(t(t(colMeans(z)))%*%t(colMeans(z))))

	if ((is.numeric(group))&(k>2)){
		out$var <- c(t(ug)%*%out$var%*%ug)
		out$information <- paste("Trend FH test for right-censored data",sep = "")
		if (method!="pclt") out$statistic <- (do.call("permTREND", list(x = out$scores, y = group, alternative="two.sided", method = "pclt")))$statistic
		if (alternative=="increasing"){
			out$alt.phrase <- paste("Alternative hypothesis: increasing survival functions (higher ", group.name ," implies later event times)", sep = "")
			}
		else if (alternative=="decreasing"){
			out$alt.phrase <- paste("Alternative hypothesis: decreasing survival functions (higher ", group.name, " implies earlier event times)", sep = "")
			}
		else {
			out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal", sep = "")
			}
		}
	else if (!(is.numeric(group))&(k>2)){
		out$information <- paste("K-sample test for right-censored data",sep = "")
		if (method!="pclt") out$statistic <- (do.call("permKS", list(x = out$scores, g = group, method = "pclt")))$statistic
		if (alternative!="different") warning("alternative ignored, group is factor with more than 2 groups")
		out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal")
		}
	else {
		out$var <- out$var[1,1]
		out$information <- paste("Two-sample test for right-censored data",sep = "")
		if (method=="pclt") out$statistic <- -out$statistic
		else {
			ug <- sort(unique(group))
			out$statistic <- -(do.call("permTS", list(x = out$scores[group==ug[1]], y = out$scores[group==ug[2]], alternative="two.sided", method = "pclt")))$statistic
			}
		if (alternative=="increasing"){
			out$alt.phrase <- paste("Alternative hypothesis: increasing survival functions (", names(out$n)[2] ," has later event times)", sep = "")
			}
		else if (alternative=="decreasing"){
			out$alt.phrase <- paste("Alternative hypothesis: decreasing survival functions (", names(out$n)[2] ," has earlier event times)", sep = "")
			}
		else{
			out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal", sep = "")
			}
		}

	out$information <- paste("\t",out$information,"\n\nParameters: rho=",as.character(rho), ", lambda=",as.character(lambda),"\nDistribution: permutation approach (",sep = "")
	if (method=="pclt") out$information <- paste(out$information,"asymptotic approximation using central limit theorem)",sep = "")
	else if (method=="exact.ce") out$information <- paste(out$information,"exact method using complete enumeration)",sep = "")
	else if (method=="exact.mc") out$information <- paste(out$information,"exact method using Monte Carlo with ",as.character(icp$nmc)," replications)",sep = "")
	else if (method=="exact.network") out$information <- paste(out$information,"exact method using a network algorithm)",sep = "")

	class(out)<-"FHtestrcp"
out
}
