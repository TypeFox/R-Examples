FHtestics.default <-
function(L, R, group, rho=0, lambda=0, alternative, tol=10^-8,icFIT=NULL,initfit=NULL,icontrol=icfitControl(),Lin=NULL,Rin=NULL,...){
	if(class(group)[1]=="labelled") group <- as.vector(group)
	if(lambda>0) stop("The scoretest needs lambda to be 0")
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
	if(missing(alternative)||all(alternative!=c("different","increasing","decreasing"))) alternative <- "different"

	k <- length(unique(group))
	ord <- order(group)
	if(length(Lin)>1) Lin <- Lin[ord]
	if(length(Rin)>1) Rin <- Rin[ord] 
	if(!is.null(icFIT)) icFIT$A <- icFIT$A[ord,]
	if(!is.null(initfit)) initfit$A <- initfit$A[ord,]
	limrho<-rho
	if(rho<10^-7) limrho<-10^-7
	icp <-ictest(L[ord],R[ord],group[ord], scores = "general",dqfunc=function(x){x*pbeta(1-x,1,limrho)*beta(1,limrho)}, method="pclt",icFIT=icFIT,initfit=initfit,icontrol=icontrol,Lin=Lin,Rin=Rin,...)		
	icp$scores[ord] <- icp$scores

	if ((is.numeric(group))&(k>2)){
		icsout <- ics(icp$fit,group[ord],rho,"less",tol)
		icp$fit$A[ord,] <- icp$fit$A
		out <- list(icp$fit,icp$scores,icsout$statistic,as.numeric(icsout$p.values[c(alternative==c("different","increasing","decreasing"),FALSE)]),c(icsout$V),icsout$d2L.dB2,icsout$d2L.dgam2,icsout$d2L.dBdgam)
		names(out[[4]]) <- NULL
		}
	else if (!(is.numeric(group))&(k>2)){
		icsout <- ics(icp$fit,group[ord],rho,"two.sided",tol)
		icp$fit$A[ord,] <- icp$fit$A
		out <- list(icp$fit,icp$scores,c(icsout$statistic),c(icsout$p.value),icsout$V,icsout$d2L.dB2,icsout$d2L.dgam2,icsout$d2L.dBdgam)
		}

	else	{
		icsout <- ics(icp$fit,as.character(group[ord]),rho,"less",tol)
		icp$fit$A[ord,] <- icp$fit$A
		out <- list(icp$fit,icp$scores,-icsout$statistic,as.numeric(icsout$p.values[c(alternative==c("different","decreasing","increasing"),FALSE)]),icsout$V,icsout$d2L.dB2,icsout$d2L.dgam2,icsout$d2L.dBdgam)
		}

	
	names(out)<-c("fit","scores","statistic","pvalue","var","d2L.dB2","d2L.dgam2","d2L.dBdgam")
	out$diff <- aggregate(out$scores,list(group),sum)$x
	out$call <- call
 	out$n = table(group)
	names(out$n)<-paste(group.name,"=",names(out$n),sep="")
	out$data.name <- paste("Data:",paste("{", L.name, ",", R.name, "}", " by ", group.name, sep = ""))

	if ((is.numeric(group))&(k>2)){
		out$information <- paste("Trend FH test for interval-censored data",sep = "")
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
		out$information <- paste("K-sample test for interval-censored data",sep = "")
		if (alternative!="different") warning("alternative ignored, group is factor with more than 2 groups")
		out$alt.phrase <- paste("Alternative hypothesis: survival functions not equal")
		}
	else {
		out$var <- out$var[1,1]
		out$information <- paste("Two-sample test for interval-censored data",sep = "")
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

	out$information <- paste("\t",out$information,"\n\nParameters: rho=",as.character(rho), ", lambda=0","\nDistribution: score vector approach",sep = "")
	class(out)<-"FHtestics"
out
}
