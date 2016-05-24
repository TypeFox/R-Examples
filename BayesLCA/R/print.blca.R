print.blca <-
function(x, ...){
z<-NULL

btau<- x$classprob
btheta<- x$itemprob
btau.se<- x$classprob.se
btheta.se<- x$itemprob.se

names(btau)<- paste(rep("Group", length(btau)), 1:length(btau))

if(is.matrix(btheta)) {
  if(is.null(colnames(btheta))) colnames(btheta)<- paste(rep("Col", length(btau)), 1:ncol(btheta))
} else   if(is.null(names(btheta))) names(btheta)<- paste(rep("Col", length(btau)), 1:ncol(btheta))

rownames(btheta)<- paste(rep("Group", length(btau)),1:nrow(btheta))

if(!is.null(btau.se)){
names(btau.se)<- names(btau)
dimnames(btheta.se)<- dimnames(btheta)
}

z$ItemProb<- btheta
z$MembershipProb<- btau
if(!is.null(btheta.se)) z$se$ItemProb<- btheta.se
if(!is.null(btau.se)) z$se$MembershipProb<- btau.se
#class(z)<-"summary.blca"
#return(z)

	cat("\nMAP Estimates:\n", "\n")
	cat("\nItem Probabilities:\n", "\n")
	print(round(z$ItemProb,3))
	
	cat("\nMembership Probabilities:\n", "\n")
	print(round(z$MembershipProb,3))

if(!is.null(btheta.se)){	
	cat("\nPosterior Standard Deviation Estimates:\n", "\n")
	cat("\nItem Probabilities:\n", "\n")
	print(round(z$se$ItemProb, 3))
	
	cat("\nMembership Probabilities:\n", "\n")
	print(round(z$se$MembershipProb, 3))
} else warning("Posterior standard deviations not returned.", call.=FALSE)

}
