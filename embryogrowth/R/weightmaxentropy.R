#' Search for the weights of the nests which maximize the entropy of nest temperatures distribution
#' @title Search for the weights of the nests which maximize the entropy of nest temperatures distribution
#' @author Marc Girondot
#' @return A named vector of weights
#' @param temperatures Timeseries of temperatures formated using FormatNests()
#' @param entropy.method Entropy function, for example entropy::entropy.empirical. See package entropy for description
#' @param weight A named vector of the initial weight search for each nest for likelihood estimation
#' @param control_optim A list with control paramaters for optim function
#' @param control_plot A list with control paramaters for plot function
#' @param control_entropy A list with control paramaters for entropy function
#' @param plot Do the plot of temperatures before and after weight must be shown ? TRUE or FALSE
#' @param col Colors for unweighted and weighted distributions
#' @description Search for the weights of the nests which maximize the entropy of nest temperatures distribution. Entropy is measured by Shanon index.\cr
#' Entropy method must be entropy.empirical because it is the only method insensitive to scaling.\cr
#' If no weight is given, the initial weight is uniformly distributed.\cr
#' Use control_optim=list(trace=0) for not show progress of search report.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' w <- weightmaxentropy(formated, control_plot=list(xlim=c(20,36)))
#' x <- structure(c(120.940334922916, 467.467455887442,  
#' 	306.176613681557, 117.857995419495),  
#' 	.Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p_weight <- searchR(parameters=x,  
#' 	fixed.parameters=pfixed, temperatures=formated,  
#' 	derivate=dydt.Gompertz, M0=1.7, test=c(Mean=39.33, SD=1.92),  
#' 	method = "BFGS", weight=w)
#' data(resultNest_4p_weight)
#' plotR(resultNest_4p_weight, ylim=c(0,0.50), xlim=c(15, 35))
#' # Standard error of parameters can use the GRTRN_MHmcmc() function
#' }
#' @export


weightmaxentropy <- function(temperatures=stop('Temperature data must be provided !'), 
	weight= NULL, entropy.method=entropy::entropy.empirical, plot=TRUE, 
	control_optim=list(trace=0, maxit=500), control_plot=NULL, 
	control_entropy=NULL, col=c("black", "red")) {
	
# 	weight= NULL; entropy.method=entropy.empirical; plot=TRUE; control_optim=list(trace=0, maxit=500); control_plot=NULL; control_entropy=NULL

  if (!requireNamespace("entropy", quietly = TRUE)) {
    warning("entropy package is necessary for this function")
    return()
  }
  
  
if (is.null(control_plot)) control_plot <- list(NULL)
if (is.null(control_entropy)) control_entropy <- list(NULL)
if (is.null(control_optim)) control_optim <- list(NULL)

NbTS <- temperatures$IndiceT["NbTS"]


if (is.null(weight)) {
	par <- rep(1, NbTS)
	names(par) <- names(temperatures)[1:NbTS]
} else {

	if (any(is.na(weight))) {
		par <- rep(1, NbTS)
		names(par) <- names(temperatures)[1:NbTS]
	} else {

		if (is.list(weight)) weight <- weight$weight

		if (length(setdiff(names(temperatures)[1:NbTS], names(weight)))==0) {
			par <- weight
		} else {
			print("Check the weights")
			return(invisible())
		}
	}
}

# j'ai le weight mais je dois le classer dans le mme ordre que les temperatures
pec <- NULL
for(i in 1:NbTS) {
	pec <- c(pec, par[which(names(par)==names(temperatures)[i])])
}

par <- pec[1:(NbTS-1)]/pec[NbTS]

vm <- setdiff(names(temperatures)[1:NbTS], names(par))

if (length(vm)!=1) {
				print("Check the weights")
				return(invisible())
}




control_optim <- modifyList(control_optim, list(fnscale=-1))


##################
# Calcul de la distribution sans weight
##################

pi <- NULL

NbTS <- temperatures[["IndiceT"]][3]
  
vec.weight <- c(abs(par), 1)

for(series in 1:NbTS) {

	nids <- temperatures[[series]]

# Je cree un tableau avec les donnees heures par heure
	tl1 <- (0:(nids[,1][length(nids[,1])]%/%60-1))*60
# je prends les vraies donnees
	tl2 <- nids[,1]


	it <- findInterval(tl1, tl2)
	tls <- nids[it, 2]
  
  ctls <- cut(tls, 0:45)
	
	pi <- rbind(pi, table(ctls))

}
  
  rownames(pi) <- names(temperatures)[1:NbTS]
  
  nrow <- dim(pi)[1]
  ncol <- dim(pi)[2]
   
  pisansweight <- NULL
  for(i in 1:ncol) pisansweight <- c(pisansweight, sum(pi[,i]))


##################
# Recherche des weights pour entropy maximale
##################



  print(paste("Initial entropy is", .weightmaxentropy_fit(par, pi=pi, entropy.method=entropy.method, pentropy=control_entropy, vm=vm)))


  result  <- optim(par, .weightmaxentropy_fit, pi=pi, entropy.method=entropy.method, pentropy=control_entropy, vm=vm, control=control_optim, method="BFGS", hessian=FALSE)
  
  if (result$convergence==1) {
  	print("Maximum number of iterations has been reached. Try increase maxit in controls")
  }
  
  par <- c(result$par, 1)
  
##################
# Calcul de la distribution avec weight
##################

  
vec.weight <- par
names(par) <- names(temperatures)[1:NbTS]


# juste pour reprendre la structure
pi2 <- pi
  
for(i in 1:nrow)
	pi2[i,] <- pi[i, 1:ncol]*(vec.weight[i]/sum(vec.weight))
  
  piavecweight <- NULL
  for(i in 1:ncol) piavecweight <- c(piavecweight, sum(pi2[,i]))

if (plot) {

  rm(plot)
L <- modifyList(list(x=0:44, y=pisansweight, type="h", 
	xlim=c(0,44), lwd=2, bty="n", xlab="Temperatures", 
	ylab="Relative frequency", col=col[1]), control_plot)
	
	do.call(plot, L)

# if (par("yaxs")=="i") {
#    y1 <- par("usr")[3]
#    y2 <- par("usr")[4]
#  } else {
#    y2 <- (par("usr")[3]+par("usr")[4]*26)/27
#    y1 <- y2*26-par("usr")[4]/0.04
#  }

par(new=TRUE)
L <- modifyList(list(x=(0:44)+0.3, y=piavecweight, type="h", 
	xlim=c(0,44), lwd=2, bty="n", xlab="", axes=FALSE,
	ylab="", col=col[2]), control_plot)

	do.call(plot, L)

legend("topleft", legend=c("Without weight", "With weight"), col=col, lty=1, lwd=2, bty="n")

}
  
  
  print(paste("Maximal entropy is", result$value))
  
  return(list(weight=abs(par), freq.obs=pisansweight, freq.max.entropy=piavecweight))

  
}
