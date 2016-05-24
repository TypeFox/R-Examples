#' plotR shows the fitted growth rate dependent on temperature
#' @title Show the fitted growth rate dependent on temperature
#' @author Marc Girondot
#' @return Nothing
#' @param result A result object or a list of result objects
#' @param ... Parameters for plot() such as main= or ylim=
#' @param parameters Indicate some parameters if the result object is not supplied
#' @param fixed.parameters Indicate some parameters if the result object is not supplied
#' @param SE The standard error for the parameters or a list of SE if several results. Use NA to force not use SE
#' @param set.par 1 or 2 or a list of 1 or 2 to designate with set of parameters to show
#' @param size If indicated, will show the growth rate for this size
#' @param legend Text to show in bottom right legend or a list of text if several results
#' @param col The color to use for a list of colors if several results
#' @param lty The type of line to use if several results as a list
#' @param ltyCI The type of line to use for confidence interval as a list
#' @param lwd The type of line to use if several results as a list
#' @param lwdCI The type of line to use for confidence interval as a list
#' @param xlim Range of values for x-axis
#' @param xlimR Range of values to be displayed for R curve; can be a list if a list of results is used
#' @param scaleY Scaling factor for y axis or "auto"
#' @param replicate.CI Number of randomizations to estimate CI
#' @param show.box If TRUE show a box with "mean" and "confidence interval"
#' @param local.box Position of the box with "mean" and "confidence interval", default="topleft"
#' @description To show the growth rate, the syntaxe is:\cr
#' plotR(result=res)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_4p_80 <- GRTRN_MHmcmc(result=resultNest_4p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' plotR(result=resultNest_4p, SE=result_mcmc_4p$TimeSeriesSE,  
#' ylim=c(0,0.3))
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' plotR(parameters=x, xlim=c(20,35))
#' pfixed <- c(rK=2.093313)
#' resultNest_6p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_6p)
#' plotR(list(resultNest_4p, resultNest_6p),  
#' col=c("black", "red"), legend=c("4 parameters", "6 parameters"))
#' ##########################################
#' # new formulation of parameters
#' data(resultNest_newp)
#' plotR(resultNest_newp)
#' }
#' @export


plotR <-
function(result=NULL, parameters=NULL, fixed.parameters=NULL, col="black", legend=NA, 
	SE=NULL, set.par=1, size=NA, xlim=c(20,35), scaleY="auto", lty=1, ltyCI=3, lwd=1, lwdCI=1, 
  xlimR=xlim, replicate.CI=100, show.box=TRUE, local.box="topleft", ...) {

  # result=NULL;parameters=NULL; fixed.parameters=NULL; lty=1; ltyCI=3; lwd=1; lwdCI=1; col="black"; legend=NA; SE=NULL; set.par=1; size=NA; xlim=c(20,35); xlimR=xlim; scaleY="auto"; replicate.CI=100; show.box=TRUE; local.box="topleft"
  # result <- resultNest_4p; parameters <- newp; SE <- result_mcmc_newp$TimeSeriesSE; ylim <- c(0,0.4)
  
if (is.null(result) & is.null(c(parameters, fixed.parameters))) {
	print("Or parameters or result from searchR or likelihoodR must be provided !")
	return()
}

afficheCI <- FALSE

if (!is.list(col)) col <- list(col)
if (!is.list(lty)) lty <- list(lty)
if (!is.list(ltyCI)) ltyCI <- list(ltyCI)
if (!is.list(lwd)) lwd <- list(lwd)
if (!is.list(lwdCI)) lwdCI <- list(lwdCI)
if (!is.list(set.par)) set.par <- list(set.par)
if (!is.list(SE)) SE <- list(SE)
if (!is.list(xlimR)) xlimR <- list(xlimR)
if (!is.list(legend)) legend <- list(legend)
if (!is.list(parameters)) parameters <- list(parameters)


if (!is.null(result)) {
		if (class(result)!="list") result <- list(result)
		
			nbr <- max(length(result), length(set.par), length(SE), length(legend), length(col), length(parameters))
# sinon il faut que je complte recycle col et mette legend et SE  NA
			result <- c(result, rep(result, nbr-length(result)))
			col <- as.list(rep(unlist(col), nbr)[1:nbr])
			lty <- as.list(rep(unlist(lty), nbr)[1:nbr])
      ltyCI <- as.list(rep(unlist(ltyCI), nbr)[1:nbr])
      lwd <- as.list(rep(unlist(lwd), nbr)[1:nbr])
      lwdCI <- as.list(rep(unlist(lwdCI), nbr)[1:nbr])
			set.par <- as.list(rep(unlist(set.par), nbr)[1:nbr])
			SE <- c(SE, rep(list(NULL), nbr-length(SE)))
			for(rs in 1:length(result)) {
			  if (is.null(SE[[rs]])) SE[[rs]] <- result[[rs]]$SE
			}
			
			xlimR <- c(xlimR, rep(xlimR, nbr-length(xlimR)))
			legend <- as.list(c(unlist(legend), rep(NA, nbr-length(legend))))
			parameters <- c(parameters, rep(list(NULL), nbr-length(parameters)))

	} else {
# j'ai des paramtres
	  
	  
	  
	  nbr <- max(length(set.par), length(SE), length(legend), length(col), length(parameters))
	  
			result <- list(NA)
			result <- as.list(rep(unlist(result), nbr)[1:nbr])
#			col <- as.list(col)
#			lty <- as.list(lty)
#			SE <- as.list(SE)
#			legend <- as.list(legend)
#			parameters <- as.list(parameters)
 # 		set.par <- as.list(set.par)
	  col <- as.list(rep(unlist(col), nbr)[1:nbr])
	  lty <- as.list(rep(unlist(lty), nbr)[1:nbr])
    ltyCI <- as.list(rep(unlist(ltyCI), nbr)[1:nbr])
    lwd <- as.list(rep(unlist(lwd), nbr)[1:nbr])
    lwdCI <- as.list(rep(unlist(lwdCI), nbr)[1:nbr])
	  set.par <- as.list(rep(unlist(set.par), nbr)[1:nbr])
	  SE <- c(SE, rep(list(NA), nbr-length(SE)))
	  xlimR <- c(xlimR, rep(xlimR, nbr-length(xlimR)))
	  legend <- as.list(c(unlist(legend), rep(NA, nbr-length(legend))))
	  parameters <- c(parameters, rep(list(NULL), nbr-length(parameters)))
}

premier <- TRUE

for (rs in 1:nbr) {

# J'introduis les paramtres fixes - 16/7/2012
if (is.na(result[rs])) {
	parssm <- c(parameters[[rs]], fixed.parameters)
	res <- SE[[rs]]
} else {
	if (is.null(parameters[[rs]])) {
		parssm <- c(result[[rs]]$par, result[[rs]]$fixed.parameters)
		if (all(!is.na(SE[[rs]]))) {res <- c(SE[[rs]], result[[rs]]$SE)} else {res <- NA}
	} else {
		parssm <- c(parameters[[rs]], result[[rs]]$fixed.parameters)
		if (all(!is.na(SE[[rs]]))) {res <- c(SE[[rs]], result[[rs]]$SE)} else {res <- NA}
	}
}

# je suis en Anchor
if (all(names(parssm)!="Rho25")) {
  xlR <- c(min(as.numeric(names(parssm[(names(parssm)!="rK") & (names(parssm)!="K") & (names(parssm)!="Scale")])), na.rm=TRUE), max(as.numeric(names(parssm[(names(parssm)!="rK") & (names(parssm)!="K") & (names(parssm)!="Scale")])), na.rm=TRUE))
  if (xlR[1]>273) xlR <- xlR-273.15
  xlR <- c(max(xlR[1], xlimR[[rs]][1]), min(xlR[2], xlimR[[rs]][2]))

  } else {
  xlR <- xlimR[[rs]]
}
  
x <- seq(xlR[1],xlR[2],by=0.1)
# if (x<273) x <- x+273.15
voutlist <- .SSM(x, parssm)
# voutlist <- embryogrowth:::.SSM(x, parssm)



if (!is.na(size) & !is.na(parssm["transition_S"]) & !is.na(parssm["transition_P"])) {
r <- voutlist[[1]]
r_L <- voutlist[[2]]
transition <- 1/(1+exp(parssm["transition_S"]*(size-parssm["transition_P"])))
vout <- r*transition+r_L*(1-transition)

} else {
vout <- voutlist[[set.par[[rs]]]]
}

if (scaleY=="auto") scaleY <- 10^(-floor(log10(max(vout)))-1)

y<- scaleY*vout

if (premier) {
	L <- modifyList(list(type = "l", las=1, col=col[[rs]], lty=lty[[rs]], lwd=lwd[[rs]], axes = TRUE, bty = "n", xlab = expression("Temperatures in " * degree * "C"), ylab = paste("r*", scaleY, sep=""), xlim=xlim), modifyList(list(x=x, y=y), list(...)))
} else {
	L <- modifyList(modifyList(list(x=x, y=y), list(...)), list(type = "l", las=1, col=col[[rs]], lty=lty[[rs]], lwd=lwd[[rs]], axes = FALSE, bty = "n", xlab = "", ylab = "", xlim=xlim, ylim=c(y1, y2), main="")) 
}

# tp <- NULL
# if (premier) {
#  L <- modifyList(list(type = "l", las=1, col=col[[rs]], lty=lty[[rs]], lwd=lwd[[rs]], axes = TRUE, bty = "n", xlab = expression("Temperatures in " * degree * "C"), ylab = paste("r*", scaleY, sep=""), xlim=xlim), list(x=x, y=y, tp)) 
# } else {
#  L <- modifyList(list(x=x, y=y, tp), list(type = "l", las=1, col=col[[rs]], lty=lty[[rs]], lwd=lwd[[rs]], axes = FALSE, bty = "n", xlab = "", ylab = "", xlim=xlim, ylim=c(y1, y2), main="")) 
# }


if (length(which(names(L)=="show.box"))!=0) {
  L <- L[-which(names(L)=="show.box")]
}

if (length(which(names(L)=="xlimR"))!=0) {
  L <- L[-which(names(L)=="xlimR")]
}


do.call(plot, L) 

y2 <- (par("usr")[3]+par("usr")[4]*26)/27
y1 <- y2*26-par("usr")[4]/0.04
# ylim=c(y1, y2)

if (!is.null(res)) {
if (!all(is.na(res))) {

## Nouvelle methode prenant beaucoup moins de memoire
  

# ess <- list(Parametre=matrix(rep(NA, 16*replicate.CI), ncol=16, dimnames=list(NULL, c("DHA", "DHL", "DHH", "T12L", "T12H", "DT", "Rho25", "DHA_L", "DHL_L", "DHH_L", "T12L_L", "T12H_L", "DT_L", "Rho25_L", "transition_P", "transition_S"))), moyenne=rep(0,length(x)), moyenne2=rep(0,length(x)))
# if (!is.na(parssm["DHA"])) ess$Parametre[,"DHA"] <- rnorm(replicate.CI,parssm["DHA"], res["DHA"])
# if (!is.na(parssm["Rho25"])) ess$Parametre[,"Rho25"] <- rnorm(replicate.CI,parssm["Rho25"], res["Rho25"])
# if (!is.na(parssm["DHL"])) ess$Parametre[,"DHL"] <- rnorm(replicate.CI,parssm["DHL"], res["DHL"])
# if (!is.na(parssm["DHH"])) ess$Parametre[,"DHH"] <- rnorm(replicate.CI,parssm["DHH"], res["DHH"])
# if (!is.na(parssm["T12L"])) ess$Parametre[,"T12L"] <- rnorm(replicate.CI,parssm["T12L"], res["T12L"])
# if (!is.na(parssm["T12H"])) ess$Parametre[,"T12H"] <- rnorm(replicate.CI,parssm["T12H"], res["T12H"])
# if (!is.na(parssm["DT"])) ess$Parametre[,"DT"] <- rnorm(replicate.CI,parssm["DT"], res["DT"])
# if (!is.na(parssm["DHA_L"])) ess$Parametre[,"DHA_L"]=rnorm(replicate.CI,parssm["DHA_L"], res["DHA_L"])
# if (!is.na(parssm["Rho25_L"])) ess$Parametre[,"Rho25_L"]=rnorm(replicate.CI,parssm["Rho25_L"], res["Rho25_L"])
# if (!is.na(parssm["DHL_L"])) ess$Parametre[,"DHL_L"]=rnorm(replicate.CI,parssm["DHL_L"], res["DHL_L"])
# if (!is.na(parssm["DHH_L"])) ess$Parametre[,"DHH_L"]=rnorm(replicate.CI,parssm["DHH_L"], res["DHH_L"])
# if (!is.na(parssm["T12L_L"])) ess$Parametre[,"T12L_L"]=rnorm(replicate.CI,parssm["T12L_L"], res["T12L_L"])
# if (!is.na(parssm["T12H_L"])) ess$Parametre[,"T12H_L"]=rnorm(replicate.CI,parssm["T12H_L"], res["T12H_L"])
# if (!is.na(parssm["DT_L"])) ess$Parametre[,"DT_L"]=rnorm(replicate.CI,parssm["DT_L"], res["DT_L"])
# if (!is.na(parssm["transition_P"])) ess$Parametre[,"transition_P"]=rnorm(replicate.CI,parssm["transition_P"], res["transition_P"])
# if (!is.na(parssm["transition_S"])) ess$Parametre[,"transition_S"]=rnorm(replicate.CI,parssm["transition_S"], res["transition_S"])


# 8/2/2014 dans parssm j'ai les paramtres
ess <- list(Parametre=matrix(rep(NA, length(parssm)*replicate.CI), ncol=length(parssm), dimnames=list(NULL, names(parssm))), moyenne=rep(0,length(x)), moyenne2=rep(0,length(x)))


for (i in 1:length(parssm)) {
  if (!is.na(res[names(parssm[i])])) {
    ess$Parametre[,i]=rnorm(replicate.CI, parssm[i], res[names(parssm[i])])
  } else {
    ess$Parametre[,i]=rep(parssm[i], replicate.CI)
  }
}



afficheCI <- TRUE

for (i in 1:replicate.CI) {
#  valeurlist <- embryogrowth:::.SSM(x+273.15, ess$Parametre[i,])
valeurlist <- .SSM(x+273.15, ess$Parametre[i,])
		
	if (!is.na(size) & !is.na(parssm["transition_S"]) & !is.na(parssm["transition_P"])) {
		r <- valeurlist[[1]]
		r_L <- valeurlist[[2]]
		transition <- 1/(1+exp(parssm["transition_S"]*(size-parssm["transition_P"])))
		valeur <- r*transition+r_L*(1-transition)

	} else {
		valeur <- valeurlist[[set.par[[rs]]]]
	}

	
	ess$moyenne <- ess$moyenne+valeur
	ess$moyenne2 <- ess$moyenne2+valeur^2
}

sdR=sqrt(ess$moyenne2/replicate.CI-(ess$moyenne/replicate.CI)^2)*scaleY


par(new=TRUE)
plot(x, y-2*sdR, type="l", col=col[[rs]], xlab="", ylab="", xlim=xlim, ylim=c(y1, y2), bty="n", axes = FALSE, lty=ltyCI[[rs]], lwd=lwdCI[[rs]], main="")
par(new=TRUE)
plot(x, y+2*sdR, type="l", col=col[[rs]], xlab="", ylab="", xlim=xlim, ylim=c(y1, y2), bty="n", axes = FALSE, lty=ltyCI[[rs]], lwd=lwdCI[[rs]], main="")

}
}

par(new=TRUE)
premier <- FALSE

# fin de la boucle des resultats
}

if (show.box) {
if (afficheCI) {
	legend(local.box, c("Mean", "Confidence interval"), lty=c(lty[[1]], ltyCI[[1]]), lwd=c(lwd[[1]], lwdCI[[1]]), bty = "n")
} else {
	legend(local.box, c("Mean"), lty=lty[[1]], lwd=lwd[[1]], bty = "n")
}
}


if (any(!is.na(unlist(legend)))) {
  legend("bottomright", unlist(legend), lty=unlist(lty), lwd=unlist(lwd), bty = "n", col=unlist(col))
}


}
