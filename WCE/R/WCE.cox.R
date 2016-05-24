 WCE.cox  <- function(data, nknots, cutoff, constrained = FALSE, int.knots = NULL, aic = FALSE, id='Id', event = 'Event',  start='Start', stop='Stop', expos ='dose', covariates = NULL, controls = NULL,  ...) {
   # rename variables for convenience
   names(data)[names(data) == id] = 'Id'
   names(data)[names(data) == event] = 'Event'
   names(data)[names(data) == start] = 'Start'
   names(data)[names(data) == stop] = 'Stop'
   names(data)[names(data) == expos] = 'dose'
   if (sum(c("Id","Event","Start","Stop","dose") %in% names(data))!= 5)  stop("ERROR: At least one of 'id','event','start','stop', or 'expos' is missing from the dataset supplied")
   maxTime <- max(daply(data, "Id", function(df).maxfu(df)))
   init <- 1

	nknots <- sort(unique(nknots))

   if (is.null(int.knots) == F & length(nknots) >1) {
   init <- length(int.knots)
	cat('Only 1 model with the interior knots specified will be estimated. \n\n')}
   if (cutoff > maxTime)  stop("ERROR: cutoff must be smaller than the longest follow-up time")
   if (constrained != FALSE & !constrained %in% c(FALSE, 'Right', 'R', 'right', 'Left', 'L', 'left')) stop ("ERROR: constrained has to be one of : FALSE, 'Right', or 'Left'.")
	if (constrained %in% c('R', 'right')) {constrained = 'Right'}
	if (constrained %in% c('L', 'left')) {constrained = 'Left'}
  PL.c <- rep(0, length(nknots))
  coef.mat.c <- list()
vcov.mat.c <- list()


  sed <- list()
  gnu.c <-  matrix(0, cutoff, nrow=length(nknots))  
  BIC.c <- rep(0, length(nknots))
  kev <- list()
  if (is.null(covariates) == F){
    covbeta <- matrix(0, ncol = length(covariates), nrow = length(nknots))
    covrobse <- matrix(0, ncol = length(covariates), nrow = length(nknots))}
  n.events <-   sum(data$Event)
  ev <- sort(unique(data[data$Event==1,]$Stop))

  for (j in 1: length(nknots)){  
 
 if (is.null(int.knots) == F) {
	kev[[paste(nknots[j],"knot(s)", sep=" ")]] <- .augm.knots(int.knots, cutoff)} else {
	kev[[paste(nknots[j],"knot(s)", sep=" ")]] <- .augm.knots(.knots.equi(nknots[j], cutoff), cutoff)
	}



  Bbasis <- splineDesign(knots = kev[[j]], x = 1:cutoff, ord=4)
  smalldat <- data[data$Stop %in% ev,]
  Id <- unique(data$Id)
  
  kal <- data.frame(do.call("rbind", lapply(1:length(Id), function(i) .wcecalc(ev, data$dose[data$Id==Id[i]],data$Stop[data$Id==Id[i]],Bbasis, cutoff))))
  kal <- kal[is.na(kal[,1])==FALSE,]
  names(kal) <- paste("D", 1:dim(kal)[2], sep="")
  smalldat <- cbind(smalldat, kal)

  if (is.null(controls) == T) { controls = coxph.control(1e-09, .Machine$double.eps^0.75, 20, sqrt(1e-09), 10)}

  if (constrained == 'Right'){
	co <- tryCatch(.EstimateSplineConstrainedC(smalldat, nknots[j], covariates, 'Right', controls), error=function(e) NULL)}
  if (constrained == 'Left'){
	co <- tryCatch(.EstimateSplineConstrainedC(smalldat, nknots[j], covariates, 'Left', controls), error=function(e) NULL)}
  if (constrained == F){
	co <- tryCatch(.EstimateSplineUnconstrainedC(smalldat, nknots[j], covariates, controls), error=function(e) NULL)}
  if (is.null(co)==FALSE) {
    PL.c[j] <- co$ll[2]
if (constrained == 'Left') { gnu.c[j,] <-  Bbasis[,3:(2+length(co$Dvar))] %*% co$coefs[c(co$Dvar)]} else {
    gnu.c[j,] <-  Bbasis[,1:length(co$Dvar)] %*% co$coefs[c(co$Dvar)]} # 
    coef.mat.c[[paste(nknots[j],"knot(s)", sep=" ")]] <- co$coefs[c(co$Dvar)]
		vcov.mat.c[[paste(nknots[j],"knot(s)", sep=" ")]] <- co$vcovmat
	sed[[j]] <- unlist(co$SE[c(co$Dvar)])
    if (is.null(covariates) == F){
      covbeta[j,1:length(covariates)] <- co$coefs[c(covariates)]
      covrobse[j,1:length(covariates)] <- unlist(co$SE[c(covariates)])} 
    BIC.c[j] <- .my.bic.c(PL.c[j], n.events ,nknots[j], constrained, aic, covariates)}
  co <- NULL
}

# rename to ease readability
rownames(gnu.c) <- paste(nknots, 'knot(s)')
colnames(gnu.c) <- paste('t',1:ncol(gnu.c), sep='')
PL.c <- matrix(PL.c, nrow=1)
colnames(PL.c) <- paste(nknots, 'knot(s)')
 if (is.null(covariates) == F){
	rownames(covbeta) <- paste(nknots, 'knot(s)')
	colnames(covbeta) <- covariates
	rownames(covrobse)<- paste(nknots, 'knot(s)')
	colnames(covrobse) <- covariates} else {covbeta <- covrobse <- NA}
BIC.c <- matrix(BIC.c, nrow=1)
colnames(BIC.c) <- paste(nknots, 'knot(s)')
   if (is.null(covariates)){
     est <- list(knotsmat = kev, WCEmat = gnu.c,  PL = PL.c, est = coef.mat.c, vcovmat = vcov.mat.c, constrained = constrained, covariates = NULL, ne = n.events, a = aic, info.criterion = BIC.c, analysis = 'Cox')} else {
     est <- list(knotsmat = kev, WCEmat = gnu.c,  PL = PL.c, est = coef.mat.c, vcovmat = vcov.mat.c, SED = sed, beta.hat.covariates=covbeta, se.covariates= covrobse, covariates = covariates, constrained = constrained, ne = n.events, a = aic, info.criterion = BIC.c, analysis = 'Cox')}
   class(est) <- "WCE"
   return(est)
 }