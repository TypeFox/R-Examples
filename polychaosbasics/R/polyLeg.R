###################################################
# polyLeg
# Build PCE Design from provided LHD
polyLeg <- function(lhs, Y, degree)
  {

   if (degree <= 1) {
     stop("The degree should be greater than 1")
   }
    
  nvx <- ncol(lhs)
     
  nlhs <- nrow(lhs)
  if (length(Y) !=  nlhs) {
    stop("lhs and Y should have the same number of rows")
  }
  
# Calibration du lhd de base sur la plage [-1,1]
# pour polynomes de Legendre
 binfl<-rep(-1,nvx) ; bsupl<-rep(1,nvx)
 lhdc<-calibDesign(lhs,binfl,bsupl)
  
#  Construction des monomes
  plan2<-Structure(nvx,degree)
  dimnames(plan2) <- list(c("0", labelmono(plan2)),
                          paste("Input", 1:nvx))

#  Construction de la matrice du modele 
  XM<-modLeg(lhdc,degree,plan2)

 # Return
 XMY <- cbind(XM,Y)
   retour <- new("PCEpoly", .Data=XMY, design=plan2, nvx=nvx,
                 call=match.call())
 return(retour)

}
###################################################
# labelmono
# labellelise les monomes

labelmono <- function(x) {
    # oter le terme constant
  planx <- x[-1,, drop=FALSE]

  label <- rep("", nrow(planx))
  for (i in 1:nrow(planx)) {
    descr <- ""
    prem <- FALSE
    for (j in 1:ncol(planx)) {
      if (planx[i,j] >0 ) {
        if (prem) {
          descr <- paste(descr, "*", sep="")
        }
        prem <- TRUE # la 1ere variable du monome est vue
        descr <- paste(descr, j, sep="")
        if (planx[i,j] >1) {
          for (k in 2:planx[i,j]) {
            descr <- paste(descr, "*", j, sep="")
          }
        } # fin (planx[i,j] >1)
      } # fin (planx[i,j] >0 )
    } # fin j
    label[i] <- descr
  } # fin i
return(label)
} # fin labelmono
