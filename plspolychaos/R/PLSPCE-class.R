#########################################################
# Class to store the result of function 'PCESI'
#########################################################
setClass("PLSPCE", slots=c(
                      indexes="matrix", # SI
                     indexes.percent="matrix", # %SI 
                     R2="matrix",
                       Q2="matrix",
                       ncopt="numeric",
                       rmsep="matrix",
                       y.hat="vector",
                       COEF="matrix", # natural beta for all the components
                     betaCR="matrix", # CR beta for all the components
                     STRUC="PCEdesign" # la structure du polynome,
                     ))

#########################################################
# print method
# "all": option to fix the display. Vector. Valid values are:
# "TRUE" : rmsep, Beta coefficients and y.hat are displayed
# " ..." : all options passed as it to print
#########################################################

print.PLSPCE <- function (x, all=FALSE, ...) {

    cat("\nExplanation level of the response (R2, percentage and cumulated percentage)\n")
    print(round(x@R2,4), ...)
    
    cat("\nExplanation-prediction level of the response (Q2 and Q2cum)\n")
    print(round(x@Q2, 4), ...)     

    cat("\nOptimal number of components: ", x@ncopt, "\n")
#############################################
# Afficher les resultats correspondants a la composante optimale
#############################################
    cat("\nExplanation level of the optimal number of components\n")
    R2opt <- matrix(x@R2[x@ncopt,], ncol=ncol(x@R2),
                dimnames=list(rownames(x@R2)[x@ncopt], colnames(x@R2))
                )

    print(round(R2opt,4), ...)   
    cat("\nExplanation-prediction level of the optimal number of components\n")
    Q2opt <- matrix(x@Q2[x@ncopt,], ncol=ncol(x@Q2),
                dimnames=list(rownames(x@Q2)[x@ncopt], colnames(x@Q2))
                )
    print(round(Q2opt, 4), ...)     

    rmsepopt <- matrix(x@rmsep[x@ncopt,], ncol=ncol(x@rmsep),
                dimnames=list(rownames(x@rmsep)[x@ncopt], colnames(x@rmsep)))

    cat("\nRoot Mean Square Prediction of the optimal number of components\n")
    print(round(rmsepopt,4), ...)

    cat("\nPLS-PCE sensivity indexes\n")
    print(round(x@indexes,4), ...)
    cat("\n%PLS-PCE sensivity indexes\n")
    print(round(x@indexes.percent, 4), ...)
    if (all) {
        cat("\n")
        print(x@STRUC, ...)
        cat("Number of rows:", length(x@y.hat), "\n")

        cat("\nAlso included:")
        cat("\n * slot 'COEF' (PLS-regression coefficients). Dimension: ")
        cat(dim(x@COEF))
        cat("\n * slot 'betaCR' (centered-reducted PLS-regression coefficients). Dimension: ")
        cat(dim(x@betaCR))
        cat("\n * slot 'y.hat' (metamodel outputs). Length: ")
        cat(length(x@y.hat))
        cat("\n * slot 'rmsep' (Root Mean Square Predictions). Length: ")
        cat(length(x@rmsep))
        cat("\n * slot 'STRUC' (matrix coding the polynomial expression). Dimension: ")
        cat(dim(x@STRUC@.Data), "\n")
    } # fin all

    return(invisible())
} # end print.PLSPCE

#########################################################
# show method
#########################################################
show.PLSPCE  <- function(object){
  print.PLSPCE(object)
  return(invisible())
} # end show.PLSPCE


setMethod("show", signature(object="PLSPCE"),
          definition=show.PLSPCE)

                
#########################################################
# plot
#########################################################
plot.PLSPCE <- function(x, pce, options =c("fit", "bar", "compo"), ...) {
  # x: return of calcPLSPCE (object PLSPCE)
  # pce: return of analyticsPolyLeg or polyLeg (object PCEpoly)

###   # Three plots
  Y <- pce@.Data[, ncol(pce@.Data), drop=FALSE] 
  nmono <- nrow(x@STRUC) -1 # -1, car il y a le terme cte
  degree <- x@STRUC@degree
  
## subtitle
  lesub <- paste("degree", degree, "-", nmono, "monomials",
            "-", length(Y), "rows -", pce@nvx, "inputs")
 
### Plot of the observed and fitted y
  if ("fit" %in% options) {
  y.hat <-x@y.hat
  reg <- lm(y.hat ~Y)
  plot(y.hat, Y,
       xlab="metamodel output", ylab="computer model output",
       sub=lesub,
       main="Scatter plot and regression line",
       cex.main=0.9)
  
  lines(reg$fitted.values, Y)
} # fin fit

  
### Bar of the indexes
  if ("bar"  %in% options) {
    if ("fit" %in% options) {
  if (dev.interactive()) {
    dev.new()
  }
} # fin fit

    
    # Tri decroissant des TPE
    TPEtri <- sort(x@indexes[,"TPE"], decreasing = TRUE, index.return = TRUE)
    indicesbar <- x@indexes[TPEtri$ix, c("PE", "TPE")]
    # comme barplot accumule les valeurs des 2 barres
    # je soustrais PE de TPE
    indicesbar[, "TPE"] <- indicesbar[, "TPE"] - indicesbar[, "PE"]

    # Ecrire les noms des inputs verticalement s'ils sont longs
    # (>3 car) ou s'il y a bcp de barres (>15)
    if ( any(nchar(rownames(indicesbar)) >3) ||
        (pce@nvx >15)) {
      las <- 3
    } else {
      las <-0
    }
      
  barplot(t(indicesbar),  legend.text=TRUE,
          sub=lesub,
          main="Polynomial (PE) and Total (TPE) PLS-PCE sensivity indexes",
          cex.main=0.9, col=c("#4086A4", "#A44040"), las=las)

  } # fin bar
  

###    Plot the TPE SI of all the components
  if ("compo"  %in% options) {
    if ( ("fit" %in% options) || ("bar" %in% options)) {
      if (dev.interactive()) {
        dev.new()
      }
    } # fin ( ("fit" %in% options) || ("bar" %in% options))
    
  # Compute the PLS-PCE SI of all the components
  beta <- x@COEF
  nctot <- ncol(beta)

  # XM: inputs without Y but with the constant term
  XM <- pce@.Data[, -ncol(pce@.Data)]
  # data.exp: inputs without Y and without the constant term
  data.exp <- pce@.Data[, -c(1, ncol(pce@.Data)), drop=FALSE]
  mu.x <-apply(data.exp, 2, mean)
  sd.x <- apply(data.exp, 2, sd)
  sd.y <- sd(Y)
  indices <- matrix(NA, nrow=nctot, ncol=pce@nvx)
  
  for (nc in 1: nctot) {
    indices[ nc,] <- calcPLSPCESI(beta[,nc, drop=FALSE], data.exp, XM, Y,
                         mu.x,  sd.x,sd.y, pce )[, "TPE"]
  }
  colnames(indices)<- colnames(pce@STRUC)
  bornesy <- range(c(indices, x@Q2[, "Q2cum"]))

   plot(c(1,nctot), bornesy, type="n",
       xlab="components", ylab="PLS-PCE indexes (TPE)",
       main="Total PLS-PCE sensivity indexes against components",
       cex.main=0.9,
       sub=paste("Optimal number of components:", x@ncopt, "-",
         lesub))


    
  for (i in 1:pce@nvx) {
    lines(1:nctot, indices[,i], col=i)
        }

     # Rajouter les Q2cum
    lines(1:nctot, x@Q2[, "Q2cum"], lty=2)
    
  legend( "topright", c(colnames(indices), "Q2cum"),
         lty=c(rep(1, ncol(indices)),2),
           cex=0.8, col=c(1:pce@nvx,1), inset=c(0,0.1))
  # marquer la composante optimale
    lines(c(x@ncopt,x@ncopt),
          c(bornesy[1], bornesy[2]),
          col = "blue")

  axis(1, at=x@ncopt, col="blue")
  } # fin compo

  
  return(invisible())
}
