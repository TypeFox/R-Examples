##==============================================================================
##                                                                            ##
## LINEAR INVERSE MODELLING INPUT   -   LIM                                   ##
##    Karline Soetaert                                                        ##
##                                                                            ##
##==============================================================================


PrintMat <- function(lim)  {
  A           <- lim$A
  colnames(A) <- lim$Unknowns
  rownames(A) <- lim$eqnames
  G           <- lim$G
  colnames(G) <- lim$Unknowns
  rownames(G) <- lim$ineqnames
  print("A,B")
  print(cbind(A,"  B = "=lim$B))

  print("G,H")
  print(cbind(G,"  H = "=lim$H))

  if (!is.null(lim$Cost)) {
    Cost           <-  matrix(ncol = lim$NUnknowns, data=lim$Cost)
    colnames(Cost) <- lim$Unknowns
    rownames(Cost) <- lim$costnames
    print("Cost")
    print(Cost)
  }

  if (!is.null(lim$Profit)) {
    Profit           <-  matrix(ncol = lim$NUnknowns, data=lim$Profit)
    colnames(Profit) <- lim$Unknowns
    rownames(Profit) <- lim$profitnames
    print("Profit")
    print(Profit)
  }
}


