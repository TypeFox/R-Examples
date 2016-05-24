corReorder <-
function (R=mycor, vars=NULL, first=0,
          heat.map=TRUE, main=NULL, bottom=3,right=3,
          pdf.file=NULL, pdf.width=5, pdf.height=5) {


  # cor matrix:  mycor as class out_all, mycor$cors, or stand-alone matrix
  cor.nm <- deparse(substitute(R))
  .cor.exists(cor.nm)  # see if matrix exists in one of the 3 locations
  if (class(R) == "out_all")
    R <- eval(parse(text=paste(cor.nm, "$cors", sep="")))  # go to $cors 


  # translate variable names into column positions
  vars.all <- as.list(seq_along(as.data.frame(R)))
  names(vars.all) <- names(as.data.frame(R))
  vars.num <- eval(substitute(vars), vars.all, parent.frame())

  NVOld <- nrow(R)
  NVC <- as.integer(NVOld)

  Label <- integer(length=NVOld)
  Label <- as.integer(as.vector(vars.num))

  # -----------------------------
  # if not specified get vars.num 
  if (is.null(vars.num)) {


    IFirst <- as.integer(first)

    NV <- NVC

  # If IFirst = 0 (default), best variable is chosen first
  # If IFirst gt 0, IFirst is the first variable chosen by user

    Label <- rep(0, NV)

  # Find max sum sq variable for IFirst if not user specified
  # IFirst gives starting value
    if (IFirst == 0) {
      VMax <- 0.0
      for (I in 1:NV) {
        VT <- 0.0
        for (J in 1:NV) {
          VT <- VT + R[I,J]**2
        }
        if (VT > VMax) {
          VMax <- VT
          IMax <- I
        }
      }
      IFirst <- IMax
    }

  # get Labels

    Label[1] <- IFirst * 1000
    Label[IFirst] <- Label[IFirst] + 1

    for (I in 2:NV) {
      RMax <- 0.0
      K <- Label[I-1] / 1000
      for (J in 1:NV) {
        Lbctad <- Label[J] - (Label[J] %/% 1000) * 1000
        if (Lbctad == 0) {
          if (abs(R[K,J]) > RMax) {
            RMax <- abs(R[K,J])
            IMax <- J
          }
        }
      }
      Label[IMax] <- Label[IMax] + I
      Label[I] <- Label[I]+1000 * IMax
    }

    for (I in 1:NV) {
      Label[I] <- Label[I] %/% 1000
    }

    #vars.num <- Label  # derived vars.num (not user specified)
  }


  # -----------------------------
  # re-order R matrix
  R <- R[Label,Label]

  if (heat.map) {
   if (is.null(main)) main <- "Reordered Item Coefficients"
   .corcolors(R, nrow(R), main, bottom, right, diag=0,
              pdf.file, pdf.width, pdf.height)
  }

  # finish
  cat("\n")
  return(R)
}





