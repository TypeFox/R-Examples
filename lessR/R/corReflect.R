corReflect <- 
function (R=mycor, vars,
          main=NULL, heat.map=TRUE, bottom=3,right=3, 
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

  NVOld <- as.integer(nrow(R))
  NVC <- as.integer(length(vars.num))

  Label <- as.integer(as.vector(vars.num))

  for (LL in 1:NVC) {
    for (J in 1:NVOld) {
      if (Label[LL] != J) {
        R[J,Label[LL]] <- -R[J,Label[LL]]
        R[Label[LL],J] <- R[J,Label[LL]]
      }
    }
  }

  if (heat.map) {
    if (is.null(main)) main <- "With Reflected Item Coefficients"
   .corcolors(R, NVOld, main, bottom, right, diag=0,
              pdf.file, pdf.width, pdf.height)
  }

  cat("\n")
  return(R)
}
