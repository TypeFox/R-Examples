regPlot <- 
function(out, type, digits.d=NULL, pred.intervals=TRUE,
         res.sort=c("cooks","rstudent","dffits","off"),
         res.rows=NULL, cooks.cut=1, scatter.coef=NULL,
         pdf=FALSE, pdf.width=5, pdf.height=5, manage.gr=FALSE, ...) {


  if (options("device") != "RStudioGD"  &&  is.null(options()$knitr.in.progress)) 
    .graphwin(1)


    # need: out$residuals
    if (type == 2) .reg3dnResidual(out, pdf, pdf.width, pdf.height, manage.gr)


    # need: out$fitted.values, out$residuals
    if (type == 3) .reg3resfitResidual(out, out$cooks.distance, cooks.cut,
                   pdf=FALSE, pdf.width=5, pdf.height=5, manage.gr=FALSE)


    # need: out$model (the data)
    if (type == 1) {

    #nm <- all.vars(out$terms)  # names of vars in the model
    nm <- all.vars(out$formula)
    n.vars <- length(nm)
    n.pred <- n.vars - 1

    in.data.frame <- TRUE
    for (i in 1:n.vars) {
      if (!(nm[i] %in% names(out$model))) {
        cat("\n\n\n>>> Note: ", nm[i], "is not in the data frame.\n")
        in.data.frame <- FALSE
      }
    }
    
    # check for all numeric vars  in.data.frame <- TRUE
    numeric.all <- TRUE
    for (i in 1:n.vars) {
      if (in.data.frame && 
        !is.numeric(out$model[1,which(names(out$model) == nm[i])])) {
          cat("\n\n\n>>> Note: ", nm[i], "is not a numeric variable.\n")
          numeric.all <- FALSE
        }
      }

    if (pred.intervals) {
      data <- out$model
      if (n.pred == 1) {
        o <- order(data[,nm[2]], decreasing=FALSE)
        data <- data[o,]
      }
      r2.out <- lm(out$formula, data=data)
      c.int <- data.frame(predict(r2.out, interval="confidence"))
      p.int <- data.frame(suppressWarnings(predict(r2.out, interval="prediction")))
    }
    else {
      c.int <- NULL
      p.int <- NULL
      pred.rows <- 0
    }

    if (!is.null(scatter.coef))
      s.c <- scatter.coef
    else
      s.c <- TRUE


    .reg5Plot(out, res.rows=NULL, pred.rows=NULL,
         scatter.coef=s.c, scatter.3D=FALSE, X1.new=NULL,
         numeric.all, in.data.frame, c.int, p.int,
         pdf=FALSE, pdf.width=5, pdf.height=5, manage.gr=FALSE, ...)

    } # end type == 1

}
