ANOVA <-
function(my.formula, data=mydata, brief=getOption("brief"), digits.d=NULL, 
         Rmd=NULL, graphics=TRUE,
         rb.points=TRUE, res.rows=NULL, res.sort=c("zresid", "fitted", "off"),
         pdf=FALSE, pdf.width=5, pdf.height=5, fun.call=NULL, ...) {  


  if (is.null(fun.call)) fun.call <- match.call()

  res.sort <- match.arg(res.sort)

  if (missing(my.formula)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify a model by listing it first or set according to:  my.formula\n\n")
  }

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      if (names(dots)[i] == "knitr.file") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "knitr.file  no longer used\n",
          "Instead use  Rmd  for R Markdown file\n\n")
      }
    }
  }

  dname <- deparse(substitute(data))
  options(dname = dname)
 
  op <- options()  # save current options to reset at end

  if (!exists(dname)) {
    txtC <- "Function ANOVA requires the data exist in a data frame\n"
    if (dname == "mydata") 
      txtA <- ", the default data frame name, " else txtA <- " "
    txtB1 <- "Either create the data frame, such as with data.frame function, or\n"
    txtB2 <- "  specify the actual data frame with the parameter: data\n"
    txtB <- paste(txtB1, txtB2, sep="")
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        txtC, "Data frame ", dname, txtA, "does not exist\n\n", txtB, "\n")
  }

  nm <- all.vars(my.formula)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1
  n.obs <- nrow(data)

  in.data.frame <- TRUE
  for (i in 1:n.vars) {
    if (!(nm[i] %in% names(data))) {
      cat("\n\n\n>>> Note: ", nm[i], "is not in the data frame.\n")
      in.data.frame <- FALSE
    }
  }

  for (i in 2:n.vars) {  # all IVs must be factors
      nms <-  which(names(data) == nm[i])
      if (in.data.frame && !is.factor(data[ , nms])) {
        cat("\n>>> Note: Converting", nm[i], "to a factor for this analysis only.\n")
        data[ ,nms] <- as.factor(data[ ,nms])
      }
    }  

  # ANOVA
  #   all analysis done on data in model construct av.out$model
  #   this model construct contains only model vars, with Y listed first
  #assign("av.out", aov(my.formula, data=data), pos=.GlobalEnv)
  av.out <- aov(my.formula, data=data)

  n.keep <- nrow(av.out$model)

  if (is.null(digits.d)) digits.d <- .getdigits(data[,nm[1]], 2)
    

# ----------
# Background
# ----------

  title_bck <- "  BACKGROUND"

  tx <- character(length=0)

  if (sys.nframe() == 1) {  # only accurate if not called from model
    tx[length(tx)+1] <- paste("Data Frame: ", dname)
    tx[length(tx)+1] <- ""
  }
  
  for (i in 1:n.vars) {
    ind <- i
    tx2 <- .varlist2(n.pred, ind, nm[i], "Factor", n.obs,
                     n.keep, levels(data[,nm[i]]))
    for (j in 1:length(tx2)) tx[length(tx)+1] <- tx2[j]
  }

  if (n.pred == 2) {
    if (!is.list(replications(my.formula, data=data))) {
      tx[length(tx)+1] <- ""
      tx[length(tx)+1] <- "The design is balanced"
    }
    else {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The design is not balanced. The results would be invalid.\n",
        "Consider function  lmer  in the  lme4  package.\n\n")
    }
  }

  txbck <- tx


# --------
# Analysis
# --------

  # keep track of generated graphics
  if (graphics) {
    plot.i <- 0L
    plot.title  <- character(length=0)
    manage.gr <- .graphman()
  }

  if (n.pred == 1)  {
    plt1 <- .ANOVAz1(av.out, av.out$model[,nm[1]], av.out$model[,nm[2]],
        nm, n.obs, digits.d, brief, graphics, pdf, pdf.width, pdf.height)
    title_des <- plt1$title_des
    txdes <- plt1$txdes
    title_basic <- plt1$title_basic
    txanv <- plt1$txanv
    txeft <- plt1$txeft
    txhsd <- plt1$txhsd
    if (graphics) {
      for (i in (plot.i+1):(plot.i+plt1$i)) plot.title[i] <- plt1$ttl[i-plot.i]
      plot.i <- plot.i + plt1$i
    }
  }

  if (n.pred == 2) {
    plt2 <- .ANOVAz2(av.out, av.out$model[,nm[1]], av.out$model[,nm[2]],
        av.out$model[,nm[3]], nm, digits.d, brief, as.character(my.formula)[3],
        rb.points, graphics, pdf, pdf.width, pdf.height)
    txbck2 <- plt2$txbck2
    for (i in 1:length(txbck2)) tx[length(txbck)+1] <- txbck2[i]
    title_des <- plt2$title_des
    txb2 <- plt2$txbck2
    for (i in 1:length(txb2)) txbck[length(txbck)+1] <- txb2[i]
    txcn <- plt2$txcn
    txcm <- plt2$txcm
    txmm <- plt2$txmm
    txgm <- plt2$txgm
    txcs <- plt2$txcs
    title_basic <- plt2$title_basic
    txanv <- plt2$txanv
    txeft <- plt2$txeft
    txhsd <- plt2$txhsd
    if (graphics) {
      for (i in (plot.i+1):(plot.i+plt2$i)) plot.title[i] <- plt2$ttl[i-plot.i]
      plot.i <- plot.i + plt2$i
    }
  }

  # residuals
  if (!brief) {
    tx <- character(length=0)

    n.keep <- nrow(av.out$model)
    if (is.null(res.rows)) if (n.keep < 20) res.rows <- n.keep else res.rows <- 20 
    if (res.rows == "all") res.rows <- n.keep  # turn off resids with res.rows=0

    title_res <- "  RESIDUALS" 
    tx[length(tx)+1] <- "Fitted Values, Residuals, Standardized Residuals"
    if (res.sort == "zresid")
      tx[length(tx)+1] <- "   [sorted by Standardized Residuals, ignoring + or - sign]"
    if (res.sort == "fitted")  
      tx[length(tx)+1] <- "   [sorted by Fitted Value, ignoring + or - sign]"
    if (res.rows < n.keep)
      txt <- "cases (rows) of data, or res.rows=\"all\"]"
    else
      txt="]"
    tx[length(tx)+1] <- paste("   [res.rows = ", res.rows, ", out of ", n.keep, " ", txt, sep="")

    fit <- fitted(av.out)
    res <- residuals(av.out)
    sres <- rstandard(av.out)
    out <- cbind(av.out$model[c(nm[seq(2,n.vars)],nm[1])], fit, res, sres)
    out <- data.frame(out)
    names(out)[n.vars+1] <- "fitted"
    names(out)[n.vars+2] <- "residual"
    names(out)[n.vars+3] <- "zresid"

    if (res.sort != "off") {
      if (res.sort == "zresid") o <- order(abs(out$zresid), decreasing=TRUE)
      if (res.sort == "fitted") o <- order(abs(out$fitted), decreasing=TRUE)
      out <- out[o,]
    }
    names(out)[n.vars+3] <- "z-resid"
    for (i in 1:(n.vars+3))
      if (is.numeric(out[,i])) if (!is.integer(out[,i])) {
        if (digits.d > 2) dec.digits <- digits.d-1 else dec.digits <- digits.d
        out[,i] <- .fmt(out[,i],dec.digits)
      }
    tx2 <- .prntbl(out[1:res.rows,], digits.d)
    for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]

    txres <- tx
  }

# pairwise.t.test(mydata$Steady, mydata$TrtA)
# power.anova.test(groups=4, n=8, between.var=16.33, within.var=2.179)

  options(op)  # restore options going into reg

  # display list of plots if more than 1
  txplt <- ""
  if (graphics) {
    if (plot.i > 1) txplt <- .plotList2(plot.i, plot.title)
    dev.set(which=2)  # reset graphics window for standard R functions
  }

  # R Markdown
  txkfl <- ""
  if (!is.null(Rmd)) {
    txt <- ifelse (grepl(".Rmd", Rmd), "", ".Rmd")
    Rmd <- paste(Rmd, txt, sep="") 
    txknt <- .av.Rmd(nm, dname, fun.call, digits.d)
    cat(txknt, file=Rmd, sep="\n")
    txkfl <- .showfile2(Rmd, "R Markdown instructions")
  }

  class(title_bck) <- "out_piece"
  class(txbck) <- "out_piece"
  class(title_des) <- "out_piece"
  if (n.pred == 1)
    class(txdes) <- "out_piece"
  if (n.pred == 2) {
    class(txcn) <- "out_piece"
    class(txcm) <- "out_piece"
    class(txmm) <- "out_piece"
    class(txcs) <- "out_piece"
  }
  class(title_basic) <- "out_piece"
  class(txanv) <- "out_piece"
  class(txeft) <- "out_piece"
  class(txhsd) <- "out_piece"
  class(title_res) <- "out_piece"
  class(txres) <- "out_piece"
  class(txplt) <- "out_piece"
  

  if (n.pred == 1)  {
    output <- list(
      call=fun.call, formula=my.formula,

      out_title_bck=title_bck, out_background=txbck,

      out_title_des=title_des,
      out_descriptive=txdes,

      out_title_basic=title_basic,
      out_anova=txanv, out_effects=txeft, out_hsd=txhsd, 

      out_title_res=title_res, out_res=txres,

      out_plots=txplt,

      n.vars=n.vars, n.obs=n.obs, n.keep=n.keep,
      residuals=res, fitted=fit
    )
  }

  if (n.pred == 2)  {
    output <- list(
      call=fun.call, formula=my.formula,

      out_title_bck=title_bck, out_background=txbck,

      out_title_des=title_des,
      out_cell.n=txcn, out_cell.means=txcm, out_marginals=txmm,
      out_gm=txgm, out_cell.sd=txcs,

      out_title_basic=title_basic,
      out_anova=txanv, out_effects=txeft, out_hsd=txhsd, 

      out_title_res=title_res, out_res=txres,

      out_plots=txplt,

      n.vars=n.vars, n.obs=n.obs, n.keep=n.keep,
      residuals=res, fitted=fit
    )
  }

  class(output) <- "out_all"

  return(output)

}
