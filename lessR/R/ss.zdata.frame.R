.ss.data.frame <-
function(x, n.cat, brief, ...)  {

  for (i in 1:ncol(x)) {
    cat("\n\n")

    nu <- length(unique(na.omit(x[,i])))

    x.name <- names(x)[i]
    options(xname = x.name)

    if (is.numeric(x[,i]) && nu > n.cat) {
      stuff <- .ss.numeric(x[,i], brief=brief, ...)
      txsts <- stuff$tx
      txotl <- .outliers(x[,i])
      class(txsts) <- "out_piece"
      class(txotl) <- "out_piece"
      output <- list(out_stats=txsts, out_outliers=txotl)
      class(output) <- "out_all"
      print(output)
    }

    else if (is.factor(x[,i]) || is.character(x[,i]) ||
             (.is.num.cat(x[,i], n.cat))) {
      gl <- .getlabels(xlab=NULL, ylab=NULL, main=NULL, cex.lab=NULL)
      x.name <- gl$xn; x.lab <- gl$xb; x.lbl <- gl$xl
      stats <- .ss.factor(x[,i], x.name=x.name, x.lbl=x.lbl, ...)
      txttl <- stats$title
      counts <- stats$counts
      chi <- stats$chi
      class(txttl) <- "out_piece"
      class(counts) <- "out_piece"
      class(chi) <- "out_piece"
      output <- list(out_title=txttl, out_counts=counts, out_chi=chi)
      class(output) <- "out_all"
      print(output)      

      if (is.numeric(x[,i]) && nu <= n.cat)
        cat("\n>>> Variable is numeric, but only has", nu, "<= n.cat =", n.cat, "levels,",
        "so treat as categorical.\n",
        "   To obtain the numeric summary, decrease  n.cat  to indicate a lower\n",
        "   number of unique values such as with function: set.\n", 
        "   Perhaps make this variable a factor with the R factor function.\n")
    }

    else cat("\n>>> The following type of variable not processed: ", 
             class(x[,i]), "\n")

  }

}
