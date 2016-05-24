flatten.xbalresult <- function (x,
                                show.signif.stars=getOption("show.signif.stars"),
                                show.pvals=!show.signif.stars,...) {
  ##Notes: right now we've decided that you can't print both signif stars and p-values. make a choice.
  theresults <- x$results
  thestrata <- dimnames(theresults)[["strata"]]
  thestats <- dimnames(theresults)[["stat"]]
  thevars <- dimnames(theresults)[["vars"]]

  DIGITS <- max(2, getOption("digits")-4)

  if ("overall" %in% names(x)) { ##Extract the omnibus chisquared test from the xbal object...
    theoverall <- x$overall
  } else {
    theoverall<-NULL ##..or set it to NULL if it does not exist.
  }

  if (length(theresults) == 0 & is.null(theoverall)) {
    stop("There is a problem. Probably all of the variables (",
         all.vars(formula(x)),
         ") are constants within strata. Or else there is some other problem, try debug(RItools:::xBalance) to see what might be going on.")
  }

  if (length(theresults)==0 & !is.null(theoverall)){##The user has requested only the omnibus test and not the tests for the individual variables
    theresults<-NULL
    thevartab<-NULL
  }

  signifier <- function(data) {
    symnum(data,
           corr = FALSE,
           na = FALSE,
           abbr.colnames=FALSE, ##from print.summary.lm
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "** ", "*  ", ".  ", "   "))
  }

  flattener <- function(x) {   ##Summarize the variable-by-variable output array as a flat contingency table
    row.vars <- match("vars", names(dimnames(x)))
    col.vars <- match(c("strata","stat"), names(dimnames(x)))
    dn <- dimnames(x)
    dx <- dim(x)
    names(dx) <- names(dn)
    y <- aperm(x, c(row.vars, rev(col.vars)))
    dim(y) <- c(prod(dx[row.vars]), prod(dx[col.vars]))
    dimnames(y) <- list(variables=dn[["vars"]],
                        statistics=rep(dn[["stat"]],
                            dx["strata"]))
    as.data.frame(y)
  }

  if (show.signif.stars & !show.pvals & !is.null(theresults)) {

    Signif <- signifier(theresults[,"p",thestrata,drop=FALSE])

    ##junk<-do.call(cbind,lapply(which.strata,function(x){cbind(as.data.frame(theresults[,,x])," "=format(Signif[,,x]))}))
    newresults <- array(dim=dim(theresults)+c(0,1,0),
                        dimnames=list(vars=dimnames(theresults)[["vars"]],
                            stat=c(dimnames(theresults)[["stat"]],"sig."),
                            strata=dimnames(theresults)[["strata"]]))

    newresults[,-grep("sig.", dimnames(newresults)[[2]]),] <- theresults


    thevartab <- flattener(
        newresults[thevars, c(thestats[!(thestats=="p")],"sig."), thestrata, drop=FALSE])
    thevartab[names(thevartab)=="sig."] <- format(Signif)
    names(thevartab)[names(thevartab)=="sig."] <- ""
  }

  if (show.pvals & ("p" %in% dimnames(theresults)[["stat"]]) & !is.null(theresults)) {
    thevartab <- flattener(
        theresults[thevars,thestats,thestrata,drop=FALSE])
  }

  if (!is.null(theoverall)) {
    nc <- length(theresults)/2
    latex.annotation <- NULL
    ## paste("\\\\ \\hline Overall",
    ##       paste("\\multicolumn{",nc,"}{c}{",preSig,"}"),
    ##       paste("\\multicolumn{",nc,"}{c}{",postSig,"}"),
    ##       sep=" & ")
    if (show.signif.stars) {
      ChiSignif <- signifier(theoverall[thestrata,"p.value"])

      theoveralltab <- cbind(format(theoverall[thestrata,],digits=DIGITS),format(ChiSignif))
      names(theoveralltab)[4]<-" "
    }
    theoveralltab<-format(theoverall[thestrata,],digits=DIGITS)
  } else {
    theoveralltab<-NULL
  }

  theprintresults <- list(vartable=thevartab,overalltable=theoveralltab)
  latex.annotation <- paste("  \\multicolumn{",length(thestats),
                            "}{c}{",thestrata,"}")
  latex.annotation <- paste(latex.annotation,  collapse=" & ")
  theclinestarts <- seq(from=2, by=length(thestats), length.out=length(thestrata))
  theclinestops <-seq(from=1+length(thestats),
                      by=length(thestats), length.out=length(thestrata))

  theclines <- paste("\\cline{",theclinestarts,"-", theclinestops,"}",sep="")
  theclines <- paste(theclines, collapse=" ")
  latex.annotation <- paste("\\hline \\multicolumn{1}{r}{Strata:} &",
                            latex.annotation," \\\\",
                            theclines)

  structure(theprintresults,latex.annotation=latex.annotation)
}
