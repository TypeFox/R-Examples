##' Make a scatter plot iwth histograms
##'
##' From an example in the R documentation
##'
##' @param x x variable
##' @param y y variable
##' @param hist.col color of histogram
##' @param trend.line how to add trend line
##' @param ... passed on
##' @return NULL
##' @export
scatter.with.hist <-
  function(x,y,hist.col=gray(.95),trend.line="lm",...) {

    ## Make a scatterplot with trendline and
    ## histograms of each distribution.
    
    on.par <- par(no.readonly = TRUE)
    on.exit(par(on.par))                # see ?par for details

    nf <- layout(matrix(c(1,0,          # which order to place graphs
                          3,2),
                        2,2,byrow=TRUE),
                 widths=c(3,1),         # 3/4 wide for col. 1
                 heights=c(1,3),        # 3/4 wide for row 2
                 respect=TRUE)          # make square
    layout.show(nf)

    n<-length(x)
    no.breaks = max(nclass.scott(x),nclass.scott(y))
    xhist <- hist(x,breaks=no.breaks, plot=FALSE)
    yhist <- hist(y,breaks=no.breaks, plot=FALSE)
    top <- max(c(xhist$counts, yhist$counts))

    ## adjust margins for better look
    par(mar=c(0,3,1,1))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top),
            space=0,col=hist.col)

    par(mar=c(3,0,1,1))
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top),
            space=0,col=hist.col, horiz=TRUE)


    par(mar=c(4,4,1,1))
    x.name = deparse(substitute(x))
    y.name = deparse(substitute(y))
    plot(x,y,xlab=x.name,ylab=y.name,...)
    if(!is.null(trend.line) && !is.na(trend.line)) {
      switch(trend.line,
             "lm" = abline(lm(y~x)),
             "supsmu" = lines(supsmu(x,y)),
             "lowess" = lines(lowess(x,y)),
             NULL
             )
    }

    invisible()                         # restores par settings
  }
