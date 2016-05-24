case <- function(fit, ...)
  UseMethod("case")

lm.case <- function(...)
  .Defunct("case.lm", package="HH")

case.lm <-  ## previously lm.case
function(fit, lms=summary.lm(fit), lmi=lm.influence(fit), ...) {
  e <- resid(fit)
  s <- lms$sigma
  xxi <- diag(lms$cov.unscaled)
  si <- lmi$sigma
  h <- lmi$hat
  bi <-  coef(lmi)
  sta.res <- e/(s*(1-h)^.5)
  stu.res <- e/(si*(1-h)^.5)  ## uses si, not s
  dfbetas <- bi/(si %o% xxi^.5)
  dimnames(dfbetas)[[2]] <- names(coef(fit))
  dffit <- h*e/(1-h)
  dffits <- h^.5*e/(si*(1-h))
  cook <- sta.res^2/length(xxi) * h/(1-h)
  tmp <- cbind(e, h, si, sta.res, stu.res, dffit, dffits, cook, dfbetas)
  class(tmp) <- c("case", class(tmp))
  tmp
}

## plot.case is based on:
## Section 4.3.3 Influence of Individual Obervations
## in Chambers and Hastie, Statistical Models in S.
plot.case <- function(x, fit,
                      which=c("stu.res","si","h","cook","dffits",
                        dimnames(x)[[2]][-(1:8)]),  ##DFBETAS
                      between.in=list(y=4, x=9),
                      cex.threshold=1.2,
                      main.in=list(
                        paste(deparse(fit$call), collapse=""),
                        cex=main.cex),
                      sigma.in=summary.lm(fit)$sigma,
                      p.in=summary.lm(fit)$df[1]-1,
##                      obs.large=".lm.case.large",
##                      obs.large.env=.GlobalEnv,
                      main.cex=NULL,
                      ...) {
  p <- p.in
  n <- dim(x)[1]

  ncs <- dimnames(x)[[2]]

  ncs.keep <- if (is.numeric(which)) which else match(which, ncs, 0)

  ncs[-(1:8)] <- paste("DFBETAS", ncs[-(1:8)])
  ncs[match("si",ncs)] <- "deleted std dev"
  ncs[match("cook",ncs)] <- "Cook's distance"
  ncs[match("stu.res",ncs)] <- "Student del resid"
  case.data.frame <-
    data.frame(y     = as.vector(unlist(x[,ncs.keep])),
               id    = rep(1:n, length(ncs.keep)),
               group = factor(
                 rep(ncs[ncs.keep], rep(n, length(ncs.keep))),
                 levels=ncs[ncs.keep]))
  case.data.frame$rownames <- rep(dimnames(x)[[1]],
                                  length(ncs.keep))
  ## 'rownames' is a character variable.  This is not 'row.names'.
  scales.y <- list(relation="free", limits=list())
  scales.y["at"] <- list(NULL)
  ## case.large <- list()
  ## class(case.large) <- "lm.case.large"
  ## assign(".lm.case.large", case.large, envir=environment())
  ## large.env <- environment(case.large)

  nn <- n
  pp <- p
  ss <- sigma.in
  panel.labels <- levels(case.data.frame$group)
  thresh <- rep(list(list( ## DFBETAS
    threshold=c(-2,0,2)/sqrt(n),
    thresh.label=expression(-2 / sqrt(n), 0, 2 / sqrt(n)),
    thresh.id=c(-2,2)/sqrt(n)
  )),
                length=length(panel.labels))
  names(thresh) <- panel.labels

  thresh$"Student del resid"=list(
           threshold=c(-3,0,3),
           thresh.label=c(-3,0,3),
           thresh.id=c(-3,3))

  thresh$"deleted std dev"=list(
           threshold=c(.95*ss,ss,1.05*ss)-ss,
           thresh.label=c(".95 s","s","1.05 s"),
           thresh.id=c(.95*ss, 1.05*ss))

  thresh$"h"=list(
           threshold=c(1/nn, c(2,3) * (pp+1)/nn),
           thresh.label=c("1 / n", "2 (p+1)/n", "3 (p+1)/n"),
           thresh.id=c(0, 2*(pp+1)/nn))

  thresh$"Cook's distance"=list(
           threshold=c(-1,0,1),
           thresh.label=c("","0","1"),
           thresh.id=c(-1,1))

  thresh$"dffits"=list(
           threshold=c(-2,0,2)*sqrt(pp/nn),
           thresh.label=expression(-2 ~ sqrt(p/n), 0, 2 ~ sqrt(p/n)),
           thresh.id=c(-2,2)*sqrt(pp/nn))

  ## the last three aren't currently displayed
  thresh$"e"=list(
           threshold=c(-3*ss,-2*ss,0,-ss,ss,2*ss,3*ss),
           thresh.label=c("-3s","-2s","0","-s","s","2s","3s"),
           thresh.id=c(-2*ss, 2*ss))

  thresh$"sta.res"=list(
           threshold=c(-3,0,3),
           thresh.label=c(-3,0,3),
           thresh.id=c(-3,3))

  thresh$"dffit"=list(
           threshold=c(-100,0,100), ## arbitrary large numbers
           thresh.label=c(-100,0,100),
           thresh.id=c(-100,100))

  case.large <- rep(list(), length=length(panel.labels))
  names(case.large) <- panel.labels
  for (i in panel.labels) {
    y <- case.data.frame[i==case.data.frame$group, "y"]
    y.compare <- (y < thresh[[i]]$thresh.id[1]) | (y > thresh[[i]]$thresh.id[2])
    y.compare[is.na(y.compare)] <- FALSE
    case.large[[i]] <- (1:nn)[y.compare]
    scales.y$limits[[i]] <- range(y, thresh[[i]]$threshold)
  }

  scales.y$limits[["deleted std dev"]] <- ## this needs adjusting
    range(case.data.frame["deleted std dev"==case.data.frame$group, "y"],
          thresh[["deleted std dev"]]$threshold + ss)

  result <-
  xyplot(y ~ id | group, data=case.data.frame,
         subscripts=TRUE,
         rownames=case.data.frame$rownames,
         group.names=levels(case.data.frame$group),
         panel=panel.case,
         main=main.in,
         xlab="",
         ylab="",
         between=between.in,
         as.table=TRUE,
         scales=list(
           y=scales.y,
           x=list(alternating=2, tck=-.005)),
         thresh=thresh, case.large=case.large,
         nn=n, pp=p, ss=sigma.in,
         cex.threshold=cex.threshold,
         par.settings = list(
           clip=list(panel="off"),
           layout.widths = list(axis.key.padding = 16)),
         ...)
  class(result$panel.args.common$case.large) <- "lm.case.large"
  class(result) <- c("trellis.case", class(result))
  result
}
## environment(plot.case) <- environment(plot.likert)


"panel.case" <- function(x, y, subscripts, rownames, group.names,
                         thresh, case.large, nn, pp, ss, cex.threshold,
                         ...) {

  cell.num <- panel.number()
  panel.label <- group.names[cell.num]
  cex.x <- trellis.par.get("axis.text")$cex
  cex.y <- cex.x

  y.plot <- y
  pretty.y <- pretty(y)
  lin.pretty.y <- pretty.y

  new.viewport <- FALSE
  switch(panel.label,
         "deleted std dev"={y.plot <- y-ss
                            lin.pretty.y <- pretty.y-ss
                            new.viewport <- TRUE
                            cvp <- current.viewport()
                            cvp$yscale <- cvp$yscale - ss
                            pushViewport(cvp)
                            ## threshold <- c(.95*ss,ss,1.05*ss)-ss
                            ## thresh.label <- c(".95 s","s","1.05 s")
                            ## thresh.id <- c(.95*ss, 1.05*ss)
                          },
         e={## threshold <- c(-3*ss,-2*ss,0,-ss,ss,2*ss,3*ss)
            ## thresh.label <- c("-3s","-2s","0","-s","s","2s","3s")
            ## thresh.id <- c(-2*ss, 2*ss)
          },
         h={pretty.y <- pretty(c(0,y))
            lin.pretty.y <- pretty.y
            new.viewport <- TRUE
            cvp <- current.viewport()
            cvp$yscale[1] <- 0
            pushViewport(cvp)
            ## threshold <- c(1/nn, c(2,3) * (pp+1)/nn)
            ## thresh.label <- c("1 / n", "2 (p+1)/n", "3 (p+1)/n")
            ## thresh.id <- c(0, 2*(pp+1)/nn)
          },
         "Cook's distance"={new.viewport <- TRUE
                            cvp <- current.viewport()
                            cvp$yscale <- c(0, max(1.05, y, na.rm=TRUE))
                            pushViewport(cvp)
                            pretty.y <- pretty(cvp$yscale)
                            lin.pretty.y <- pretty.y
                            ## threshold <- c(-1,0,1)
                            ## thresh.label <- c("","0","1")
                            ## thresh.id <- c(-1,1)
                          },
         ## ## if we ever want to restore the 4/(nn-pp-1) criterion
         ## "Cook's distance"={threshold <- c(-1,0,4/(nn-pp-1), 1)
         ##                    thresh.label <- c("","0","4 / (n-p-1)", "1")
         ##                    thresh.id <- c(-1, 4/(nn-pp-1))
         ##                  },
         dffits={## threshold <- c(-2,0,2)*sqrt(pp/nn)
                 ## thresh.label <- c("-2 sqrt(p/n)", "0", "2 sqrt(p/n)")
                 ## thresh.id <- c(-2,2)*sqrt(pp/nn)
               },
         sta.res=,
         "Student del resid"={## threshold <- c(-3,0,3)
                              ## thresh.label <- c(-3,0,3)
                              ## thresh.id <- c(-3,3)
                            },
         dffit={## threshold <- c(-100,0,100) ##arbitrary large numbers
                ## thresh.label <- c(-100,0,100)
                ## thresh.id <- c(-100,100)
              },
         ## DFBETAS
         {## threshold <- c(-2,0,2)/sqrt(nn)
          ## thresh.label <- c("-2 / sqrt(n)", "0", "2 / sqrt(n)")
          ## thresh.id <- c(-2,2)/sqrt(nn)
        }
         )
  threshold <- thresh[[panel.label]]$threshold
  thresh.label <- thresh[[panel.label]]$thresh.label
  thresh.id <- thresh[[panel.label]]$thresh.id


  panel.xyplot(x, y.plot, type="h", ...)

  ## axis 2, "left"
  if (length(lin.pretty.y) > 0)
    panel.axis(
      side="left",
      at=lin.pretty.y,
      labels=pretty.y,
      text.cex=1.0*cex.y, outside=TRUE)

  ## axis 4, "right"
  if (length(threshold) > 0)
    ## We turned off clipping to get outside axes to print
    yscale <- current.viewport()$yscale
  thresh.inside.yscale <-
    (yscale[1] <= threshold) & (threshold <= yscale[2])
  panel.abline(h=threshold[thresh.inside.yscale], lty=2, err=-1)
  panel.axis(
    side="right",
    at=threshold[thresh.inside.yscale],
    labels=thresh.label[thresh.inside.yscale],
    text.cex=cex.threshold*cex.y, tck=-.02, outside=TRUE)

  ## axis 1, "bottom"
  ## use the for() loop to force all items to print
  subs <- case.large[[panel.label]]
  names.x <- rownames[subscripts]
  if (length(subs) > 0)
    for (i in subs)
      panel.axis(
        side="bottom",
        at=i, labels=names.x[i],
        text.cex=1.2*cex.x,
        ticks=TRUE, tck=-.06, outside=TRUE, rot=0)
  if (new.viewport) popViewport()
}

print.lm.case.large <- function(x, dimnames2="Noteworthy Observations", ...) {
  if (length(x) > 0) {
    res <- as.matrix(sapply(x, function(xi) paste(xi, collapse=" ")))
    dimnames(res)[[2]] <- dimnames2
    class(res) <- "noquote"
    print(res)
  }
invisible(x)
}

print.trellis.case <- function(x, ...) {
  xx <- x
  print(x$panel.args.common$case.large)
  class(xx) <- class(x)[-1]
  print(xx)
  invisible(x)
}


## if (FALSE) {
## tmp <- list(a=c(19,23), b=c(13,40,49,56,89), cc=c(5,19), d=c(), e=c(12, 14, 18))
## tmp
## tmp2 <- sapply(tmp, function(inlist) paste(inlist, collapse=" "))
## class(tmp2) <- "noquote"
## tmp2
## tmp3 <- as.matrix(tmp2)
## tmp3
## dimnames(tmp3)[[2]] <- "observations"
## tmp3


## ## from regc.tex
## data(rent)
##
## rent.lm12m <- aov(alf.till ~ lime * cow.dens, data=rent)
## anova(rent.lm12m)
## summary.lm(rent.lm12m)
##
## rent.case12m.trellis <-
##    plot(case(rent.lm12m), rent.lm12m, par.strip.text=list(cex=1.2),
##         layout=c(3,3), main.cex=1.6, col=likertColor(2)[2], lwd=4)
## rent.case12m.trellis
##
## }
