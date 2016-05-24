# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"xp.plot" <-
  function(plot.ids=TRUE,
           idscex=0.7,
           ptscex=0.7,
           prompt=TRUE,

           gamobj=NULL,
           ##main=NULL,
           ##object,
           ##main = NULL,
           ##xlb  = NULL,
           ##ylb  = NULL,
           ##onlyfirst=TRUE,
           ##inclZeroWRES=FALSE,
           ##subset=xsubset(object),
           ## abline=c(0,1),
           ##smooth=TRUE,
           ##abllwd=2,
           ...) {

    if(is.null(gamobj)){
      gamobj <- check.gamobj()
      if(is.null(gamobj)){
        return()
      } else {
      }
    } else {
      c1 <- call("assign",pos=1, "current.gam", gamobj,immediate=T)
      eval(c1)
    }


        #current.gam$terms
    if(length(attr(eval(parse(text="current.gam$terms")),"factors"))==0){
      ##if (length(current.gam$terms)==0){
      cat("\nNo covariates found for this parameter\n")
      return()
    }


    #assign(pos=1,"gamdata",current.gam$data)
    #assign(pos=1,"form",current.gam$form)

    final.gam <- gam(eval(parse(text="current.gam$form")),
                     weights=eval(parse(text="current.gam$weights")),
                     data=eval(parse(text="current.gam$data")))

    #pre.obj <- preplot.gam(final.gam)
    ## HERE is the problem browser()
    pre.obj <- preplot(final.gam)

    ## Significant terms
    trms <- names(pre.obj)
    numplots <- length(trms)

    ## Partial residuals
    parts <- predict(eval(parse(text="current.gam")),type="terms") +
      residuals(eval(parse(text="current.gam")),type="pearson")
    ylm <- range(parts)

    ## add 10% to range
    ylmm <- diff(ylm)*0.05
    ylm[1]= ylm[1]-ylmm
    ylm[2]= ylm[2]+ylmm

    ## plot using the gam.plot function
    ##ylmm <- diff(range(parts))
    ##plot(final.gam,residuals=TRUE,rugplot=FALSE,scale=ylmm)

    ## Get the idlabs
    if(any(is.null(eval(parse(text="current.gam$data$ID"))))){
      ids <- "n"
    } else {
      ids <- eval(parse(text="current.gam$data$ID"))
    }


    ## create enpty list for plots
    plotList <- vector("list",length(trms))


    ## Loop over the terms
    for(i in 1:length(trms)) {

      ##for testing
      ##i=3

      pres <- parts[,trms[i]]

      x <- pre.obj[trms[i]][[1]]$x
      y <- pre.obj[trms[i]][[1]]$y


      ## for testing
      ##idscex=0.7
      ##ptscex=0.7
      ##plot.ids = TRUE
      ##main=NULL

      ##main <- paste("GAM results for \n", trms[i], " on ", current.gam$pars, " (Run ",
      ##                current.gam$runno, ")",sep="")

      main <- NULL

      if(!is.factor(x)) {
        xplot <- xyplot(y~x,res=pres,ids=ids,
                        ylim=ylm,
                        xlab= list(pre.obj[trms[i]][[1]]$xlab,cex=1),
                        ylab= list("Residuals",cex=1),
                        scales=list(cex=1,tck=-0.01),
                        main=main,
                        panel =
                        function(x,y,res,ids,...) {
                          xord <- order(x)
                          panel.xyplot(x[xord],y[xord],type="l",...)
                          if(!any(ids == "n")&& plot.ids==TRUE) {
                            addid(x,res,ids=ids,
                                  idsmode=TRUE,
                                  idsext =0.05,
                                  idscex = idscex,
                                  idsdir = "both")
                          } else {
                            panel.xyplot(x,res,cex=ptscex,col="black")
                          }
                        }
                        )
      } else {
        xplot <-bwplot(y~x,ylim=ylm,res=pres,ids=ids,
                       scales=list(cex=1,tck=-0.01),
                       xlab= list(pre.obj[trms[i]][[1]]$xlab,cex=1),
                       ylab= list("Residuals",cex=1),
                       main=main,
                       panel=
                       function(x,y,res,ids,...) {
                         if(!any(ids == "n")&& plot.ids==TRUE) {
                           addid(jitter(as.numeric(x)),res,ids=ids,
                                 idsmode=TRUE,
                                 idsext =0.05,
                                 idscex = idscex,
                                 idsdir = "both")
                           panel.bwplot(x,y,...)
                         } else {
                           panel.xyplot(jitter(as.numeric(x)),res,cex=ptscex,col="black")
                           panel.bwplot(x,y,...)
                         }
                       }
                       )
      }



      plotList[[i]] <- xplot
}

    plotTitle <- paste("GAM results for covariates on ", eval(parse(text="current.gam$pars")),
                        " (Run ",
                        eval(parse(text="current.gam$runno")), ")",
                        sep="")

    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

  }

