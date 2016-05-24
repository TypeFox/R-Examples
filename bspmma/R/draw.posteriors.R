draw.post <- function(mcout,burnin=1000,ind.par=NULL,adjust=1,
                            ...)
  {
    # mcout: is a list with nM components.  Each component is a matrix of
    #   MCMC output from dirichlet.c with nr rows and nstudies+2
    #   columns. The number of simulations was nr -1; each row gives the
    #   output from a MCMC run, except the first row which is the
    #   initial values.
    #   Cols. 1-nstudies are the individual study effects,
    #   col. nstudies+1 is the overall mu, col. nstudies+2 is tau.
    # burnin: number of rows that will be dropped from estimation of the
    #   posterior.
    # ind.par: integer vector, indices of which par.'s to graph
    #
    nM <- length(mcout) # number of M values, no. of overlapping plots
    dOut <- dim(mcout[[1]])# (no. MC runs +1) x (nstudies + 2)
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(mcout[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in mcout must have same size\n")}
      }  
    if (ncycles < 10)
      {
        stop("need at least ten iterations in each chain\n")
      }  
    if (burnin > ncycles - 10)
      {
        burnin.input <- burnin
        if (ncycles < 20) burnin <- 0
        else if (ncycles >=20) burnin <- floor(ncycles/2)
        warning("burnin=",burnin.input,"is too high, burnin changed
                to", burnin)
      }
    param.labels <- dimnames(mcout[[1]])[[2]] #not null if mcout is from
                                    #dirichlet.c or dirichlet.o
    if (is.null(param.labels))
       {
         stop("mcout should have non-null dimnames 2nd component\n")
       }
    if (!any(param.labels == "mu") || !any(param.labels == "tau"))
      {
        stop("need parameters called mu and tau")
      }
    oldpar <- par(no.readonly=TRUE)                
    if (is.null(ind.par))
      {
       ind.par <- (1:nparam)[param.labels=="mu" | param.labels=="tau"
                              | param.labels=="eta"]
      }  
    names(ind.par) <- param.labels[ind.par]            
    nplot <- length(ind.par)            
    Mlabels <- names(mcout)
    arrayx <- array(dim=c(nplot,512,nM))
    arrayy <- array(dim=c(nplot,512,nM))
    param.to.plot <- param.labels[ind.par]
    par.plot.seq <- seq(along=param.to.plot)
    dimnames(arrayx) <- list(param.to.plot,NULL,Mlabels)
    dimnames(arrayy) <- list(param.to.plot,NULL,Mlabels)
    par(mar=c(4,4,3.6,2)+.1,...)  #mar=c(xbot,xlef,xtop,xrig)
    cex.legend <- par("cex")
    rowind <- seq(burnin+1,dOut[1])#1st get indices of rows to extract
    # from the matrices of MCMC output
    for (ii in 1:nplot)
      {
      for (jj in 1:nM)
        {
        dens.obj <- density(mcout[[jj]][rowind,ind.par[ii]],
                            adjust=adjust) 
               arrayx[ii, ,jj] <- dens.obj$x
               arrayy[ii, ,jj] <- dens.obj$y        
        }
      }
    nn <- nplot+100
    if (any(param.to.plot == "mu")) {is.mu <- TRUE;
    imu <- which(param.to.plot == "mu")} else {is.mu <- FALSE; imu <-nn}
    if (any(param.to.plot == "tau")) {is.tau<-TRUE;
    itau<-which(param.to.plot == "tau")} else {is.tau<-FALSE;itau<-nn}
    if (any(param.to.plot == "eta")) {is.eta<- TRUE;
    ieta <- which(param.to.plot == "eta")} else {is.eta<-FALSE;ieta<-nn}   
    ind.rem <- par.plot.seq[-c(imu,itau,ieta)]
    par(mfrow=c(1,2))
    if (is.mu & is.tau) {
     if (is.eta) par(mfrow=c(1,3))
    }
    if (is.mu)
     {
      xrange <- range(arrayx[imu,,])
      yrange <- range(arrayy[imu,,])*1.05
      matplot(arrayx[imu,,],arrayy[imu,,],
       type="l",xlab="",ylab="", #lty 1:5, col 1:6 by default in matplot
       xlim=xrange,ylim=yrange,
       col=c("red","blue","magenta","brown","purple"),
       lty=c(1,2,4,5,6),       
       ...)
      legend(x="topright",
           legend=paste("M=",Mlabels,sep=""),
           lty=c(1,2,4,5,6),
           col=c("red","blue","magenta","brown","purple"),
           bty="n",cex=cex.legend)
      title(main="mu")
     }
    if (is.tau)
     { 
      xrange <- range(arrayx[itau,,])
      yrange <- range(arrayy[itau,,])*1.05        
      matplot(arrayx[itau,,],arrayy[itau,,],
       type="l",xlab="",ylab="",
       xlim=xrange,ylim=yrange,
       col=c("red","blue","magenta","brown","purple"),
       lty=c(1,2,4,5,6),
       ...)
      legend(x="topright",
           legend=paste("M=",Mlabels,sep=""),
           lty=c(1,2,4,5,6),
           col=c("red","blue","magenta","brown","purple"),
           bty="n",cex=cex.legend)
      title(main="tau")
     } 
    if (is.eta)
      {
       xrange <- range(arrayx[ieta,,])
       yrange <- range(arrayy[ieta,,])*1.05
       matplot(arrayx[ieta,,],arrayy[ieta,,],
        type="l",xlab="",ylab="", #lty 1:5, col 1:6 by default in matplot
        xlim=xrange,ylim=yrange,
        col=c("red","blue","magenta","brown","purple"),
        lty=c(1,2,4,5,6),       
        ...)
       legend(x="topright",
           legend=paste("M=",Mlabels,sep=""),
           lty=c(1,2,4,5,6),
           col=c("red","blue","magenta","brown","purple"),
           bty="n",cex=cex.legend)       
      }
    nplot.rem <- nplot - (is.mu+is.tau+is.eta)
    if (nplot.rem < 1)
        {
         par(oldpar) 
         return()
        }         
    if (nplot.rem >= 1 & nplot.rem <= 3) par(mfrow=c(1,nplot.rem))
    if (nplot.rem > 3) 
      {
        n.in.row <- min(floor(nplot.rem/2),4)
        par(mfrow=c(2,n.in.row))
      }
    xrange <- range(arrayx[ind.rem,,])
    yrange <- range(arrayy[ind.rem,,])*1.05
    for (ii in ind.rem)
      {
      matplot(arrayx[ii,,],arrayy[ii,,],
      type="l",xlab="",ylab="",xlim=xrange,ylim=yrange,...)
      title(main=param.to.plot[ii])
    }
    par(oldpar)
  }
