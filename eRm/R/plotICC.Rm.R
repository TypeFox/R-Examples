# produces ICC plots
# object of class Rm
plotICC.Rm <- function(
  object,
  item.subset = "all",
  empICC = NULL,
  empCI = NULL,
  mplot = NULL,    # ask,mplot added rh 2007-12-01
  xlim = c(-4,4),
  ylim = c(0,1),
  xlab = "Latent Dimension",
  ylab = "Probability to Solve",
  main = NULL,       # main rh 2010-03-06
  col = NULL,
  lty = 1,
  legpos = "left",
  ask = TRUE,
  ...)
{
  # save and reset original graphics parameters
  old_par <- par(mar = c(4,4,3,1)+.25, no.readonly = TRUE)
  on.exit(par(old_par))
  
  if(item.subset != "all" && length(item.subset) == 1L) ask <- FALSE

  X <- object$X
  if(is.null(col)) col <- 1:(max(apply(X, 2L, max, na.rm = TRUE))+1)
  main.arg <- main # rh added 2010-11-23 otherwise always same item in title if NULL

  # some sanity checks

  if(is.null(empICC)){
    emp.plot <- FALSE
  } else if(!(empICC[[1L]] %in% c("raw", "loess", "tukey", "kernel"))) {
    emp.plot <- FALSE
    warning('empICC must be one of "raw", "loess", "tukey", "kernel"!\n')
  } else if(object$model != "RM"){
    warning("Empirical ICCs can only be plotted for a dichotomous Rasch model!\n")
    emp.plot <- FALSE
  } else {
    th.est   <- person.parameter(object)
    thetapar <- th.est$thetapar
    if(length(thetapar) != 1) {   # Too complicated with NA'groups (for each NAgroup separate plots...)
      warning("Empirical ICCs are not produced for different NA groups!\n")
      emp.plot <- FALSE
    } else {
      thetapar.u <- unique(round(unlist(thetapar), 5))
      if(length(thetapar.u) < 4){
        warning("No empirical ICCs for less the 4 different person parameters!\n")
        emp.plot <- FALSE
      } else {
        emp.plot <- TRUE
      }
    }
  }

  theta  <- seq(xlim[1], xlim[2], length.out = 201L)   # x-axis
  p.list <- plist.internal(object, theta)              # matrix of probabilities
  th.ord <- order(theta)

  if(any(item.subset=="all")){
    textlab <- colnames(object$X)
    ivec <- seq_along(p.list)
  } else {
    if(is.character(item.subset)){                         #item names specified
      ivectemp <- matrix(seq_along(p.list), nrow = 1L)
      colnames(ivectemp) <- colnames(object$X)
      ivec <- ivectemp[, item.subset]
      textlab <- item.subset
      textlab[ivec] <- textlab
      it.legend <- item.subset
    } else {                                                    #numeric vector specified
      textlab <- colnames(object$X)[item.subset]
      textlab[item.subset] <- textlab
      ivec <- item.subset
    }
  }

  if(object$model=="RM"){                                       # Rasch model
    p.list <- lapply(p.list,function(x) {x[,-1]})               # Delete 0-probabilites
    p.mat <- matrix(unlist(p.list),ncol=length(p.list))         # matrix with solving probabilities
    text.ylab <- p.mat[(1:length(theta))[theta==median(theta)],]
  }

  ## plot for non RMs #################
  if(object$model != "RM"){
    if(ask) par("ask" = TRUE)                               # added rh 2007-12-01
    if(is.null(mplot)) mplot <- FALSE
    if(mplot) par(mfrow = c(2L, 2L))

    for(j in seq_along(ivec)){                                 # loop for items
      i <- ivec[j]

      yp <- as.matrix(p.list[[i]])
      yy <- yp[th.ord,]

      if(is.null(main.arg)) main <- paste0("ICC plot for item ", textlab[i])    # rh 2010-03-06
      matplot(sort(theta),yy,type="l",lty=lty,col=col,
              #main=paste("ICC plot for item ",textlab[i]),xlim=xlim,  # replaced to allow for user titles rh 2010-03-06
              main=main, xlim=xlim,
              ylim=ylim,xlab=xlab,ylab=ylab,...)
      if(is.character(legpos)) legend(legpos, legend = paste0(c("Category "), 0:(dim(yp)[2]-1)), col = col, lty = lty, ...)  # added rh 2007-12-01
    }

  ## plot for  RMs #####################
  } else {
    if(is.null(mplot)) mplot <- TRUE       ### FIX MM 2012-03-18
    if(length(ivec) == 1) mplot <- FALSE   ### FIX MM 2012-03-18
    if(mplot) par(mfrow = c(2L, 2L))
    if(ask) par("ask" = TRUE)                       # added rh 2007-12-01

    for(j in seq_along(ivec)){                                 #runs over items
      i <- ivec[j]

      yp <- as.matrix(p.list[[i]])
      yy <- yp[th.ord,]
      if(is.null(main.arg)) main<-paste("ICC plot for item ",textlab[i])    # rh 2010-03-06
      matplot(sort(theta),yy,type="l",lty=lty,col=col,
              #main=paste("ICC plot for item ",textlab[i]),xlim=xlim,  # replaced to allow for user titles rh 2010-03-06
              main=main, xlim=xlim,
              ylim=ylim,xlab=xlab,ylab=ylab,...)
              ##ylim=ylim,xlab=xlab,ylab=ylab,"ask"=TRUE,...)

      ## empirical ICC
      if(emp.plot){
         freq.table <- as.matrix(table(rowSums(X), X[,i]))
         rel.freq   <- freq.table[,2]/rowSums(freq.table)
         idx        <- as.numeric(rownames(freq.table))
         xy         <- cbind(th.est$pred.list[[1]]$y[idx+1], rel.freq)

         if(empICC[[1]]=="loess")  if(!is.null(empICC$smooth)) smooth <- empICC$smooth else smooth <- 0.75
         if(empICC[[1]]=="kernel") if(!is.null(empICC$smooth)) smooth <- empICC$smooth else smooth <- 0.5

         nn <- rowSums(freq.table)
         switch(empICC[[1]],
           "raw"={},
           "loess"={xy[,2]<-loess(xy[,2]~xy[,1],span=smooth)$fitted},#+;cyf<-cbind(xy[,2] * nn, nn)},
           "tukey"={xy[,2]<-smooth(xy[,2])},#;cyf<-cbind(xy[,2] * nn, nn)}
           "kernel"={xy[,2]<-ksmooth(xy[,1],xy[,2],bandwidth=smooth,x.points=xy[,1])[[2]]}
         )
         xy[,2] <- ifelse(xy[,2] > 1, 1, ifelse(xy[,2] < 0, 0, xy[,2])) # bounding p in [0,1]

         if(is.null(empICC$type)) empICC$type <- "p"
         if(is.null(empICC$pch)) empICC$pch <- 1
         if(is.null(empICC$col)) empICC$col <- "black"
         if(is.null(empICC$lty)) empICC$lty <- "solid"

         # confidence intervals for empirical ICC
         if(!is.null(empCI)) {
           # functions from prop.test()
           p.L <- function(x, n, alpha){ if (x <= 0) 0 else qbeta(alpha, x, n - x + 1) }
           p.U <- function(x, n, alpha){ if (x >= n) 1 else qbeta(1 - alpha, x + 1, n - x) }
           CINT <- function(x, n, conf.level){
             alpha <- (1 - conf.level)/2
             c(p.L(x,n, alpha), p.U(x,n, alpha))
           }

           if(is.null(empCI$clevel)) empCI$clevel <- 0.95
           if(is.null(empCI$col)) empCI$col <- "red"
           if(is.null(empCI$lty)) empCI$lty <- "dotted"

           cyf <- cbind(xy[,2L]*nn, nn)
           cy  <- apply(cyf, 1L, function(x){ CINT(x[1L], x[2L], empCI$clevel) })

           apply(cbind(xy[,1L], t(cy)), 1L, function(x){ segments(x[1L],x[2L],x[1L],x[3L],lty=empCI$lty,col=empCI$col) })
         }

         # plots the point estimates of the empirical ICC
         lines(xy[,1], xy[,2], type = empICC$type, pch = empICC$pch, col = empICC$col, lty = empICC$lty, ...)
      } # end if(emp.plot)
    }
  }
}
