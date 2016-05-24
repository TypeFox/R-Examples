## confounding plot without conditioning (not exported to users)
ConfoundingPlot.None <- function(theta00, theta01, theta10, theta11,
                                  PrY1.setX0,
                                  PrY1.setX1,
                                  epsilon,
                                  color,
                                  legend){

  PrY0.setX0 <- 1 - PrY1.setX0
  PrY0.setX1 <- 1 - PrY1.setX1

  K <- length(color)
  
  ## helper functions
  absdif1.fun <- function(myvec, theta00, theta01, theta10, theta11, 
                          PrY0.setX0, PrY1.setX0, epsilon){

    psi10 <- myvec[1]
    psi11 <- myvec[2]
    
    dif1 <- PrY0.setX0 -
      (theta10 * (1-psi10) + theta11*(1-psi11) + theta00)
    
    dif2 <- PrY1.setX0 -
      ( theta10 * psi10 + theta11 * psi11 + theta01  ) 
    
    bigdif <- max(abs(dif1),  abs(dif2))
    
    return(bigdif < epsilon)
    
  }
  
  absdif2.fun <- function(myvec, theta00, theta01, theta10, theta11, 
                          PrY0.setX1, PrY1.setX1, epsilon){

    psi00 <- myvec[1]
    psi01 <- myvec[2]
    
    dif3 <- PrY0.setX1 -
      ( theta00 * (1 - psi00) + theta01 * (1 - psi01) + theta10 ) 
    
    dif4 <- PrY1.setX1 -
      ( theta00 * psi00 + theta01 * psi01 + theta11 ) 
    
    bigdif <- max(abs(dif3),  abs(dif4))
    
    return(bigdif < epsilon )
  }



  if (legend){
    layout(matrix(c(1,2,3), 1, 3), widths=c(4,4,1.75), heights=c(4,4,4),
           respect=TRUE)
    par(mar=c(5.5, 5.5, 5.5, 5.5))
  }
  else{
    par(mfrow=c(1,2), mar=c(5.5, 5.5, 5.5, 5.5), pty="s")  
  }


  if (legend){
    local.cex=1.75
    axis.cex=1.5
  }
  else{
    local.cex=1.2
    axis.cex=1.0
  }

  
  
  plot(0:1, 0:1, xlim=c(0,1), ylim=c(0,1),
       xlab=expression(paste(psi["00"],": Fraction Helped in (X=0, Y=0)",
           sep="")),
       ylab=expression(paste(psi["01"],
           ": Fraction Always Succeed in (X=0, Y=1)", sep="")),
       cex.lab=local.cex, type="n", cex.axis=axis.cex)
  abline(h=c(0, .2, .4, .6, .8, 1), col="lightgray")
  abline(v=c(0, .2, .4, .6, .8, 1), col="lightgray")  
  
  psimat <- rbind(cbind(0, seq(from=0, to=1, by=.001)),
                  cbind(1, seq(from=0, to=1, by=.001)),
                  cbind(seq(from=0, to=1, by=.001), 0),
                  cbind(seq(from=0, to=1, by=.001), 1))
  
  for (k in 1:K){
        
    indic2 <- apply(psimat, 1, absdif2.fun, theta00=theta00, theta01=theta01,
                    theta10=theta10, theta11=theta11, PrY0.setX1=PrY0.setX1,
                    PrY1.setX1=PrY1.setX1, epsilon=epsilon[k])
    
    psimat2 <- psimat
    psimat2[!indic2,] <- NA
    psimat2 <- na.omit(psimat2)
    polygon(psimat2[chull(psimat2),], border=NA, col=color[k])

  } ## end first K loop

  if (legend){
    local.cex=1.1
  }
  else{
    local.cex=1.2
  }
  
  axis(side=4, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["01"], "): Fraction Hurt in (X=0, Y=1)",
      sep="")), side=4, line=3, cex=local.cex) 
  axis(side=3, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["00"],
      "): Fraction Never Succeed in (X=0, Y=0)", sep="")), side=3, line=3,
        cex=local.cex) 
  

  
  if (legend){
    local.cex=1.75
    axis.cex=1.5
  }
  else{
    local.cex=1.2
    axis.cex=1.0
  }


  plot(0:1, 0:1, xlim=c(0,1), ylim=c(0,1),
       xlab=expression(paste(psi["10"], ": Fraction Hurt in (X=1, Y=0)",
           sep="")),
       ylab=expression(paste(psi["11"],
           ": Fraction Always Succeed in (X=1, Y=1)", sep="")),
       cex.lab=local.cex, type="n", cex.axis=axis.cex)
  
  abline(h=c(0, .2, .4, .6, .8, 1), col="lightgray")
  abline(v=c(0, .2, .4, .6, .8, 1), col="lightgray")  
  
  for (k in 1:K){

    indic1 <- apply(psimat, 1, absdif1.fun, theta00=theta00, theta01=theta01,
                    theta10=theta10, theta11=theta11,  PrY0.setX0=PrY0.setX0,
                    PrY1.setX0=PrY1.setX0, epsilon=epsilon[k])

    psimat1 <- psimat
    psimat1[!indic1,] <- NA
    psimat1 <- na.omit(psimat1)

    polygon(psimat1[chull(psimat1),], border=NA, col=color[k])
    
  } ## end K loop


  if (legend){
    local.cex=1.1
  }
  else{
    local.cex=1.2
  }

  axis(side=4, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["11"],
      "): Fraction Helped in (X=1, Y=1)", sep="")),
        side=4, line=3, cex=local.cex) 
  axis(side=3, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["10"],
      "): Fraction Never Succeed in (X=1, Y=0)", sep="")),
        side=3, line=3, cex=local.cex) 


  if(legend){
    par(mar=c(2,1,5.5,2.1))
    plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
    legendtext <- paste("Abs. Diff. <", round(epsilon, 4))
    legend(0,1, legend=legendtext, fill=color, bty="n", cex=1.75)
    text(x=0.5, y=.2, labels=paste("Assumed Pr(Y(X=0)=0) =",
                     round(PrY0.setX0,2)), cex=1.5)
    text(x=0.5, y=.15, labels=paste("Assumed Pr(Y(X=1)=0) =",
                     round(PrY0.setX1,2)), cex=1.5)
    text(x=0.5, y=.1, labels=paste("Assumed Pr(Y(X=0)=1) =",
                     round(PrY1.setX0,2)), cex=1.5)
    text(x=0.5, y=0.05, labels=paste("Assumed Pr(Y(X=1)=1) =",
                     round(PrY1.setX1,2)), cex=1.5)
  }
  
} ## end ConfoundingPlot.None










## confounding plot within Treated (not exported to users)
ConfoundingPlot.Treated <- function(theta00, theta01, theta10, theta11,
                                    PrY1.setX0.withinTreated,
                                    epsilon,
                                    color,
                                    legend){

  PrY0.setX0.withinTreated <- 1 - PrY1.setX0.withinTreated
  PrY1.setX1.withinTreated <- theta11 / (theta11 + theta10)
  PrY0.setX1.withinTreated <- 1 - PrY1.setX1.withinTreated
  
  K <- length(color)
  
  ## helper functions
  absdif.fun <- function(myvec, theta00, theta01, theta10, theta11, 
                          PrY0.setX0.withinTreated, PrY1.setX0.withinTreated,
                          epsilon){

    psi10 <- myvec[1]
    psi11 <- myvec[2]
    
    dif1 <- PrY0.setX0.withinTreated -
      ((1-psi10) * (theta10/(theta11 + theta10)) +
       (1-psi11) * (theta11/(theta11 + theta10)))
      
    
    dif2 <- PrY1.setX0.withinTreated -
      (psi10 * (theta10/(theta11 + theta10)) +
       psi11 * (theta11/(theta11 + theta10)))
    
    bigdif <- max(abs(dif1),  abs(dif2))
    
    return(bigdif < epsilon)
    
  }



  if (legend){
    layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), widths=c(4,1.75),
           heights=c(1, 4), respect=TRUE)
  }
  else{
    layout(matrix(c(1,1,2,2), 2, 2, byrow=TRUE), widths=c(2, 2),
           heights=c(1, 4.0), respect=TRUE)
  }


  if (legend){
    local.cex=1.4
    axis.cex=1.2
  }
  else{
    local.cex=1.4
    axis.cex=1.2
  }

  par(mar=c(1.1, 4.1, 1.1, 2.1))
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
  text(x=.5, y=0.15, labels="Confounding Plot for Post-Intervention Distribution Within Treated", cex=1.75)
  par(mar=c(5.5, 5.5, 5.5, 5.5))
  

  plot(0:1, 0:1, xlim=c(0,1), ylim=c(0,1),
       xlab=expression(paste(psi["10"], ": Fraction Hurt in (X=1, Y=0)",
           sep="")),
       ylab=expression(paste(psi["11"],
           ": Fraction Always Succeed in (X=1, Y=1)", sep="")),
       cex.lab=local.cex, type="n", cex.axis=axis.cex)

  abline(h=c(0, .2, .4, .6, .8, 1), col="lightgray")
  abline(v=c(0, .2, .4, .6, .8, 1), col="lightgray")  
  
  psimat <- rbind(cbind(0, seq(from=0, to=1, by=.001)),
                  cbind(1, seq(from=0, to=1, by=.001)),
                  cbind(seq(from=0, to=1, by=.001), 0),
                  cbind(seq(from=0, to=1, by=.001), 1))
  
  for (k in 1:K){
        
    indic <- apply(psimat, 1, absdif.fun, theta00=theta00, theta01=theta01,
                    theta10=theta10, theta11=theta11,
                    PrY0.setX0.withinTreated=PrY0.setX0.withinTreated,
                    PrY1.setX0.withinTreated=PrY1.setX0.withinTreated,
                    epsilon=epsilon[k])
    
    psimat[!indic,] <- NA
    psimat <- na.omit(psimat)
    polygon(psimat[chull(psimat),], border=NA, col=color[k])

  } ## end first K loop

  if (legend){
    local.cex=1.2
  }
  else{
    local.cex=1.2
  }
  
  
  axis(side=4, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["11"],
      "): Fraction Helped in (X=1, Y=1)", sep="")),
        side=4, line=3, cex=local.cex) 
  axis(side=3, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["10"],
      "): Fraction Never Succeed in (X=1, Y=0)", sep="")),
        side=3, line=3, cex=local.cex) 



  if(legend){
    par(mar=c(2,1,5.5,2.1))
    plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
    legendtext <- paste("Abs. Diff. <", round(epsilon,4))
    legend(0,1, legend=legendtext, fill=color, bty="n", cex=1.4)
    text(x=.5, y=.25, labels=paste("Assumed Pr(Y(X=0)=0 | X=1) =",
                       round(PrY0.setX0.withinTreated,2)))
    text(x=.5, y=.2, labels=paste("Observed Pr(Y(X=1)=0 | X=1) =",
                       round(PrY0.setX1.withinTreated,2)))
    text(x=.5, y=.15, labels=paste("Assumed Pr(Y(X=0)=1 | X=1) =",
                      round(PrY1.setX0.withinTreated,2)))
    text(x=.5, y=.1, labels=paste("Observed Pr(Y(X=1)=1 | X=1) =",
                      round(PrY1.setX1.withinTreated,2)))
  }
  
  
} ## end ConfoundingPlot.Treated










## confounding plot within Control (not exported to users)
ConfoundingPlot.Control <- function(theta00, theta01, theta10, theta11,
                                    PrY1.setX1.withinControl,
                                    epsilon,
                                    color,
                                    legend){



  PrY1.setX0.withinControl <-  theta01 / (theta00 + theta01)
  PrY0.setX0.withinControl <- 1 - PrY1.setX0.withinControl
  PrY0.setX1.withinControl <- 1 - PrY1.setX1.withinControl
  
  K <- length(color)
  
  ## helper functions
  absdif.fun <- function(myvec, theta00, theta01, theta10, theta11, 
                          PrY0.setX1.withinControl, PrY1.setX1.withinControl,
                          epsilon){

    psi00 <- myvec[1]
    psi01 <- myvec[2]
    
    dif1 <- PrY0.setX1.withinControl -
      ((1-psi00) * (theta00/(theta00 + theta01)) +
       (1-psi01) * (theta01/(theta00 + theta01)))
      
    
    dif2 <- PrY1.setX1.withinControl -
      (psi00 * (theta00/(theta00 + theta01)) +
       psi01 * (theta01/(theta00 + theta01)))
    
    bigdif <- max(abs(dif1),  abs(dif2))
    
    return(bigdif < epsilon)
    
  }



  if (legend){
    layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), widths=c(4,1.75),
           heights=c(1, 4), respect=TRUE)
  }
  else{
    layout(matrix(c(1,1,2,2), 2, 2, byrow=TRUE), widths=c(2, 2),
           heights=c(1, 4.0), respect=TRUE)
  }


  if (legend){
    local.cex=1.4
    axis.cex=1.2
  }
  else{
    local.cex=1.4
    axis.cex=1.2
  }

  par(mar=c(1.1, 4.1, 1.1, 2.1))
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
  text(x=.5, y=0.15, labels="Confounding Plot for Post-Intervention Distribution Within Control", cex=1.75)
  par(mar=c(5.5, 5.5, 5.5, 5.5))
  
  plot(0:1, 0:1, xlim=c(0,1), ylim=c(0,1),
       xlab=expression(paste(psi["00"],": Fraction Helped in (X=0, Y=0)",
           sep="")),
       ylab=expression(paste(psi["01"],
           ": Fraction Always Succeed in (X=0, Y=1)", sep="")),
       cex.lab=local.cex, type="n", cex.axis=axis.cex)
  
  abline(h=c(0, .2, .4, .6, .8, 1), col="lightgray")
  abline(v=c(0, .2, .4, .6, .8, 1), col="lightgray")  
  
  psimat <- rbind(cbind(0, seq(from=0, to=1, by=.001)),
                  cbind(1, seq(from=0, to=1, by=.001)),
                  cbind(seq(from=0, to=1, by=.001), 0),
                  cbind(seq(from=0, to=1, by=.001), 1))
  
  for (k in 1:K){
        
    indic <- apply(psimat, 1, absdif.fun, theta00=theta00, theta01=theta01,
                    theta10=theta10, theta11=theta11,
                    PrY0.setX1.withinControl=PrY0.setX1.withinControl,
                    PrY1.setX1.withinControl=PrY1.setX1.withinControl,
                    epsilon=epsilon[k])
    
    psimat[!indic,] <- NA
    psimat <- na.omit(psimat)
    polygon(psimat[chull(psimat),], border=NA, col=color[k])

  } ## end first K loop

  if (legend){
    local.cex=1.2
  }
  else{
    local.cex=1.2
  }
  
  axis(side=4, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["01"], "): Fraction Hurt in (X=0, Y=1)",
      sep="")), side=4, line=3, cex=local.cex) 
  axis(side=3, at=c(1, .8, .6, .4, .2, 0),
       labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=axis.cex)
  mtext(expression(paste("(1 - ", psi["00"],
      "): Fraction Never Succeed in (X=0, Y=0)", sep="")), side=3, line=3,
        cex=local.cex) 
  


  if(legend){
    par(mar=c(2,1,5.5,2.1))
    plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
    legendtext <- paste("Abs. Diff. <", round(epsilon,4))
    legend(0,1, legend=legendtext, fill=color, bty="n", cex=1.4)
    text(x=.5, y=.25, labels=paste("Observed Pr(Y(X=0)=0 | X=0) =",
                       round(PrY0.setX0.withinControl,2)))
    text(x=.5, y=.2, labels=paste("Assumed Pr(Y(X=1)=0 | X=0) =",
                       round(PrY0.setX1.withinControl,2)))
    text(x=.5, y=.15, labels=paste("Observed Pr(Y(X=0)=1 | X=0) =",
                      round(PrY1.setX0.withinControl,2)))
    text(x=.5, y=.1, labels=paste("Assumed Pr(Y(X=1)=1 | X=0) =",
                      round(PrY1.setX1.withinControl,2)))
  }
  
  
} ## end ConfoundingPlot.Control












## the main ConfoundingPlot function (exported to users)
ConfoundingPlot <- function(theta00, theta01, theta10, theta11,
                            conditioning = c("None", "Treated", "Control"),
                            PrY1.setX0 = NULL,
                            PrY1.setX1 = NULL,
                            PrY1.setX0.withinTreated = NULL,
                            PrY1.setX1.withinControl = NULL,
                            epsilon = 0.025,
                            color = "black",
                            legend=FALSE){
  
  ## error checking
  if (min(theta00, theta01, theta10, theta11) < 0){
    stop("theta value less than 0 in ConfoundingPlot()\n")
  }
  if (max(theta00, theta01, theta10, theta11) > 1){
    stop("theta value greater than 1 in ConfoundingPlot()\n")
  }
  if (!isTRUE(all.equal(1, sum(theta00, theta01, theta10, theta11)))){
    stop("theta values do not sum to 1 in ConfoundingPlot()\n")
  }
  if (min(epsilon) < 0){
    stop("epsilon value less than 0 in ConfoundingPlot()\n")
  }
  
  conditioning <- match.arg(conditioning)

  if (conditioning == "None"){
    if (!is.null(PrY1.setX0.withinTreated)){
      warning("'PrY1.setX0.withinTreated' is non null while\n    'conditioning' is equal to 'None'\n") 
    }
    if (!is.null(PrY1.setX1.withinControl)){
      warning("'PrY1.setX1.withinControl' is non null while\n    'conditioning' is equal to 'None'\n") 
    }
  }

  if (conditioning == "Treated"){
    if (!is.null(PrY1.setX0)){
      warning("'PrY1.setX0' is non null while\n    'conditioning' is equal to 'Treated'\n") 
    }
    if (!is.null(PrY1.setX1)){
      warning("'PrY1.setX1' is non null while\n    'conditioning' is equal to 'Treated'\n") 
    }
    if (!is.null(PrY1.setX1.withinControl)){
      warning("'PrY1.setX1.withinControl' is non null while\n    'conditioning' is equal to 'Treated'\n") 
    }    
  }

  if (conditioning == "Control"){
    if (!is.null(PrY1.setX0)){
      warning("'PrY1.setX0' is non null while\n    'conditioning' is equal to 'Control'\n") 
    }
    if (!is.null(PrY1.setX1)){
      warning("'PrY1.setX1' is non null while\n    'conditioning' is equal to 'Control'\n") 
    }
    if (!is.null(PrY1.setX0.withinTreated)){
      warning("'PrY1.setX0.withinTreated' is non null while\n    'conditioning' is equal to 'Control'\n") 
    }    
  }

  if (length(epsilon) != length(color)){
    stop("length(epsilon) not equal to length(color) in ConfoundingPlot()\n")
  }
  

  ## order the colors so that the widest bands are plotted first
  eps.ord <- order(epsilon)
  epsilon <- rev(epsilon[eps.ord])
  color <- rev(color[eps.ord])
  


  ## functions to actually do the plotting
  if (conditioning == "None"){

    if (is.null(PrY1.setX0)){
      PrY1.setX0 <-  theta01 / (theta00 + theta01)
    }
    if (is.null(PrY1.setX1)){
      PrY1.setX1 <-  theta11 / (theta10 + theta11)
    }
    
    ConfoundingPlot.None(theta00, theta01, theta10, theta11,
                         PrY1.setX0,
                         PrY1.setX1,
                         epsilon,
                         color,
                         legend)
  }
  
  else if (conditioning == "Treated"){

    if (is.null(PrY1.setX0.withinTreated)){
      PrY1.setX0.withinTreated <-  theta01 / (theta00 + theta01)
    }
    
    ConfoundingPlot.Treated(theta00, theta01, theta10, theta11,
                            PrY1.setX0.withinTreated,
                            epsilon,
                            color,
                            legend)
  }

  else if (conditioning == "Control"){

    if (is.null(PrY1.setX1.withinControl)){
      PrY1.setX1.withinControl <-  theta11 / (theta10 + theta11)
    }

    ConfoundingPlot.Control(theta00, theta01, theta10, theta11,
                            PrY1.setX1.withinControl,
                            epsilon,
                            color,
                            legend)
  }


  
  
} ## end ConfoundingPlot


