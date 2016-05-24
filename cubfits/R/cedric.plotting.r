plotPTraces <- function(pMat, ...)
{
  input_list <- as.list(list(...))
  pMat <- do.call("rbind", pMat)
  n.traces <- min(dim(pMat)[2], 4)
  
  if("log" %in% names(input_list)){
    logY <- input_list$log
    input_list$log <- NULL
  }else{
    logY <- NULL
  }  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- "Hyperparameter Traces"
  }
  if("xlab" %in% names(input_list)){
    xlab <- input_list$xlab
    input_list$xlab <- NULL
  }else{
    xlab <- "Samples"
  } 
  if("ylab" %in% names(input_list)){
    ylabs <- input_list$ylab
    input_list$ylab <- NULL
  }else{
    if(n.traces == 2)
    {
      ylabs <- c(expression(mu[phi]), expression(sigma[phi]))
    }else{
      ylabs <- c(expression(sigma[epsilon]), expression(mu[phi]), expression(sigma[phi]), "K")
    }
  }
  if("type" %in% names(input_list)){
    type <- input_list$type
    input_list$type <- NULL
  }else{
    type <- "l"
  }
  
  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, ceiling(n.traces/2) ))
  for(i in 1:n.traces)
  {
    if(!is.null(logY) & min(pMat[, i]) > 0)
    {
      do.call(plot, c(input_list, list(x=pMat[, i]), list(xlab=xlab), list(ylab=ylabs[i]), list(type=type), list(log=logY)) )
    }else{
      do.call(plot, c(input_list, list(x=pMat[, i]), list(xlab=xlab), list(ylab=ylabs[i]), list(type=type)) )
    }
    #plot(pMat[, i], xlab=xlab, ylab=ylabs[i], type=type)
  }
  title(main=main, outer=T)
}

plotExpectedPhiTrace <- function(phiMat, ...)
{
  input_list <- as.list(list(...))
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- expression(paste("Trace E[", phi, "]", sep=""))
  }
  if("xlab" %in% names(input_list)){
    xlab <- input_list$xlab
    input_list$xlab <- NULL
  }else{
    xlab <- "Iteration"
  }
  if("ylab" %in% names(input_list)){
    ylab <- input_list$ylab
    input_list$ylab <- NULL
  }else{
    ylab <- expression(paste("E[", phi, "]", sep=""))
  }
  if("type" %in% names(input_list)){
    type <- input_list$type
    input_list$type <- NULL
  }else{
    type <- "l"
  }
  phiMat <- do.call("cbind", phiMat)
  phiMat <- colMeans(phiMat)
  
  
  do.call( plot, c(input_list, list(x=phiMat), list(xlab=xlab), list(ylab=ylab), list(main=main), list(type=type)) )
  abline(h=1, col="red")
}

plotCUB <- function(reu13.df.obs, bMat=NULL, bVec=NULL, phi.bin, n.use.samples=2000,
                     main="CUB", model.label=c("True Model"), model.lty=1, weightedCenters=TRUE)
{
  ### Arrange data.
  aa.names <- names(reu13.df.obs)
  #phi.bin <- phi.bin * phi.scale
  phi.bin.lim <- range(phi.bin)#range(c(phi.bin, phiMat))
  if(is.null(bMat))
  {
    Eb <- bVec  
  }else{
    lbound <- max(0, length(bMat)-n.use.samples)
    ubound <- length(bMat)
    b.mat <- do.call(cbind, bMat[lbound:ubound]) 
    Eb <- rowMeans(b.mat)  
  }
  Eb <- convert.bVec.to.b(Eb, aa.names)
  
  ### Compute.
  ret.phi.bin <- prop.bin.roc(reu13.df.obs, phi.bin, weightedCenters=weightedCenters)
  prediction <- prop.model.roc(Eb, phi.bin.lim)
  
  
  ### Fix xlim at log10 scale. 
  lim.bin <- range(log10(ret.phi.bin[[1]]$center))
  lim.model <- range(log10(prediction[[1]]$center))
  xlim <- c(lim.bin[1] - (lim.bin[2] - lim.bin[1]) / 4,
            max(lim.bin[2], lim.model[2]))
  
  
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = TRUE)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  ### Plot results.
  for(i.aa in 1:length(aa.names))
  {
    tmp.obs <- ret.phi.bin[[i.aa]]
    tmp.roc <- prediction[[i.aa]]
    
    plotbin(tmp.obs, tmp.roc, main = "", xlab = "", ylab = "",
            lty = model.lty, axes = FALSE, xlim = xlim)
    box()
    main.aa <- oneLetterAAtoThreeLetterAA(aa.names[i.aa])
    text(0, 1, main.aa, cex = 1.5)
    if(i.aa %in% c(1, 5, 9, 13, 17)){
      axis(2, las=1)
    }
    if(i.aa %in% 17:19){
      axis(1)
    }
    if(i.aa %in% 1:4){
      axis(3)
    }
    if(i.aa %in% c(4, 8, 12,16)){
      axis(4, las=1)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)
  }
  
  ## adding a histogram of phi values to plot
  hist.values <- hist(log10(phi.bin), plot=FALSE, nclass=30)
  plot(hist.values, axes = FALSE, main="", xlab = "", ylab = "")
  axis(1)
  axis(4, las=1)
  

  ### Add label.
#  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
#        xlim = c(0, 1), ylim = c(0, 1))
#  legend(0.1, 0.8, model.label, lty = model.lty, box.lty = 0)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))  
  #text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
}


plotBMatrixPosterior <- function(bMat, names.aa, interval, param = c("logmu", "deltaeta", "deltat"), main="AA parameter posterior", nclass=100, center=F)
{
  bmat <- convert.bVec.to.b(bMat[[1]], names.aa)
  bmat <- convert.b.to.bVec(bmat)
  names.b <- names(bmat)
  #id.intercept <- grep("Intercept", names.b)
  id.intercept <- grep("log", names.b)
  id.slope <- 1:length(names.b)
  id.slope <- id.slope[-id.intercept]
  
  
  id.plot <- rep(0, length(names.b))
  if(param[1] == "logmu"){
    xlab <- expression(paste("log ( ", mu, " )"))
    id.plot[id.intercept] <- id.intercept
  } else if(param[1] == "deltat"){
    xlab <- expression(paste(Delta, "t"))
    id.plot[id.slope] <- id.slope
  } else if(param[1] == "deltaeta"){
    xlab <- expression(paste(Delta, eta))
    id.plot[id.slope] <- id.slope
  } else if(param[1] == "deltaomega"){
    xlab <- expression(paste(Delta, omega))
    id.plot[id.slope] <- id.slope
  }

  nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  #   nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
  #                rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
  
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  
  ### Plot by aa.
  postMeans <- NULL
  totalncodons <- 1
  for(i.aa in names.aa){
    #id.tmp <- grepl(paste(i.aa,".",sep=""), names.b, fixed=T) & id.plot
    id.tmp <- grepl(paste(i.aa, i.aa, sep="."), names.b, fixed=T) & id.plot
    trace <- lapply(1:length(bMat), function(i){ bMat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)
    if(length(trace) == 0) next
    ncodons <- sum(id.tmp)
    
    ## find x and y limits
    ymax <- vector(mode = "numeric", length = ncodons)
    for(i in 1:ncodons) {
      ymax[i] <- max(hist(trace[interval, i], plot=F, nclass=nclass)$counts)
    }
    if(center){
      if(ncodons > 1){means <- colMeans(trace[interval,])}else{means <- mean(trace[interval,])}
      trace <- trace - mean(means)
    }
    if(ncodons > 1){means <- colMeans(trace[interval,])}else{means <- mean(trace[interval,])}
    postMeans <- c(postMeans, means)
    
    xlim <- range(trace[interval, ])
    ylim <- c(0, max(ymax))
    
    # create empty plot
    main.aa <- oneLetterAAtoThreeLetterAA(i.aa)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = "Frequency", main = main.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
       
    ## Fill plots
    for(i.codon in plot.order){
      hist(trace[interval, i.codon], add=T, nclass=nclass, col=.CF.PT$color[i.codon], lty=0)
    }
    stddev <- format(sd(trace[interval, ]), digits = 3)
    text(x=(xlim[2]+xlim[1])/2, y=ylim[2]-0.1*ylim[2], label=paste("sd =", stddev))
  }

  dens <- density(postMeans)
  xlim <- range(dens$x)
  ylim <- range(dens$y)
  plot(dens, xlim = xlim, ylim = ylim+c(0, 0.2*ylim[2]),
       xlab = xlab, ylab = "Density", main = "Posterior mean distribution")
  stddev <- format(sd(postMeans), digits = 3)
  text(x=(xlim[2]+xlim[1])/2, y=ylim[2]+0.1*ylim[2], label=paste("sd =", stddev))
}

plotTraces <- function(bMat, names.aa, param = c("logmu", "deltaeta", "deltat"), main="AA parameter trace")
{  
  
  bmat <- convert.bVec.to.b(bMat[[1]], names.aa)
  bmat <- convert.b.to.bVec(bmat)
  names.b <- names(bmat)
  id.intercept <- grep("log", names.b)
  id.slope <- 1:length(names.b)
  id.slope <- id.slope[-id.intercept]

  
  id.plot <- rep(0, length(names.b))
  if(param[1] == "logmu"){
    ylab <- expression(paste("log ( ", mu, " )"))
    id.plot[id.intercept] <- id.intercept
  } else if(param[1] == "deltat"){
    ylab <- expression(paste(Delta, "t"))
    id.plot[id.slope] <- id.slope
  } else if(param[1] == "deltaeta"){
    ylab <- expression(paste(Delta, eta))
    id.plot[id.slope] <- id.slope
  } else if(param[1] == "deltaomega"){
    ylab <- expression(paste(Delta, omega))
    id.plot[id.slope] <- id.slope
  }
  
  x <- 1:length(bMat)
  xlim <- range(x)
  
  ### Trace plot.
  nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
#   nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
#                rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
  
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  ### Plot by aa.
  for(i.aa in names.aa){
    id.tmp <- grepl(paste(i.aa, i.aa, sep="."), names.b, fixed=T) & id.plot
    trace <- lapply(1:length(bMat), function(i){ bMat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)
    if(length(trace) == 0) next
    ylim <- range(trace, na.rm=T)
    
    main.aa <- oneLetterAAtoThreeLetterAA(i.aa)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = trace[, i.codon], col = .CF.PT$color[i.codon])
    } 
  }
}

oneLetterAAtoThreeLetterAA <- function(letter)
{
  letters <- list(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", 
                  G="Gly", H="His", I="Ile", L="Leu", K="Lys", F="Phe", P="Pro", 
                  S=expression("Ser"[4]), Z=expression("Ser"[2]), T="Thr", Y="Tyr", V="Val", M="Met", W="Trp")
  return(letters[[letter]])
}
