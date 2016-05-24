overview2 <- function (my.inds, cols = 1, n.inds = NULL, xlim = NULL, ylim = NULL, ladder, channel.ladder = NULL, ploidy = 2, ci.upp = 1.96, 
                       ci.low = 1.96, dev = 50, method="iter", init.thresh=200, ladd.init.thresh=200, lwd=.25, warn=TRUE, min.panel=100, env = parent.frame()) 
{
  
  
  if(method == "ci"){
    print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
  }
  ## channel where the ladder is located
  if(is.null(channel.ladder )){
    channel.ladder <- dim(my.inds[[1]])[2]
  }else{ channel.ladder <- channel.ladder}
  if(dim(my.inds[[1]])[2] < channel.ladder){
    print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
    stop
  }
  ## number of samples to do
  if(is.null(n.inds )){
    n.inds <- c(1:length(my.inds))
  }else{n.inds <- n.inds}
  ## x limits
  if(is.null(xlim )){
    xlim <- c(min(ladder), max(ladder))
  } else{ xlim <- xlim}
  ####################
  ## initialize the progress bar
  count <- 0
  tot <- length(n.inds)
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  #####################
  
  my.inds2 <- list(NA)
  for (i in 1:length(n.inds)) {
    v1 <- n.inds[i]
    my.inds2[[i]] <- my.inds[[v1]]
    names(my.inds2)[i] <- names(my.inds)[i]
  }
  my.inds <- my.inds2
  ncfp <- c("COL1", "COL2", "COL3", "COL4", "COL5")
  cfp <- c("cornflowerblue", "chartreuse4", "gold2", "red", 
           "orange", "purple")
  col.list <- list(NA)
  att1 <- numeric()
  ############################################
  list.data <- list(NA)
  if(exists("list.data.covarrubias")){
    list.data <- env$list.data.covarrubias
  }else{
    list.ladders <- lapply(my.inds, function(x){y <- x[,channel.ladder]; return(y)})
    # extract ladder channels for all plants
    list.data <- lapply(list.ladders, find.ladder, ladder=ladder, ci.upp=ci.upp, ci.low=ci.low, draw=F, dev=dev, warn=warn, method=method,init.thresh=ladd.init.thresh)
  } # this models uses indexes and predicts base airs
  list.models <- lapply(list.data, function(da){y <- da[[3]]; x <- da[[1]];mod <- lm(y~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5), data=da); return(mod)})
  # this models uses pairs and predicts indexes
  list.models.inv <- lapply(list.data, function(da){x <- da[[3]]; y <- da[[1]];mod <- lm(y~ x, data=da); return(mod)})
  ##############################################
  xx <- lapply(my.inds2, function(x, cols) {
    1:length(x[, cols])
  }, cols = cols)
  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)
  for (h in 1:length(xx)) {
    h1 <- n.inds[h]
    ###################
    count <- count + 1
    ###################
    newxx <- as.vector(predict(list.models[[h1]], newdata = data.frame(x = xx[[h]])))
    newyy <- my.inds2[[h]][, cols]
    new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
    ################################
    setTxtProgressBar(pb, (count/tot)*.5)### keep filling the progress bar
    ################################
  }
  common <- lapply(list.data, function(x, xlim) {
    mins <- abs(x$wei - xlim[1])
    y <- x$pos[which(mins == min(mins))][1]
    return(y)
  }, xlim = xlim)
  heii <- lapply(my.inds2, function(x) {
    max(x[, cols])[1]
  })
  
  ## ---------------------------------------
  ## provide inital guesses of which should be the panel peaks
  my.panel <- lapply(new.whole.data,
         function(popo){
           pann <- big.peaks.col(popo$yy, tre=init.thresh)
           pann2 <- popo$xx[pann$pos]
           pann3 <- list(pos=pann$pos, hei=pann$hei, wei=pann2)
           pkpn <- separate(pann3, type="bp", shift=1)
           return(list(wei=pkpn$wei, hei=pkpn$hei))
         }
  )
  #plot(unlist(my.panel))
  ## unlist all the peaks found for all the plants
  allpan <- unlist(lapply(my.panel, function(x){x$wei}))
  allhei <- unlist(lapply(my.panel, function(x){x$hei}))
  ## create a vector to store a the good peaks
  panel1.1 <- numeric()
  heis1.1 <- numeric()
  for(za in seq(1,500, by=1)){
    step1 <- abs(za- allpan)
    good <- which(step1 < 0.48) # peak present at with minumum error of x bp
    if(length(good) > (length(n.inds)*.05)){ # more than 20% of the times present
      panel1.1[za] <- mean(allpan[good])
      heis1.1[za] <- mean(allhei[good])
    }else{panel1.1[za] <- NA; heis1.1[za] <- NA}
  }
  if(is.null(xlim)){
    panel.sugg <- panel1.1[-which(panel1.1 < min.panel | is.na(panel1.1))]
    heis.sugg <- heis1.1[-which(panel1.1 < min.panel | is.na(panel1.1))]
  }else{
    prov <- panel1.1[which(panel1.1 > xlim[1] & panel1.1 < xlim[2])]
    bad <- which(is.na(prov))
    if(length(bad) > 0){panel.sugg <- prov[-bad]}else{panel.sugg <- prov}
    #--
    prov2 <- heis1.1[which(panel1.1 > xlim[1] & panel1.1 < xlim[2])]
    bad2 <- which(is.na(prov2))
    if(length(bad2) > 0){heis.sugg <- prov[-bad2]}else{heis.sugg <- prov2}
  }
  
  ## ------------------------------------------
  ## parameters for plots and lines
  tot.heii <- max(unlist(heii), na.rm = T)
  ## ylims defaults
  if(is.null(ylim )){
    ylim <- c(0, tot.heii)
  }else{ylim <- ylim}
  ##
  layout(matrix(1, 1, 1))
  nn <- n.inds
  plot(new.whole.data[[1]]$xx[-c(1:common[[1]])], y = new.whole.data[[1]]$yy[-c(1:common[[1]])], 
       type = "l", xlim = c(xlim[1], xlim[2]), 
       ylim=c(ylim[1],ylim[2]), yaxt = "n", col = transp(cfp[cols],0.6), xlab = "Size in base pairs", 
       ylab = "DNA intensity in RFU", xaxt = "n", 
       lwd = lwd)
  axis(1, at = seq(xlim[1], xlim[2], by = 2), labels = seq(xlim[1], 
                                                             xlim[2], by = 2), cex.axis=0.7)
  axis(2, at = seq(0, tot.heii, by = 500), labels = seq(0, tot.heii, by = 500), las=1, cex.axis=0.4)
  if(length(n.inds) == 1){
    count <- count + 50
    setTxtProgressBar(pb, (count/tot)*.5)
  }else{count <- count + 1}
  
  #b <- sum(unlist(heii)[1])
  #legend(x = xlim[1], y = b, legend = paste("Plant", nn[1]), 
  #      bty = "n")
  if (length(n.inds) > 1) {
    for (i in 2:length(my.inds2)) {
      ###################
      count <- count + 1
      ###################
      a <- sum(unlist(heii)[1:(i - 1)])
      b <- sum(unlist(heii)[1:i])
      yy <- new.whole.data[[i]]$yy
      lines(new.whole.data[[i]]$xx[-c(1:common[[i]])], 
            y = yy[-c(1:common[[i]])], type = "l", col = transp(cfp[cols],0.6), 
            lwd = lwd)
      #legend(x = xlim[1], y = b, legend = paste("Plant", 
      #                                          nn[i]), bty = "n")
      ################################
      setTxtProgressBar(pb, (count/tot)*.5)### keep filling the progress bar
      ################################
    }
  }
  legend("topright", legend="Peaks suggested", pch=20, col="red", bty="n", cex=0.75)
  
  points(x=panel.sugg, y=heis.sugg, pch=20, cex=0.7, col="red")
  points(x=panel.sugg, y=heis.sugg, cex=0.9, col="black")
  
  close(pb) # close the progress bar
  cat("\n THE PEAKS RETURNED ARE SUGGESTIONS. \n   My suggestion: \n a) Use the locator function, i.e. ''my.panel <- locator(type='p', pch=20, col='red')$x'' \n b) Click over the peaks you want to include in your panel \n c) Press the 'esc' key when done selecting peaks \n d) Make sure to provide the panel vector in the score.easy() function \n \n")
  
  return(panel.sugg)
}