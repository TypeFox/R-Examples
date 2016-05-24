overview <-
function (my.inds, cols = 1, n.inds = c(1:length(my.inds)), xlimi=c(min(ladder),max(ladder)), ladder, channel.ladder=dim(my.inds[[1]])[2], ploidy=2, ci.upp=1.96, ci.low=1.96, dev=50, method="iter", init.thresh=200, ladd.init.thresh=200, warn=TRUE, env = parent.frame()) 
{
  if(dim(my.inds[[1]])[2] < channel.ladder){
    print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
    stop
  }
  if(method == "ci"){
    print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
  }
  
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
  #if (!require("zoom")) {
   # install.packages("zoom")
  #  require("zoom")
  #}
  ncfp <- c("COL1", "COL2", "COL3", "COL4", "COL5")
  cfp <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
  col.list <- list(NA)
  att1 <- numeric()
  #####################################################################################################
  # this part of the code finds the average height in all samples to be able to plot more uniformly
  #max.lis <- function(l1, colo){
  # takes a list containing a data frmae and returns the maximum value of the column provided in colo argument
  #  y <- max(l1[-c(1:1500),colo]); return(y)
  #}
  #maxi <- lapply(my.inds, max.lis, colo=cols)
  #pro <- mean(unlist(maxi))
  ######################################################################################################
  ### this part extracts all the models for each single plant in my plants
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
  
  ######################################################################
  #########################################################################
  xx <- lapply(my.inds2, function(x, cols){1:length(x[,cols])}, cols=cols)
  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)
  for(h in 1:length(xx)){
    h1 <- n.inds[h]
    ###################
    count <- count + 1
    ###################
    newxx <- as.vector(predict(list.models[[h1]], newdata=data.frame(x=xx[[h]])))
    newyy <- my.inds2[[h]][,cols]
    new.whole.data[[h]] <- list(xx=newxx, yy=newyy) 
    ################################
    setTxtProgressBar(pb, (count/tot)*.5)### keep filling the progress bar
    ################################
  }
  # list.data constains the information for the ladder and picks the closest to lineup the plants
  common <- lapply(list.data, function(x, xlimi){mins <- abs(x$wei - xlimi[1]); y <- x$pos[which(mins == min(mins))]; return(y)}, xlimi=xlimi)
  # gets the maximum height for each plant
  heii <- lapply(my.inds2, function(x){max(x[,cols])[1]})
  # get total height for the plot
  tot.heii <- sum(unlist(heii), na.rm=T)
  # stablish a 1x1 layout matrix
  layout(matrix(1,1,1))
  # start the plot with the first plant
  nn <- n.inds
  plot(new.whole.data[[1]]$xx[-c(1:common[[1]])],y=new.whole.data[[1]]$yy[-c(1:common[[1]])], type="l", yaxt="n",
       xlim=c(xlimi[1],xlimi[2]), ylim=c(0, tot.heii), col=cfp[cols], xlab="Size in base pairs", ylab="Plants selected from bottom to top", xaxt="n", lwd=2)
  axis(1, at=seq(xlimi[1],xlimi[2], by=2), labels=seq(xlimi[1],xlimi[2], by=2))
  b <- sum(unlist(heii)[1]) 
  legend(x=xlimi[1],y=b, legend=paste("Plant",nn[1]), bty="n")
  count <- count + 1
  # make a loop for adding the lines of the other plants
  if(length(n.inds) > 1){
    for(i in 2:length(my.inds2)){
      hh1 <- n.inds[i]
      ###################
      count <- count + 1
      ###################
      a <- sum(unlist(heii)[1:(i-1)]) # maximum height of prevoius plant, bottom
      b <- sum(unlist(heii)[1:i]) # maximim height adding the new plant
      yy <- new.whole.data[[i]]$yy + a
      lines(new.whole.data[[i]]$xx[-c(1:common[[i]])], y=yy[-c(1:common[[i]])], type="l", col=cfp[cols],lwd=2)
      legend(x=xlimi[1],y=b, legend=paste("Plant",nn[hh1]), bty="n")
      #plot(cc[-c(1:common[[i]])], , ylim=c(0, tot.heii))
      ################################
      setTxtProgressBar(pb, (count/tot)*.5)### keep filling the progress bar
      ################################
    }
  }
  close(pb) # close the progress bar
  return(names(my.inds2))
}
