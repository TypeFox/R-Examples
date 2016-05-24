detect.ladder <-
  function(stored, ind=1, ladder, channel.ladder=dim(stored[[1]])[2], ci.upp=1.96, ci.low=1.96, draw=TRUE, dev=50, warn=TRUE, init.thresh=250, sep.index=8, method="cor", avoid=1500, who="sample"){
    ##### stablish a minimum threshold
    x <- stored[[ind]][,channel.ladder]
    thresh=init.thresh
    roxy <- big.peaks.col(x, thresh) # find all peaks with that threshold
    nnn <- length(roxy$pos)
    while(nnn < length(ladder)){
      warning(paste("Warning: Your ladder for individual", who, "is bad, very low intensity was found, we are reducing the threshold 2x trying to find it, don't worry the analysis will keep going anyways"))
      thresh=thresh/2
      roxy <- big.peaks.col(x, thresh) # find all peaks with that threshold
      nnn <- length(roxy$pos)
    }
    ####
    #use the indicator peak at the beggining which is usually the tallest peak in the ladder
    qw <- which(roxy$hei == max(roxy$hei))
    roxy <- list(pos=roxy$pos[-qw], hei=roxy$hei[-qw])
    # get rid of very bad peaks noising the ladder
    roxy <- separate(roxy, shift=sep.index, type="pos") 
    ####################
    if(method == "cor"){
      if(length(roxy$pos) > (length(ladder)+10)){
        print(paste("WOOW too many peaks in this",who,"!! low thresholds throw too many noisy peaks, consider increasing the initial thresold for your ladder, the number of possible combinations is too high to be computed, we will have to do 10,000 samples and get the most likely sizing, you better double check this sample"))
        
        #thresh = init.thresh
        #roxy <- big.peaks.col(x, thresh)
        nono <- which(roxy$pos < avoid)
        roxy <- list(pos = roxy$pos[ -nono], hei = roxy$hei[-nono])
        roxy <- separate(roxy, shift = sep.index, type = "pos")
        pos.mod <- matrix(0, ncol=15000, nrow=length(ladder))
        for(k in 1:15000){
          pos.mod[,k] <- sort(sample(roxy$pos, size=length(ladder), replace=FALSE), decreasing=FALSE)
          
        }
        dd <- apply(pos.mod, 2, function(x5, ladder) {
          cor(x5, ladder)
        }, ladder)
        v <- which(dd == max(dd))[1]
        v2 <- which(roxy$pos %in% pos.mod[, v])
        roxy <- list(pos = pos.mod[, v], hei = roxy$hei[v2], 
                     wei = ladder)
        
      }else{
        mi=length(ladder)
        if(mi > length(roxy$pos)){
          print("ERROR!! using the initial threshold you specified we did not find enough peaks, please try reducing the 'init.thresh' argument and run the analysis again :)")
        }else{
          pos.mod <- combn(roxy$pos, m=mi)
          dd <- apply(pos.mod, 2, function(x5,ladder){cor(x5,ladder)}, ladder)
          v <- which(dd == max(dd))
          v2 <- which(roxy$pos %in% pos.mod[,v])
          roxy <- list(pos=pos.mod[,v], hei=roxy$hei[v2], wei=ladder) 
        } 
      }
    }
    ###################
    if(method == "ci"){
      z <- which(roxy$hei ==  max(roxy$hei))
      mm <- median(roxy$hei[-z]) # get the median height of the peaks
      se2.low <- (sd(roxy$hei[-z])/sqrt(length(roxy$hei[-z]))) * ci.low # produce the confidnce interval
      se2.upp <- (sd(roxy$hei[-z])/sqrt(length(roxy$hei[-z]))) * ci.upp 
      v <- which( (roxy$hei > (mm-se2.low)) & (roxy$hei < (mm+se2.upp))) # reduce the ladder to the peaks inside the 95% confidence interval
      roxy <- list(pos=roxy$pos[v], hei=roxy$hei[v])
      ## keep looking
      vv <- which(diff(roxy$pos) < dev) 
      vv2 <- vv + 1
      # start condition
      while(length(vv) > 0){
        keep <- numeric()
        for(h in 1:length(vv)){
          a1 <- vv[h]
          a2 <- vv2[h]
          a3 <- c(roxy$hei[a1],roxy$hei[a2])
          a4 <- c(a1,a2)
          keep[h] <- (a4[which(a3 == max(a3))])[1]
        }
        keep <- unique(keep)
        '%!in%' <- function(x,y)!('%in%'(x,y))
        keep2 <- unique(c(vv,vv2)[which(c(vv,vv2) %!in% keep)])
        roxy <- list(pos=roxy$pos[-keep2], hei=roxy$hei[-keep2])
        # check again
        vv <- which(diff(roxy$pos) < dev) 
        vv2 <- vv + 1
      }
      
      s1 <- length(roxy$pos)- (length(ladder) - 1)
      s2 <- length(roxy$pos)
      if(s1 <= 0){
        if(warn==TRUE){
          print("Are you sure this is a ladder channel, we did not find a clear pattern, stop a minute to check the plot")
        }
        if(draw != F){
          plot(x, ylim=c(0,mm+se2.upp), type="l")
        }
        #roxy <- list(pos=roxy$pos, hei=roxy$hei, wei= ladder)
        roxy <- list(pos=seq(1,length(ladder)) + rnorm(length(ladder),0,1), hei=seq(1,length(ladder))+ rnorm(length(ladder),0,1), wei= ladder)
      }else{roxy <- list(pos=roxy$pos[s1:s2], hei=roxy$hei[s1:s2], wei= ladder)}
    }
    ################
    # once ladder is found
    ##############################################################
    if(draw == TRUE){
      xx <- roxy$wei
      yy <- roxy$pos
      mod <- lm(yy~xx)
      xlabels <- as.vector(predict(mod, newdata=data.frame(xx=seq(0,max(ladder),by=25))))
      plot(x, type="l", main=paste("Ladder in channel provided",sep=""), xlab="", ylab="Intensity", xaxt="n")
      if(method == "cor"){
        legend("topleft", legend=paste("Correlation=",round(max(dd)), sep="")) 
      }
      axis(1, at=xlabels, labels=seq(0,max(ladder),by=25), cex.axis=0.9)
      abline(v=roxy[[1]], col="red", lty=3)
      if(method == "ci"){
        abline(h=mm, col="blue", lty=3)
        abline(h=(mm+se2.upp), col="blue", lty=3)
        abline(h=(mm-se2.low), col="blue", lty=3)
        legend("topright", legend=c("90% CI", "Peaks found"), col=c("blue", "red"), bty = "n", lty=c(3,3), cex=1)
      }else{legend("topright", legend=c("Peaks found"), col=c("red"), bty = "n", lty=c(3), cex=1)}
      text(x=roxy[[1]], rep(-200, length(ladder)), labels=ladder, cex=0.6)
      #text(x=0, mm, labels="90% CI", cex=0.6)
    }
    return(roxy)
    #########
  }
