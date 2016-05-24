find.ladder <-
  function(x, ladder, ci.upp=1.96, ci.low=1.96, draw=TRUE, dev=50, warn=TRUE, init.thresh=250, sep.index=8, method="iter", avoid=1000, who="sample", attempt=10, cex.title=0.8){
    #roundUP <- function(x) 10^ceiling(log10(x))
    
    MSE <- function(x,y){
      X <- cbind(1, x)
      qr.X <- qr(X)
      b <- t(qr.Q(qr.X)) %*% y
      R <- qr.R(qr.X)
      beta <- as.vector(backsolve(R, b))
      fit <- X%*%beta
      res <- list(mse=sum((y-fit)^2), beta=beta[2])
      return(res)
    }
    MSE2 <- function(x,y){
      X <- cbind(1, x)
      qr.X <- qr(X)
      b <- t(qr.Q(qr.X)) %*% y
      R <- qr.R(qr.X)
      beta <- as.vector(backsolve(R, b))
      fit <- X%*%beta
      mse <- sum((y-fit)^2)# is actually SSE
      sst <- sum((y-mean(y))^2)
      r2 <- 1 - (mse/sst) #cor(x,y)^2#mse/sum((y-mean(y,na.rm=TRUE))^2)
      res <- list(mse=mse, beta=beta, r2=r2)
      return(res)
    }
    ##### stablish a minimum threshold
    thresh=init.thresh
    roxy <- big.peaks.col(x[1:length(x)], thresh) # find all peaks with that threshold
    nnn <- length(roxy$pos)
    while(nnn < length(ladder)){
      #print(paste("Warning: Your ladder for the sample",attr(x,"name") ,"is bad, very low intensity in RFU was found, we are reducing the threshold 2x trying to find it, don't worry too much the analysis will keep going anyways"))
      cat("\nReducing threshold 2x to find ladder \n")
      thresh=thresh/2
      roxy <- big.peaks.col(x[1:length(x)], thresh) # find all peaks with that threshold
      nnn <- length(roxy$pos)
    }
    ####
    #use the indicator peak at the beggining which is usually the tallest peak in the ladder
    # get rid of very bad peaks noising the ladder
    whot <- length(roxy$pos)*.2
    what <- which(roxy$hei == max(roxy$hei))
    # if you tallest peak is far beyond don't adjust
    #if(what > whot){
    #  roxy <- roxy
    #}else{
    roxy <- separate(roxy, shift=sep.index, type="pos")
    ii <-  which(roxy$hei == max(roxy$hei)) + 1#which(roxy$hei < 2*se(roxy$hei))[1]#
    iii <- length(roxy$hei) 
    roxy <- list(pos=roxy$pos[ii:iii], hei=roxy$hei[ii:iii])
    #}
    #####################
    if(method == "iter"){
      #caller
      #doc <- length(ladder)
      #popo <- (length(roxy$pos) - length(ladder))
      #if( popo > 9){ # we only allow for 6 extra peaks
      #  popo2 <- sort(roxy$hei, decreasing=FALSE)
      #  ggg <- which(roxy$hei > popo2[popo - 9])
      #  roxy <- list(pos=roxy$pos[ggg], hei=roxy$hei[ggg])
      #}
      # by default attempt=10 making all combinations of 1st 10 peaks
      step1 <- combn(roxy$pos[1:attempt],3)
      step2 <- apply(step1/10, 2, MSE, y=ladder[1:3]) # function(x,y){sum((summary(lm(I(x)~y))$residuals)^2)}create a custom function with matrices to extract residuals and calculate MSE
      mse <- unlist(lapply(step2, function(x){x$mse}))
      covs <- apply(step1, 2, function(x,y){cov(x,y)}, y=ladder[1:3])
      #beta <- unlist(lapply(step2, function(x){x$beta}))
      step2 <- mse * covs
      step3 <- step1[,which(step2 < sort(step2, decreasing=FALSE)[20])] # 15 models with least MSE
      #86,176
      
      #step1 <- combn(roxy$pos[1:attempt],3)
      #step2 <- apply(step1, 2, function(x,y){cor(x,y)}, y=ladder[1:3])
      #step3 <- step1[,which(step2 >= sort(step2, decreasing = TRUE)[10])]#
      
      step4 <- apply(step3,2,function(x,y){which(y %in% x)}, y=roxy$pos) # which peaks are
      #caller  6
      ############
      #-----------
      ############
      caller <- function(roxy,www, ladder.call,x){
        threshold <- length(x)
       
        posi <- numeric()
        fact2 <- length(ladder.call)
        # www <- step4[,5]; ladder.call=ladder
        ############################################################
        ## short initial model to avoid people using a ladder.call too long, we cut the ladder.call if necessary
        expect <- roxy$pos[www]
        xxx <- ladder.call[c(1:3)]
        modx <- lm(expect~poly(xxx, degree=1))
        expecto <- predict(modx, data.frame(xxx=ladder.call))# + facto
        ladder.call <- ladder.call[which(expecto < threshold*.85)]
        
        available <- length(roxy$pos) - length(ladder.call)
        ava2 <-  length(ladder.call) - abs(available)
        if(available > 0){
          if((length(ladder.call)-1) < 3){
            tope <- length(ladder.call)
          }else{tope <- length(ladder.call)-1}
        }else{tope <- ava2-2}
        #if(available > 0){tope <- length(ladder.call)-1}
        #}else{tope <- abs(available)-1}
        #####
        expect <- rep(NA,tope+1)# (length(ladder.call)-maxo)
        for(i in 3:tope){ #(length(ladder.call)-maxo)
          if(i == 3 & i != tope){#length(ladder.call)-maxo
            expect[1:3] <- roxy$pos[www]
            xxx <- ladder.call[c(1:3)]
            #mod1 <- lm(expect[1:3]~poly(xxx, degree=1))
            #expecto <- predict(mod1, data.frame(xxx=ladder.call))# + facto
            mod <- MSE2(xxx, expect[1:3])
            beta <- (mod)$beta
            expecto <- as.vector(beta[1] + matrix(ladder.call) %*% beta[-1])
            # now that was predicted get the real one
            act <- roxy$pos[-which(roxy$pos %in% expect)]
            yoyo <- abs(expecto[i+1] - act)
            good <- which(yoyo == min(yoyo))
            expect[i+1] <- act[good]
            # --------------------------------
            #if(summary(mod1)$r.squared < .9){
            if(mod$r2 < .9){
              i= tope #length(ladder.call) - maxo
            }
            # --------------------------------
          }
          if(i > 3 & i <= 5 ){
            xx <- ladder.call[c(1:i)]
            #mod1 <- lm(expect[1:i]~poly(xx, degree=1))
            #expecto <- predict(mod1, data.frame(xx=ladder.call))
            
            mod <- MSE2(xx, expect[1:i])
            beta <- (mod)$beta
            expecto <- as.vector(beta[1] + matrix(ladder.call) %*% beta[-1])
            #now that was predicted get the real one
            act <- roxy$pos[-which(roxy$pos %in% expect)]
            yoyo <- abs(expecto[i+1] - act)
            good <- which(yoyo == min(yoyo, na.rm = TRUE))
            #error[i] <- abs(expecto[i] - posi[i-1])
            expect[i+1] <- act[good]
            # --------------------------------
            #if(summary(mod1)$r.squared < .9){
            if(mod$r2 < .9){
              i= tope #length(ladder.call) - maxo
            }
            # --------------------------------
          }
          if(i > 5 ){
            #xx <- ladder.call[c(1:i)]
            #mod1 <- lm(expect[1:i]~poly(xx, degree=4, raw=TRUE))
            #expecto <- predict(mod1, data.frame(xx=ladder.call))
            xx <- cbind(ladder.call[c(1:i)],ladder.call[c(1:i)]^2,ladder.call[c(1:i)]^3,ladder.call[c(1:i)]^4)
            mod <- MSE2(xx, expect[1:i])
            beta <- (mod)$beta
            if(length(which(is.na(beta))) > 0){
              beta[which(is.na(beta))] <- 0
            }
            toto <- cbind(matrix(ladder.call),matrix(ladder.call)^2, matrix(ladder.call)^3,matrix(ladder.call)^4)
            expecto <- cbind(rep(1,dim(toto)[1]),toto) %*% beta
            #now that was predicted get the real one
            act <- roxy$pos[-which(roxy$pos %in% expect)]
            yoyo <- abs(expecto[i+1] - act)
            good <- which(yoyo == min(yoyo))
            #error[i] <- abs(expecto[i] - posi[i-1])
            expect[i+1] <- act[good]
            # --------------------------------
            #if(summary(mod1)$r.squared < .9){
            if(is.na(mod$r2)){ # if model is too bad just assign a zero value
              mod$r2 <- 0.1
            }
            if(mod$r2 < .9){
              i= tope #length(ladder.call) - maxo
            }
            # --------------------------------
          }
          if(i == tope & i != 3){ #length(ladder.call)-maxo
            if(i < 5){
              expect[1:3] <- roxy$pos[www]
              xx <- ladder.call[c(1:i)]
            }else{xx <- cbind(ladder.call[c(1:i)],ladder.call[c(1:i)]^2,ladder.call[c(1:i)]^3,ladder.call[c(1:i)]^4)}
            mod <- MSE2(xx, expect[1:i])
            beta <- (mod)$beta
            if(length(which(is.na(beta))) > 0){
              beta[which(is.na(beta))] <- 0
            }
            if(i < 5){
              toto <- cbind(matrix(ladder.call))
            }else{toto <- cbind(matrix(ladder.call),matrix(ladder.call)^2, matrix(ladder.call)^3,matrix(ladder.call)^4)}
            expecto <- cbind(rep(1,dim(toto)[1]),toto) %*% beta
            #now that was predicted get the real one
            act <- roxy$pos[-which(roxy$pos %in% expect)]
            yoyo <- abs(expecto[i+1] - act)
            good <- which(yoyo == min(yoyo))
            #error[i] <- abs(expecto[i] - posi[i-1])
            expect[i+1] <- act[good]
          }
          if(i == tope & i == 3){
            expect[1:3] <- roxy$pos[www]
          }
          
          ##
        } # end of for loop
        
        ###########################################
        posi <- expect
        ## get rid of selected peaks after reaching the maximum values
        tutu <- abs(length(x) - posi)
        #if(length(posi) > 3){
        posi <- posi[1:which(tutu == min(tutu,na.rm = TRUE))]
        #}
        heii <- roxy$hei[which(roxy$pos %in% posi)]
        
        fact3 <- length(posi)/fact2
        
        if(length((posi)) < 6){
          fact <- summary(lm(ladder.call[1:length(posi)]~poly(posi, degree=length((posi))-1)))$r.squared * fact3
        }else{
          fact <- summary(lm(ladder.call[1:length(posi)]~poly(posi, degree=5)))$r.squared * fact3
        }
        #plot(ladder.call~posi)
        roxy2 <- list(pos=posi, hei=heii, wei=ladder.call[1:length(posi)], corr=abs(cor(ladder.call[1:length(posi)],posi)), error=fact)#sum(error, na.rm=TRUE))
        return(roxy2)
      }
      ############
      #-----------
      ############
      ## end of caller
      rt <- apply(data.frame(step4), 2, FUN=caller, roxy=roxy, ladder.call=ladder,x=x)#
      corrs3 <- unlist(lapply(rt, function(x){x$error})) #; dis[which(dis == Inf)] <- 1
      #corrs <- unlist(lapply(rt, function(x){x[[4]]})) # extract correlations
      roxy3 <- rt[[which(corrs3 == max(corrs3))]]
      if(draw == TRUE){
        limi <- sort(roxy3$hei, decreasing = TRUE)
        plot(x, type="l", xaxt="n", ylim=c(0,(limi[3]+1000)), cex.axis=0.6, las=2, xlim=c((min(roxy3$pos)-100),(max(roxy3$pos)+100)), col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, main=attributes(x)$mycomm, cex.main=cex.title)
        axis(1, at=roxy3$pos, labels=roxy3$wei, cex.axis=0.6)
        points(x=roxy3$pos, y=roxy3$hei,cex=1.1, col=transp("black",0.85))
        points(x=roxy3$pos, y=roxy3$hei, pch=20, col=transp("red",0.7))
        legend("topleft", legend=paste("Correlation:",round(roxy3$corr, digits=4), sep=""), bty="n")
        legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
        
      }
      roxy <- roxy3
    }
    #####################
    #####################
    #####################
    if(method == "cicor"){
      roxy <- separate(roxy, shift=40, type="pos")
      nnn <- length(ladder)
      mm <- median(roxy$hei)
      vvv <- which(roxy$hei == max(roxy$hei))
      se <- sd(roxy$hei[-vvv])/sqrt(length(roxy$hei[-vvv]))
      reduced <- roxy$pos[which(roxy$hei < mm+(1*se) & roxy$pos > mm-(1*se) & roxy$hei > init.thresh)]
      reduced2 <- roxy$pos[which(roxy$pos >= min(reduced))]
      
      ## if there's still just toomany roxy after the reduced search
      if(length(reduced2) >= nnn & length(reduced2) < nnn+8){
        all.combs <- combn(reduced2, m=length(ladder))
        cors <- apply(all.combs, 2, function(x,y){cor(x,y)}, y=ladder)
        # positions of roxy found
        found <- all.combs[,which(cors == max(cors))]
        roxy <- list(pos=found, hei=x[found], wei=ladder)
        
        nono <- length(x) - max(roxy$pos)
        
        if(draw == TRUE){
          plot(x, type="l", xaxt="n", ylab="Intensity", col=transp("gray29",0.6), main=attributes(x)$mycomm, cex.main=cex.title)
          axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
          abline(v=found, col="red", lty=3)
        }
        # reduced search to maximum correlations
      }else{# there's actually no ladder so no good correlation was found
        roxy <- list(pos=seq(1,nnn) + rnorm(nnn,0,1), hei=seq(1,nnn)+ rnorm(nnn,0,1), wei= ladder)
        print("Friend I was not able to find a ladder in this sample or you used the wrong ladder, look the plot")
        plot(x, type="l", col=transp("black",0.5), main=attributes(x)$mycomm, cex.main=cex.title)
      }
    }
    #####################
    if(method == "red"){
      nnn <- length(roxy$pos) - length(ladder)
      
      if(nnn > 15){
        yoy <- length(roxy$pos) - length(ladder)
        yoy2 <- rev(seq(1,yoy,by=1))
        corres <- numeric()
        for(i in 1:yoy){
          corres[i] <- cor(roxy$pos[i:(length(roxy$pos)-yoy2[i])], ladder)
        }
        vv <- which(corres == max(corres))[1]
        roxy <- list(pos=roxy$pos[vv:((vv-1)+length(ladder))], hei=roxy$hei[vv:((vv-1)+length(ladder))], wei=ladder)
        #lapply(roxy,length)
        if(draw == TRUE){
          plot(x, type="l", xaxt="n", ylim=c(0,(max(roxy$hei)+1000)), cex.axis=0.6, las=2, xlim=c((min(roxy$pos)-100),(max(roxy$pos)+100)), col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, main=attributes(x)$mycomm,cex.main=cex.title)
          axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
          points(x=roxy$pos, y=roxy$hei,cex=1.1, col=transp("black",0.85))
          points(x=roxy$pos, y=roxy$hei, pch=20, col=transp("red",0.7))
          legend("topleft", legend=paste("Correlation found:",round(max(cors), digits=4), sep=""), bty="n")
          legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
        }
      }else{
        nnn <- length(ladder)
        #### get the length of the peak region
        le <- length(roxy$pos[1]:(roxy$pos[length(roxy$pos)]))
        #### create a vector of expected indexes according to our ladder
        tra <- (le * ladder) / max(ladder)
        #### number of possible tries
        tries <- length(seq(tra[length(tra)], roxy$pos[length(roxy$pos)], 1))+100
        
        # define function returning absolute values
        abso <- function(test1, pos){
          test1 <- matrix(test1, nrow=1)
          res <- apply(test1,2,function(x, y){
            xx1 <- abs(as.vector(x) - y)
            z <- y[which(xx1 == min(xx1))[1]]
            return(z)}, y=pos)
          return(res)
        }
        
        # get all possible tests
        all.tests <- apply(matrix(1:tries,nrow=1), 2, function(q1, q2){q1+q2},q2=tra)
        # get all possible absolute differences
        all.abso <- apply(all.tests, 2, abso, pos=roxy$pos)
        # get all possible correlations
        all.cors <- apply(all.abso, 2, function(x,y){cor(x,y)}, y=ladder)
        # get all sum of squares
        step1 <- (all.tests - all.abso)^2
        all.ss2 <- apply(step1, 2, sum)
        # get all standarized variances and final response
        sss2 <- abs ( (all.ss2 - mean(all.ss2) )/ sd(all.ss2)) # the smaller the better
        response <- all.cors/sss2
        
        vv4 <- which(response >= min(sort(response, decreasing=TRUE)[1:5]))
        #vv4 <- which(response >= .99)
        reduced <- sort(unique(as.vector(all.tests[,vv4])))
        reduced2 <- roxy$pos[which(roxy$pos >= min(reduced))]
        
        ## if there's still just toomany roxy after the reduced search
        if(length(reduced2) >= nnn & length(reduced2) < nnn+8){
          all.combs <- combn(reduced2, m=length(ladder))
          cors <- apply(all.combs, 2, function(x,y){cor(x,y)}, y=ladder)
          # positions of roxy found
          found <- all.combs[,which(cors == max(cors))]
          roxy <- list(pos=found, hei=x[found], wei=ladder)
          if(draw == TRUE){
            #facto <- roundUP(max(roxy$hei) / 10)
            #yyline <- seq(0,2*max(roxy$hei), by=facto)
            
            #if((length(x) - min(roxy$pos)) > 500){nono <- min(roxy$pos) - 200}else{nono <- 0}
            plot(x, type="l", xaxt="n", ylim=c(0,(max(roxy$hei)+1000)), cex.axis=0.6, las=2, xlim=c((min(roxy$pos)-100),(max(roxy$pos)+100)), col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, main=attributes(x)$mycomm, cex.main=cex.title)
            #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
            
            #abline(h=yyline, col="white", lwd=1.5, lty=3)
            #lines(x=x, lwd=2, col=transp("grey25",0.6))
            #polygon(x=x, col=transp("grey35",0.6))
            axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
            #abline(v=found, col="red", lty=3)
            points(x=roxy$pos, y=roxy$hei,cex=1.1, col=transp("black",0.85))
            points(x=roxy$pos, y=roxy$hei, pch=20, col=transp("red",0.7))
            legend("topleft", legend=paste("Correlation found:",round(max(cors), digits=4), sep=""), bty="n")
            legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
            
          }
          #which(sss2 == min(sss2))
          # reduced search to maximum correlations
        }else{# there's actually no ladder so no good correlation was found
          roxy <- list(pos=seq(1,nnn) + rnorm(nnn,0,1), hei=seq(1,nnn)+ rnorm(nnn,0,1), wei= ladder)
          print("Friend I don't think you have ladder in this sample or you used the wrong ladder, look the plot")
          plot(x, type="l", main=attributes(x)$mycomm, cex.main=cex.title)
        }
      }
    }
    ### end of method "red"
    #####################
    #####################
    #####################
    if(method == "cor"){ #EXHAUSTIVE CORRELATION
      if(length(roxy$pos) > (length(ladder)+10)){
        print(paste("WOOW too many peaks in this",who,"!! low thresholds throw too many noisy peaks, consider increasing the initial thresold for your ladder, the number of possible combinations is too high to be computed, we will have to do 15,000 samples and get the most likely sizing, you better double check this sample"))
        
        #thresh = init.thresh
        #roxy <- big.peaks.col(x, thresh)
        nono <- which(roxy$pos < avoid)
        roxy <- list(pos = roxy$pos[ -nono], hei = roxy$hei[-nono])
        roxy <- separate(roxy, shift=sep.index, type="pos")
        pos.mod <- matrix(0, ncol=15000, nrow=length(ladder))
        for(k in 1:15000){pos.mod[,k] <- sort(sample(roxy$pos, size=length(ladder), replace=FALSE), decreasing=FALSE)}
        dd <- apply(pos.mod, 2, function(x5, ladder) {cor(x5, ladder)}, ladder)
        v <- which(dd == max(dd))[1]
        v2 <- which(roxy$pos %in% pos.mod[, v])
        roxy <- list(pos = pos.mod[, v], hei = roxy$hei[v2], wei = ladder)
      }else{
        mi=length(ladder)
        if(mi > length(roxy$pos)){
          print(paste("ERROR!! using the initial threshold you specified we did not find enough peaks for", who,", please try reducing the 'init.thresh' or check the plot for this plant using 'detect.ladder' function, maybe this sample did not have ladder :)"))
          roxy <- list(pos=seq(1,length(ladder)) + rnorm(length(ladder),0,1), hei=seq(1,length(ladder))+ rnorm(length(ladder),0,1), wei= ladder)
        }else{
          pos.mod <- combn(roxy$pos, m=mi)
          dd <- apply(pos.mod, 2, function(x5,ladder){cor(x5,ladder)}, ladder)
          v <- which(dd == max(dd))
          v2 <- which(roxy$pos %in% pos.mod[,v])
          roxy <- list(pos=pos.mod[,v], hei=roxy$hei[v2], wei=ladder) 
        }
      }
    }
    # END OF METHOD "cor"
    #######################
    #####################
    #####################
    #####################
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
          plot(x, ylim=c(0,mm+se2.upp), type="l", main=attributes(x)$mycomm, cex.main=cex.title)
        }
        #roxy <- list(pos=roxy$pos, hei=roxy$hei, wei= ladder)
        roxy <- list(pos=seq(1,length(ladder)) + rnorm(length(ladder),0,1), hei=seq(1,length(ladder))+ rnorm(length(ladder),0,1), wei= ladder)
      }else{roxy <- list(pos=roxy$pos[s1:s2], hei=roxy$hei[s1:s2], wei= ladder)}
    }
    ################
    # once ladder is found
    ##############################################################
    if(method == "cor" | method == "ci"){
      if(draw == TRUE){
        xx <- roxy$wei
        yy <- roxy$pos
        mod <- lm(yy~xx)
        
        xlabels <- as.vector(predict(mod, newdata=data.frame(xx=seq(0,max(ladder),by=25))))
        plot(x, type="l", main=attributes(x)$mycomm,cex.main=cex.title, xlab="", ylab="Intensity", xaxt="n", col=transp("grey39",0.6), xlim=c(0,(max(roxy$pos)+100)), ylim=c(-100,max(roxy$hei)))
        if(method == "cor"){
          legend("topleft", legend=paste("Correlation=",round(max(dd)), sep=""), bty="n") 
        }
        axis(1, at=xlabels, labels=seq(0,max(ladder),by=25), cex.axis=0.9)
        abline(v=roxy[[1]], col="red", lty=3)
        if(method == "ci"){
          abline(h=mm, col="blue", lty=3)
          abline(h=(mm+se2.upp), col="blue", lty=3)
          abline(h=(mm-se2.low), col="blue", lty=3)
          legend("topright", legend=c("90% CI", "Peaks selected"), col=c("blue", "red"), bty = "n", lty=c(3,3), cex=1)
        }else{legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", lty=c(3), cex=1)}
        text(x=roxy[[1]], rep(-200, length(ladder)), labels=ladder, cex=0.6)
        #text(x=0, mm, labels="90% CI", cex=0.6)
      }
    }
    return(roxy)
    #########
  }
