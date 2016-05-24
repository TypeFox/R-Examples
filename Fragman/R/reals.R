reals <-
  function(x, panel=c(100:400), shi=1, ploidy=2, left.cond=c(0.4,3), right.cond=0.2, window=0.5){
    # where x is a list with positions, heights and weights
    # panel is a vector with possible values
    # get rid of the negative base pair sizes
    neg <- which(x$wei <= 0)
    if(length(neg) >0){
      x <- list(pos=x$pos[-neg], hei=x$hei[-neg], wei=x$wei[-neg])
    }
    ##############
    picos <- numeric()
    for(d in 1:length(x$wei)){
      respi <- which(abs(panel - x$wei[d]) <= window)
      if(length(respi) > 0){
        picos[d] <- 1
      }else{picos[d] <- 0}
    }
    z1 <- which(picos == 1)
    # z1 <- which(x$wei < panel[2] & x$wei > panel[1])
    if(length(z1) >0){
      x2 <- list(pos=x$pos[z1], hei=x$hei[z1], wei=x$wei[z1]) 
      ###############################################
      # if peaks are closer than 1 bp just select the tallest one and get rid of the other one
      x3 <- separate(x2, shi, type="bp")
      ################################################
      ### eliminate peaks that are not at least half the size of the tallest peak
      # but this should only happen for peaks before the tallest peak
      highest <- which(x3$hei == max(x3$hei))
      che <- x3$hei[1:highest[1]]
      ha <- which(che >= (x3$hei[highest] * left.cond[1]))
      cha <- x3$wei[1:highest[1]] 
      ha2 <- unique(c(which(abs(cha - x3$wei[highest]) >= left.cond[2]), highest)) # are further than, greater than
      # rigth side condition
      chu <- x3$hei[highest[1]:length(x3$hei)]
      hu <- (highest[1]:length(x3$hei))[which(chu > (max(x3$hei) * right.cond))]
      # intersection
      ss1 <- intersect(ha, ha2)
      ha3 <- unique(c(ss1, hu))
      if(length(ha3) > 0){
        x3 <- list(pos = x3$pos[ha3], hei = x3$hei[ha3], wei = x3$wei[ha3])
      }else{x3 <- x3}
      ### now lets select the peaks that are only after the tallest
      #ho <- which(x3$hei > (max(x3$hei)*right.cond))
      #if(length(ho) > 0){
      #  x3 <- list(pos = x3$pos[ho], hei = round(x3$hei[ho]), wei = x3$wei[ho])
      #}else{x3 <- x3}
      ########
      #he <- which(x3$hei > (max(x3$hei)*whenheight))
      #x3 <- list(pos = x3$pos[he], hei = round(x3$hei[he]), 
      #          wei = x3$wei[he])
      
      ##########################################################################
      ### BASED ON PLOIDY SELECT HOW MANY PEAKS YOU WANT
      ###############################################
      # if more than 2 peaks select the bigger peaks
      z2 <- length(x3$pos)
      ## one peak find, then duplicate
      if(z2 == 1){
        x4 <- list(pos=rep(x3$pos,ploidy), hei=rep(x3$hei,ploidy), wei=rep(x3$wei,ploidy)) 
      }
      # more than one peak found, we need to find 
      if(z2 > 1){
        # find all different combinations possible
        toget <- sort(x3$hei, decreasing=T)[1:ploidy]
        z3 <- which(x3$hei %in% toget)
        x4 <- list(pos=x3$pos[z3], hei=round(x3$hei[z3]), wei=x3$wei[z3]) 
      }#else{x4 <- x3}
    } else{x4 <- list(pos=rep(0, ploidy), hei=rep(0, ploidy), wei=rep(0,ploidy)) }
    ###############################################
    return(x4)
    
  }
