###########################################################
# plotting tam expected scores curves
#..........................................................
plot.tam <- function(x, items=1:x$nitems, type="expected" ,
                     low=-3, high=3, ngroups=6, 
                     wle=NULL, export=TRUE, export.type="png", 
                     export.args=list(), observed=TRUE, overlay=FALSE , 
                     ask=FALSE, package="lattice" , 
					 fix.devices=FALSE , ...) {
#  requireNamespace("plyr")
time1 <- NULL
if ( fix.devices ){
  old.opt.dev <- base::getOption("device")
  old.opt.err <- c( base::getOption("show.error.messages"))
  old.par.ask <- graphics::par("ask")
  # remember new pars' values
  old.par.xpd <- graphics::par("xpd")
  old.par.mar <- graphics::par("mar")
  
  base::on.exit( base::options("device"=old.opt.dev))
  base::on.exit( base::options("show.error.messages"=old.opt.err), add=TRUE)
  base::on.exit( graphics::par("ask"=old.par.ask), add=TRUE)
  # restore new pars' values
  base::on.exit( graphics::par("xpd"=old.par.xpd), add=TRUE)
  base::on.exit( graphics::par("mar"=old.par.mar), add=TRUE)
}  
  
  tamobj <- x
  ndim <- tamobj$ndim
  tammodel <- "mml"
  if(is.null(ndim)) {
    ndim <- 1
    tammodel <- "jml"
  }
  if (ndim > 1 ) {
    if ( type=="expected"){
      stop ("Expected scores curves are only available for uni-dimensional models")
    }
  }
  
  nitems <- tamobj$nitems
  
  nnodes <- 100
  if (ndim == 1 ){
    theta <- matrix(seq(low, high, length=nnodes), nrow=nnodes, ncol=ndim)
  } else {
    #	theta <- tamobj$theta
    nnodes <- 40
    nodes <- seq(low, high, length=nnodes)
    theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )	
    nnodes <- nrow(theta)	
    B <- tamobj$B
  }
  
  iIndex <- 1:nitems
  A <- tamobj$A
  B <- tamobj$B
  if (tammodel == "mml") {
    xsi <- tamobj$xsi$xsi
  }  else {
    xsi <- tamobj$xsi
  }
  maxK <- tamobj$maxK
  resp <- tamobj$resp
  resp.ind <- tamobj$resp.ind
  resp[resp.ind==0] <- NA
  AXsi <- matrix(0,nrow=nitems,ncol=maxK )
  res <- calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                      nnodes=nnodes , maxK=maxK , recalc=TRUE )  
  rprobs <- res[["rprobs"]]
  AXsi <- res[["AXsi"]]
  cat <- 1:maxK - 1
  
  if ( type == "expected" ){
    
    expScore <- sapply(1:nitems, function(i) colSums(cat*rprobs[i,,], na.rm=TRUE))
    
    if (is.null(wle)) {
      if (tammodel == "mml") {
        wleobj <- tam.wle(tamobj)
        wle <- wleobj$theta
      } 
      else {
        wle <- tamobj$WLE    # model is jml
      }
    }
    wleSorted <- sort(wle, na.last=FALSE)
    ncases <- length(wleSorted)
    groupnumber <- round(seq(1:ncases) / (ncases/ngroups) + 0.5)
    
    aggr <- stats::aggregate(wleSorted, list(groupnumber), mean)
    theta2 <- aggr$x
    
    d <- data.frame(wle, resp)
    d1 <- d[order(wle),]
    d2 <- d1[-1]
    obScore <- apply(d2,2, function(x) stats::aggregate(x, list(groupnumber), mean, na.rm=TRUE))
  }
  # adds observed score for type="items"
  if (type == "items") {
    if (is.null(wle)) {
      if (tammodel == "mml") {
        wleobj <- tam.wle(tamobj)
        wle <- wleobj$theta
      } 
      else {
        wle <- tamobj$WLE    # model is jml
      }
    }
    wleSorted <- sort(wle, na.last=FALSE)
    ncases <- length(wleSorted)
    # sometimes this returns one more group than requested
    groupnumber <- round(seq(1:ncases) / (ncases/ngroups) + 0.5)  
    groupnumber[groupnumber > ngroups] = ngroups # prevents it
    
    aggr <- stats::aggregate(wleSorted, list(groupnumber), mean)
    theta2 <- aggr$x
    
    d <- data.frame(wle, resp)
    d1 <- d[order(wle),]
    d2 <- d1[-1]
    obScore <- lapply(d2, function(item) {
      comp_case = stats::complete.cases(item)
      item = item[comp_case]
      uniq_cats = sort(unique(item))
      plyr::ldply(split(item, groupnumber[comp_case]), .id = "group", 
                  function (group) {
        ngroup = length(group)
        cat_freq = list()
        for (catt in uniq_cats) {
          cat_freq[[paste0("cat_", catt)]] = sum(group == catt)/ngroup
        }
        data.frame(cat_freq)
      })
    })
  }
  
  #*************************************************
  # begin plot function
  
  for (i in (1:nitems)[items]) {
    #***********************************************************
    #** expected item response curves	
    if ( type=="expected"){
      if (i==1 || !overlay) {
        ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )
        graphics::plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability",
             #         main=paste("Expected Scores Curve - Item ", i)
             main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] )	 ,
             ylim=ylim2 , ...
        )
      } else {
        graphics::lines(theta, expScore[,i],type="l", col=i, lwd=3, pch=1) 
      }
      if (observed) {
        graphics::lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
      }
    }
    #***********************************************************
    if ( type=="items"){	
      
      rprobs.ii <- rprobs[i,,]
      rprobs.ii <- rprobs.ii[ rowMeans( is.na(rprobs.ii) ) < 1 , ]		
      K <- nrow(rprobs.ii)		
      if ( ndim == 1 ){ theta0 <- theta }
      dat2 <- NULL
      #************
      if ( ndim > 1 ){
        B.ii <- B[i,,]	
        ind.ii <- which( colSums( B.ii ) > 0 )[1]
        rprobs0.ii <- rprobs.ii
        rprobs0.ii <- stats::aggregate( t(rprobs0.ii) , list( theta[,ind.ii] ) , mean )
        theta0 <- rprobs0.ii[,1,drop=FALSE]
        rprobs.ii <- t( rprobs0.ii[,-1] )						
      }
      #**************
      for (kk in 1:K){
        # kk <- 1
        dat2a <- data.frame( "Theta" = theta0[,1] , "cat" = kk , "P" = rprobs.ii[kk,] )
        dat2 <- rbind(dat2 , dat2a)
      }
      main <- paste("Item", colnames(x$resp)[i] )				
      auto.key <- NULL				
      simple.key <- paste0("Cat" , 1:K -  1)
      auto.key <- simple.key
      dat2$time <- dat2$cat
      dat2$time1 <- paste0("Cat" , dat2$time )
      
      simple.key <- FALSE
      Kpercol <- K
      # floor(K/Kpercol)+1
      #**************************************
      # package lattice
      if ( package=="lattice"){
        auto.key <- list(  lines=TRUE , points=FALSE , rows=2)
        h1 <- lattice::xyplot(P ~ Theta, dat2, group = time1, type = 'l', auto.key = auto.key,
                              main = main, ylim = c(-0.1,1.1), simple.key = simple.key , 
                              xlim=c(low,high) , 
                              ylab = expression(P(theta)), xlab = expression(theta), ... ) 
        graphics::plot(h1)	
      }	
      #**************************************
      # package graphics
      if ( package == "graphics" ){
        kk <- 1
        dfr <- dat2
        dfr1a <- dfr[ dfr$cat == kk , ]
        # setting larger margins on the right to fit the legend
        if(graphics::par("mar")[4] < 6.1) {
          graphics::par("mar" = c(old.par.mar[1:3], 6.1))
        }
        
        graphics::plot( dfr1a$Theta , dfr1a$P , ylim=c(-.1,1.1) , 	
              ylab = expression(P(theta)), xlab = expression(theta) ,
              col=kk+1 , type="l" , main=main, xpd = TRUE, ...
        )
        for (kk in seq(2,K) ){
          dfr1a <- dfr[ dfr$cat == kk , ]		
          graphics::lines( dfr1a$Theta , dfr1a$P , col=kk+1 )					
          # graphics::points( dfr1a$Theta , dfr1a$P , pch=kk , col=kk+1 )	
          
        }
        if(observed) {
          obScore_it = obScore[[i]]
          for(kk in seq(1:K)) {
            graphics::lines(theta2, obScore_it[, kk + 1], col = kk + 1, lty = 2)
            graphics::points(theta2, obScore_it[, kk + 1], col = kk + 1)
          }
        } 
        # puts legend outside the plot region
        # only expected
        legend_entry = paste0("Cat" , seq(0,K-1), "(exp.)")
        pch = NA_integer_
        col = 1 + 1:K
        lty = 1
        # adds observed
        if (observed) {
          legend_entry = c(legend_entry, paste0("Cat" , seq(0,K-1), "(obs.)"))
          legend_entry = sort(legend_entry)
          col = rep(1 + 1:K, each = 2)
          pch = c(NA_integer_, 1)
          lty = c(1, 2)
        }
        graphics::legend(high + 0.3, 1.1, legend_entry,
               cex = 0.7,
               pch = pch, col= col,
               horiz = FALSE, lty= lty, bty = "n", xpd = TRUE)
      }
      
      
      #***************************************
      
    }
    #***************		
    graphics::par(ask=ask)	
  }             # end item ii
  #*************************************************
  
  #*****
  # export item plots
  if (export) {
    
    if(!file.exists("Plots")) dir.create( "Plots" )
    export.type.dev <- switch(export.type,
                              "ps"="postscript",
                              "emf"=if (.Platform$OS.type == "windows") "win.metafile" else "x11",
                              "wmf"=if (.Platform$OS.type == "windows") "win.metafile" else "x11",
                              export.type)
    export.type.ff <- switch(export.type,
                             "postscript"="ps",
                             "win.metafile"="wmf",
                             "x11"="wmf",
                             export.type)
    
    
    options(show.error.messages = FALSE)
    options("device"=export.type.dev)
    
    for (i in (1:nitems)[items]) {
      
      
      itemlab <- colnames(tamobj$resp)[i]
      dev.err <- try({ 
        do.call("dev.new", 
                args=list("filename"=file.path("Plots", paste("Item_", itemlab, ".", export.type.ff, sep="")), 
                          export.args))
      })
      
      if(!is.null(dev.err)){
        warning( dev.err[1], "  --> No file created."  )
      }else{
        
        #***************************************************
        # expected response functions
        if (type=="expected"){
          ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )        
          graphics::plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability", 
               main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] ) ,
               ylim=ylim2 , ... )        
          if (observed ) {
            graphics::lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
          }
        }
        if ( type=="items" ){
          
          rprobs.ii <- rprobs[i,,]
          rprobs.ii <- rprobs.ii[ rowMeans( is.na(rprobs.ii) ) < 1 , ]		
          K <- nrow(rprobs.ii)		
          if ( ndim == 1 ){ theta0 <- theta }
          dat2 <- NULL
          #************
          if ( ndim > 1 ){
            B.ii <- B[i,,]	
            ind.ii <- which( colSums( B.ii ) > 0 )[1]
            rprobs0.ii <- rprobs.ii
            rprobs0.ii <- stats::aggregate( t(rprobs0.ii) , list( theta[,ind.ii] ) , mean )
            theta0 <- rprobs0.ii[,1,drop=FALSE]
            rprobs.ii <- t( rprobs0.ii[,-1] )						
          }
          #**************
          for (kk in 1:K){
            # kk <- 1
            dat2a <- data.frame( "Theta" = theta0[,1] , "cat" = kk , "P" = rprobs.ii[kk,] )
            dat2 <- rbind(dat2 , dat2a)
          }
          main <- paste("Item", colnames(x$resp)[i] )				
          auto.key <- NULL				
          simple.key <- paste0("Cat" , 1:K -  1)
          auto.key <- simple.key
          dat2$time <- dat2$cat
          dat2$time1 <- paste0("Cat" , dat2$time )
          
          simple.key <- FALSE
          Kpercol <- K
          # floor(K/Kpercol)+1
          if ( package=="lattice"){
            auto.key <- list(  lines=TRUE , points=FALSE , rows=2)
            h1 <- lattice::xyplot(P ~ Theta, dat2, group = time1, type = 'l', auto.key = auto.key,
                                  main = main, ylim = c(-0.1,1.1), simple.key = simple.key , 
                                  xlim=c(low,high) , 
                                  ylab = expression(P(theta)), xlab = expression(theta), ... ) 
            graphics::plot(h1)  
          }	
          #**************************************
          # package graphics
          if ( package == "graphics" ){
            kk <- 1
            dfr <- dat2
            dfr1a <- dfr[ dfr$cat == kk , ]
            # setting larger margins on the right to fit the legend
            if(graphics::par("mar")[4] < 6.1) {
              graphics::par("mar" = c(old.par.mar[1:3], 6.1))
            }
            
            graphics::plot( dfr1a$Theta , dfr1a$P , ylim=c(-.1,1.1) , 	
                  ylab = expression(P(theta)), xlab = expression(theta) ,
                  col=kk+1 , type="l" , main=main, xpd = TRUE, ...
            )
            for (kk in seq(2,K) ){
              dfr1a <- dfr[ dfr$cat == kk , ]		
              graphics::lines( dfr1a$Theta , dfr1a$P , col=kk+1 )					
              # graphics::points( dfr1a$Theta , dfr1a$P , pch=kk , col=kk+1 )	
              
            }
            if(observed) {
              obScore_it = obScore[[i]]
              for(kk in seq(1:K)) {
                graphics::lines(theta2, obScore_it[, kk + 1], col = kk + 1, lty = 2)
                graphics::points(theta2, obScore_it[, kk + 1], col = kk + 1)
              }
            } 
            # puts legend outside the plot region
            # only expected
            legend_entry = paste0("Cat" , seq(0,K-1), "(exp.)")
            pch = NA_integer_
            col = 1 + 1:K
            lty = 1
            # adds observed
            if (observed) {
              legend_entry = c(legend_entry, paste0("Cat" , seq(0,K-1), "(obs.)"))
              legend_entry = sort(legend_entry)
              col = rep(1 + 1:K, each = 2)
              pch = c(NA_integer_, 1)
              lty = c(1, 2)
            }
            graphics::legend(high + 0.3, 1.1, legend_entry,
                   cex = 0.7,
                   pch = pch, col= col,
                   horiz = FALSE, lty= lty, bty = "n", xpd = TRUE)
          }          
        }
        
        grDevices::dev.off(dev.cur())
        
      }
      
             options("device"=old.opt.dev)
             options(show.error.messages = as.character(old.opt.err))        
    }
    
    #*****
    # Print path
    if(is.null(dev.err)){ cat("....................................................\n",
                              "Plots exported in", export.type, "format into folder:\n", 
                              file.path(getwd(), "Plots")) ; 
							  utils::flush.console() }
  }
  
}

plot.tam.mml <- plot.tam
plot.tam.jml <- plot.tam

