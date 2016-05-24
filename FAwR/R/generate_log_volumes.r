#	$Id: generate_log_volumes.r 4005 2010-06-19 19:08:00Z hamannj $	


## this script will generate the log volumes
## and return a data.frame that contains the
## original stem list and the colums of the
## log volumes appended on to the end

## > sample.13$plants
##    plot sp.code         d6          dbh       tht        cr n.stems   expf
## 1     1      PP 10.6968804  8.032913840 30.188408 0.7053927       1 100.00
## 2     1      DF  6.4655181  4.455463746 30.157367 0.9636109       1 100.00
## 3     2      IC  2.6718796  0.653022095  6.602279 0.8462617       1 100.00

## this function is 200 lines long.
## too long for direct insertion in book page. 
## this function will not replace already existing columns
## in a data.frame, but will append them
generate.log.vols <- function( x,
                              log.breaks=c(2,5,12,18,32,999),
                              log.grades=c("pulp","s4","s3","s2","s1","peeler"),
##                               log.breaks=c(5,15,30,75,999), ## cm
##                               log.grades=c("pulp","s3","s2","s1","peeler"),
                              display.stems=FALSE )
{
  
  if( !is.data.frame( x ) ) {
    stop( "x needs to be a data.frame" )
  }

  ## these are the diameter breaks for the log products
  ## for pulp, cns, small saw, large saw
##   ## you can change these to shrink the number of grades...
##   log.grades <- c("pulp","a","b","c","d","e","z" )
##   log.bks <- c(10,20,30,45,60,80,999)

  ## allocate the grade volumes data frame that will be
  ## cbind'd to the stems data frame (n1 -- for now)
  grade.vols <- as.data.frame( matrix( rep( 0, nrow(x) ),
                                      nrow(x),
                                      length(log.grades),
                                      byrow=TRUE ) )
  names( grade.vols ) <- log.grades

  ## code up the taper function in here...
  ## function is orginally in inches, but need to work in metric
  ## WAY TOO MUCH BUTT SWELL, MAYBE...
  ## inbound values are in inches
  ## need to convert to centimeters (dbh) and meters (H,hi)
  bcmof.diDBH <- function( dbh, tht, cr, hi ) {

    ##dbh <- dbh / 2.54 ## convert from cm to inches
    ##tht <- tht / 3.2808 ## convert to feet
    ##hi <- hi / 3.2808 ## convert to feet  

    dbh <- dbh * 2.54 ## convert from inches to cm
    tht <- tht * 0.3048 ## convert to meters
    hi <- hi * 0.3048 ## convert to meters 

    ## these numbers should be in what units again?
    ##print( sprintf( "dbh = %f, tht=%f, hi=%f", dbh, tht, hi ) )

    if( hi >= tht ) {
      retval <- 0.0
    } else {
      p <- 0.25    
      dib <- 1.02453 * dbh^0.88809 * 1.00035^dbh
      X <- ( 1.0 - sqrt( hi / tht ) ) / ( 1.0 - sqrt( p ) )
      Z = hi / tht
      a = 0.95086 * Z * Z;
      b = -0.18090 * log( Z + 0.001 );
      c = 0.61407 * sqrt( Z ) +  -0.35105 * exp( Z );
      d = 0.05686 * ( dbh / tht ); 
      retval <- ( dib * X^( a + b + c + d ) )
    }
    retval / 2.54 ## convert from cm to inches
  }
  
  
  ## the md is the merchantable diameter
  merch.height.func <- function( hi, dbh, tht, cr, md ) {
    dib <- bcmof.diDBH( dbh, tht, cr, hi )
    diff <- dib - md
    diff
  }
  
  ## smalian volume
  ## d1, d2, l in meters or feet
  ## return value is in cubic meters or feet
  smal.vol <- function( d1, d2, l ) {
    ##smal.vol <- pi * l * (d1*d1 + d2*d2 ) / 8.0 
    smal.vol <- 0.0054541539 * l * (d1*d1 + d2*d2 ) / 2
    smal.vol
  }
  
  logs <- NULL

  ##print( "starting to buck stem..." )  
  for( s in 1:nrow(x) ) {  
    ##print( x[s,] )    
    ## if dbh == 0, don't attempt to merch
    if( x[s,]$dbh <= 5.0 ) next
    
    ## for each log.brk
    ## determine if the diameter break is smaller
    ## than the diameter inside bark at the stump
    mh.bks <- rep( 0, length( log.breaks ) )
    for( i in 1:length(log.breaks) ) {      

      ##dbh <- x[s,]$DBH / 10.0
      ##tht <- x[s,]$THT

      ## do not convert. 
      dbh <- x[s,]$dbh 
      tht <- x[s,]$tht

      if( log.breaks[i] <= bcmof.diDBH( dbh, tht, cr, 0.0 ) ) {    
        mh <- uniroot( merch.height.func,
                      c(0,tht),
                      dbh=dbh,
                      tht=tht,
                      cr=0.60,
                      md=log.breaks[i])
        ##print( mh )
        mh.bks[i] <- mh$root
      } else {
        mh.bks[i] <- 0.0
      }
    }
    
    ##print( mh.bks )  

    ## create the final bucked stem
    ## with a stump height as the last entry
    vol.bks <- data.frame( log.grades, log.breaks, mh.bks )
    vol.bks.temp <- vol.bks[mh.bks > 0.3,] ## truncate to only those that have non-zero values
    dib <- bcmof.diDBH(dbh,tht,0.60,0.30 ) ## 1/3 feet stump height
    last.grade <- log.grades[nrow(vol.bks.temp)+1]
    last.grade.df <- data.frame(log.grades=last.grade,log.breaks=dib,mh.bks=0.30)
    
    ## prepare the "top" half of the grade data frame
    vol.bks <- rbind( vol.bks.temp, last.grade.df )
    
    ## now append those sorts that are zero
    vol.bks.null <- data.frame( log.grades, log.breaks, mh.bks )
    n.sorts <- nrow( vol.bks )
    diff.rows <- nrow( vol.bks.null ) - nrow( vol.bks )
    vol.bks <- rbind( vol.bks, tail( vol.bks.null, diff.rows ) )
    
    sed <- c(0,vol.bks[1:(n.sorts-1),]$log.breaks,rep(0,diff.rows))
    led <- c(vol.bks[1:n.sorts,]$log.breaks,rep(0,diff.rows))
    vol.bks <- cbind( vol.bks, sed, led )
    
    ## compute the smalian volumes for the 
    log.lens <- c( tht - vol.bks[1,]$mh.bks,
                  abs( diff( vol.bks$mh.bks ) )[1:n.sorts-1],
                  rep(0,diff.rows) )
    vol.bks <- cbind( vol.bks, log.lens )
    
    ##len <- abs( c( diff( vol.bks$mh.bks ) ) )
    ##vol.bks$sm.vol <- smal.vol( vol.bks$sed*0.01, vol.bks$led*0.01, vol.bks$log.lens )
    vol.bks$sm.vol <- smal.vol( vol.bks$sed, vol.bks$led, vol.bks$log.lens )
    rownames(vol.bks) <- 1:nrow(vol.bks)
    

    if( display.stems ) {
      ##print( "display.stems..." )
      hi <- 0:tht
      dbh <- rep(dbh,length(hi))
      tht <- rep(tht,length(hi))
      cr <- rep(0.6,length(hi))
      stem.tpr <- data.frame( cbind( dbh, tht, cr, hi ) )
      stem.tpr$dib <- NA
      for( i in 1:nrow(stem.tpr) ){  
        stem.tpr[i,]$dib <- bcmof.diDBH( 
                                        stem.tpr[i,]$dbh, 
                                        stem.tpr[i,]$tht, 
                                        stem.tpr[i,]$cr, 
                                        stem.tpr[i,]$hi )
        
      }
      
      ## plot the stem profile
      plot( stem.tpr$dib ~ stem.tpr$hi, type="l" )
      abline( v=vol.bks$mh.bks, lty=2, lwd=2 )
      abline( h=vol.bks$log.breaks, lty=3 )
      abline( h=0 ) 
      text( x=vol.bks$mh.bks + 1.5,
           y=vol.bks$log.breaks+1,
           labels=paste( vol.bks$log.breaks, "@", round( vol.bks$mh.bks, 2 ) ) )
    }    
      
        ## now, assign the volumes for each of the grades (sorts) to the temp data frame
    grade.vols[s,] <- vol.bks$sm.vol
  }
  
  ## add the grade.volumes to the stems (n1 - for now)
  ## you can't use merge because the initial data managers were sloppy
  ##grade.vols$CRUISEID <- n1$CRUISEID
  vol <- rowSums( grade.vols )
  x <- cbind( x, vol, grade.vols )
  
  x
}


################################################################################
## this is the new function for the species summaries.
## use the existing rconifers function definition.
## these grades are a little large and should be scaled down to include pulp 2-5"
sp.sums.2 <- function( x,
                      log.breaks=c(2,5,12,18,32,999),
                      log.grades=c("pulp","s4","s3","s2","s1","peeler") ) {

  if( class( x ) != "sample.data" ) {
    stop( "Error: x is not a sample.data object." )
    return
  }

  ## this will merchandise the plant list (sample.data$plants)
  ## and return a data.frame that contains a set of metrics
  x$plants <- generate.log.vols( x$plants, 
                                     log.breaks, 
                                     log.grades,
                                     display.stems=FALSE )
  
  x.by.sp <- split( x$plants, f=list(x$plants$sp.code ) )
  
  ## this local function computes the summaries each species.
  sp.sums.f <- function( x, npts ) {
    expf <- sum( x$expf * x$n.stems ) / npts
    tht <- sum( x$expf * x$n.stems * x$tht ) / sum( x$expf * x$n.stems )
    ba <- sum( x$expf * x$n.stems * x$dbh^2*0.0054541539 ) / npts
    expgtbh<-sum(x$n.stems[x$tht>4.5]*x$expf[x$tht>4.5])/npts
    qmd <- sqrt( ba / expgtbh / 0.0054541359 )

    ## add the number of trees per acre >= 7.0 inches
    tpa.7.plus <- sum(x$n.stems[x$dbh>=7]*x$expf[x$dbh>=7.0])/npts

    ## add the volume calculation here
    w <- as.matrix( x[,log.grades]  )
    x <- as.matrix( x$expf * x$n.stems )
    v <- t(x) %*% w / npts
    v.p <- v / sum( v )

    ## todo: add the total volume 
    ##total.vol <- sum( x$expf * x$n.stems ) / npts
    total.vol <- sum(v)
    
    sp.sums <- c(qmd,tht,ba,expf,tpa.7.plus,total.vol,v,v.p)
  }
  
  ##print( "made it here." )
  df <- as.data.frame( t(sapply( x.by.sp, sp.sums.f, nrow(x$plots) )) )
  ##names( df ) <- c("qmd","tht","ba","expf","expf.7.plus",log.grades)
  names( df ) <- c("qmd","tht","ba","expf","tpa.7.plus","sm.vol",
                   log.grades,paste( rep("p",length(log.grades)), log.grades, sep="." ))

  df
}




