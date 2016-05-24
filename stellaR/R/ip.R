# Isochrones construction by mean of tracks interpolation at given ages

makeIso <- function(age, z=NULL, y=NULL, ml=NULL, afe=NULL, log=FALSE, linear=TRUE, tr=NULL, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
    
    if(linear) {
        fun <- approx
        meth <- "linear"
    } else {
        stop("sorry, LINEAR=FALSE is not yet implemented")
#        fun <- spline
#        meth <- "natural"
    } 
    
    if(is.null(tr)) {
        if( any(c(is.null(z), is.null(y), is.null(ml), is.null(afe)))) {
            stop("please provide z, y, ml, afe")
        }
        
        v.mass <- seq(0.30, 1.10, by=0.05)
        nmass <- length(v.mass)
        
                                        # get all the set of masses
                                        # for given composition
        tr <- getTrkSet(v.mass, z, y, ml, afe, baseURL)
        if(any(is.na(tr))) {
          warning("CDS is unavailable; please try later")
          return(NULL)
        }          
    }
    else {
        if( !all(c(is.null(z), is.null(y), is.null(ml), is.null(afe)))) 
            warning("discarding z, y, ml, afe from arguments list")
        if( ! all(class(tr) == c("trkset", "stellar"))) 
            stop("tr must be of class \"trkset\"")
        
        nmass <- length(tr)
        
        ztmp <- sapply(tr, function(x) x$z)
        ytmp <- sapply(tr, function(x) x$y)
        mltmp <- sapply(tr, function(x) x$ml)
        afetmp <- ifelse(sapply(tr, function(x) x$alpha.enh) > 0,1,0)
        
        test <- length(unique(ztmp)) + length(unique(ytmp)) + length(unique(mltmp)) + length(unique(afetmp))
        
        if(test == 4) {
            z <- ztmp[1]
            y <- ytmp[1]
            ml <- mltmp[1]
            afe <- afetmp[1]
        }
        else {
            stop("input tracks are disomogeneous in (z, y, ml, afe)")
        }
    }
    
                                        # points on the tracks
    length.d <- sapply(tr, function(x) dim(x$data)[1])
    max.length <- max(length.d)

    if( min(age) < 7 | max(age) > 15 ) {
      stop("the allowed age range is [7 - 15] Gyr")
    }
                                        # sort age in increasing order.
                                        # Convert to yr
    age <- sort(age) * 1e9
    
    time <- NULL
    logTe <- NULL
    logL <- NULL
    
                                        # extract information about mass and
                                        # sort in increasing order
    mass <- sapply(tr, function(x) x$mass)
    ii <- order(mass)
    mass <- sort(mass)
    
                                        # extract information about logTe and
                                        # logL (fill short tracks with NA)
    for(i in ii) {
        time <- cbind(time, c(tr[[i]]$data$time, rep(NA, max.length-length.d[i])))
        logTe <- cbind(logTe, c(tr[[i]]$data$logTe, rep(NA, max.length-length.d[i])))
        logL <- cbind(logL, c(tr[[i]]$data$logL, rep(NA, max.length-length.d[i])))
    }
    
                                        # initial point for isochrone
    min.age <- min(age)
    initial <- which(time[,1] > log10(min.age))[1]
    
                                        # Use log yr or yr for interpolation
    if(!log) {
        time <- 10^time
        age <- age
        ageiso <- age
    } else {
        time <- time
        ageiso <- age
        age <- log10(age)
    }
    
                                        # first point (on the track with
                                        # lower mass)
    Mi <- rep(mass[1], length(age))
    Li <- approx(time[,1], logL[,1], xout=age)$y
    Ti <- approx(time[,1], logTe[,1], xout=age)$y
    
                                        # interpolate tracks at given ages
    iso.tmp <- NULL
    iso.tmp <- rbind(iso.tmp, c(Li, Ti, Mi))
   
    Mi <- sapply(initial:max.length, function(i)
      approx(x=time[i,], y=mass, xout=age, method=meth)$y)
    Li <- sapply(initial:max.length, function(i)
      approx(x=time[i,], y=logL[i,], xout=age, method=meth)$y)
    Ti <- sapply(initial:max.length, function(i)
      approx(x=time[i,], y=logTe[i,], xout=age, method=meth)$y)

    if(length(age) > 1)
      iso.tmp <- rbind(iso.tmp, cbind(t(Li), t(Ti), t(Mi)))
    else
      iso.tmp <- rbind(iso.tmp, cbind(Li, Ti, Mi))
    
                                        # fill the result list, discarding NA
    iso <- list()
    nage <- length(age)
    for(i in 1:nage) {
        data <- cbind(iso.tmp[,i],iso.tmp[,nage+i], iso.tmp[,2*nage+i])
        data <- as.data.frame(data)
                                        # calculate R and log g
        data$radius <- 3.3338e7*sqrt(10^data[,1])/((10^data[,2])^2)
        data$logg <- log10(2.74496e4*data[,3]/(data$radius^2))

                                        # column names
        colnames(data) <- c("logL", "logTe", "mass", "radius", "logg")

                                        # discard NA
        data <- data[!is.na(data$mass),]
                                        # fill iso object
        iso[[i]] <- list(age=ageiso[i]*1e-9, z=z, y=y, ml=ml, alpha.enh=ifelse(afe>0,0.3,0), data=data)
        class(iso[[i]]) <- c("iso", "stellar")
    }
    class(iso) <- c("isoset", "stellar")

                                        # return the appropriate object
    if(nage > 1) 
      return(iso)
    else
      return(iso[[1]])
}


###########################################################

interpTrk <- function(z, y, ml, afe, vmass=seq(0.30,1.10, by=0.05), baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
    
                                        # check boundary
    test <- z >= 0.0001 & z <= 0.01
    test <- c(test, y >= 0.25 & y <= 0.42)
    test <- c(test, ml >= 1.70 & ml <= 1.90)
    test <- c(test, afe == 0 | afe == 1)
    verify <- all(test)
    
    if(!verify) 
        stop("please check inputs parameters (z, y, ml, afe)")

    vz <- c( (1:9)*1e-4, (1:9)*1e-3, 1e-2 )
    vy <- c(0.25, 0.27, 0.33, 0.38, 0.42)
    vml <- c(1.7, 1.8, 1.9)

    # check if ml is in the DB
    is.ml <- ml %in% vml
    if(is.ml) {
        ml.ip <- vml[which(vml == ml)]
    } else {
      # if not in the DB, trap the value
        pos.ml <- which(ml < vml)[1]
        ml.ip <- c(vml[pos.ml-1], vml[pos.ml])
    }

    # check if y is in the DB
    is.y <- y %in% vy
    if(is.y) {
        y.ip <- vy[which(vy == y)]
    } else {
      # if not in the DB, trap the value
        pos.y <- which(y < vy)[1]
        y.ip <- c(vy[pos.y-1], vy[pos.y])
    }

    # check if z is in the DB
    is.z <- z %in% vz  
    if(is.z) {
        z.ip <- vz[which(vz == z)]
    } else {
      # if not in the DB, trap the value
        pos.z <- which(z < vz)[1]
        z.ip <- c(vz[pos.z-1], vz[pos.z])
    }

    # grid of cases selected for interpolation
    grid <- expand.grid(ml=ml.ip, y=y.ip, z=z.ip)
    # number of row (it can be 1, 2, 4, 8)
    n.grid <- dim(grid)[1]
    
    if(n.grid == 1) {
      # no interpolation!
        tr <- getTrkSet(vmass, grid[1,3], grid[1,2], grid[1,1], afe, baseURL)
        return(tr)
    }
    
                                        # get tracks ...
    tr <- lapply(1:n.grid, function(i) with(grid[i, ], getTrkSet(vmass, z, y, ml, afe, baseURL)))
    nmass <- length(vmass)

    if(all(is.na(tr[[1]]))) {
      warning("CDS is unavailable; please try later")
      return(NULL)
    }
                                            # transform z to log10(z)
    z <- log10(z)
    z.ip <- log10(z.ip)
    
                                        # selctor of variable for interpolation
                                        # TRUE if interpolation occour on the var
    sel <- c(!is.ml, !is.y, !is.z)
    n.ip <- sum(sel)
    v.ip <- c(ml, y, z)[sel]
    which.v <- (1:3)[sel]
    
                                        # interpolation procedure
    for(run in 1:n.ip) {
                                        # index of first set
        ii <- seq(1, n.grid, by=2^run)
        for(i in ii) {
                                        # shift: second case from first one
            shift <- 2^(run-1)
            
            X.tgt <- v.ip[1]
            X.0 <- grid[i,which.v[run]]
            X.1 <- grid[i+shift,which.v[run]]

                                        # interpolate for all the masses
            for(j in 1:nmass) {
                V0 <- cbind(tr[[i]][[j]]$data$time, tr[[i]][[j]]$data$logTe, tr[[i]][[j]]$data$logL)
                V1 <- cbind(tr[[i+shift]][[j]]$data$time, tr[[i+shift]][[j]]$data$logTe, tr[[i+shift]][[j]]$data$logL)
                
                tgt <- V0 + (V1-V0)/(X.1-X.0)*(X.tgt-X.0)
                tr[[i]][[j]]$data$time <- tgt[,1]
                tr[[i]][[j]]$data$logTe <- tgt[,2]
                tr[[i]][[j]]$data$logL <- tgt[,3]
            }
        }
        v.ip <- v.ip[-1]
    }
    
                                        # set (z, y, ml) for interpolated tracks
    tro <- lapply(tr[[1]], function(i) {
      i$z <- 10^z
      i$y <- y
      i$ml <- ml
      i$data <- i$data[,1:5]
      return(i);
    } )
    class(tro) <- c("trkset", "stellar")
    return(tro)
    
}
