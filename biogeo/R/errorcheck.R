errorcheck <-
function (world, dem, dat, countries = "", countryfield = "NAME", 
            vars = c("bio1", "bio12", "bio5", "bio6"), res = 10,elevc="",diff=50) 
  {
    nr <- nrow(dat)
    rid <- 1:nr
    cn <- names(dat)
    cnw <- names(world)
    fieldsmissing(dat, fields = c("ID", "x", "y", "x_original", 
                                  "y_original", "Correction", "Exclude","Reason"))
    xn <- coord2numeric(dat$x)
    yn <- coord2numeric(dat$y)
    xy <- cbind(xn, yn)
    f0 <- which(xn == 0 & yn == 0)
    x1 <- (abs(xn) > 180) * 1
    y1 <- (abs(yn) > 90) * 1
    if(any(is.na(x1))){x1[is.na(x1)]<-1}
    if(any(is.na(y1))){y1[is.na(y1)]<-1}
    if (any(x1 + y1 > 0)) {
      ff <- which(x1 + y1 > 0)
      er <- str_c(dat$ID[ff], sep = "", collapse = ";")
      msg <- paste("coordinate values impossible see ID:", 
                   as.character(er))
      stop(msg)
    }
    g <- SpatialPoints(xy)
    world@proj4string = CRS(as.character(NA))
    s1 <- over(g, world)
    wc <- match(countryfield, names(s1))
    country_ext <- s1[, wc]
    rstVals <- extract(dem, g)
    f <- which(is.na(rstVals))
    wrongEnv <- rep(0, nr)
    wrongEnv[f] <- 1
    if (nchar(countries) == 0) {
      CountryMismatch <- rep(0, nr)
    } else {
      fcnt <- match(countries, cn)
      CountryMismatch <- rep(0, nr)
      country_orig <- dat[, fcnt]
      m <- rep(-1, nr)
      for (i in 1:nr) {
        m[i] <- match(country_orig[i], country_ext[i], nomatch = 2)
      }
      CountryMismatch <- (m - 1)
    }
    xydat <- data.frame(x = xn, y = yn)
    z3 <- precisioncheck(xydat, x = "x", y = "y", 10, 30)
    lowprec <- z3$preci
    r <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = res/60)
    r[] <- 1:ncell(r)
    cid <- cellFromXY(r, xy)
    cid1 <- data.frame(Species=dat$Species, cid)
    dups <- duplicated(cid1) * 1
    cid2 <- data.frame(Species=dat$Species, cid, dups)
    us <- unique(cid2$Species)
    nv <- length(vars)
    vrs <- rep(0, nv)
    for (i in 1:nv) {
      vrs[i] <- match(vars[i], cn)
      }
    if (any(is.na(vrs))){
      msg<-"variable is missing from dataset"
      stop(msg)
    }
    env <- dat[, vrs]
    ee <- {
    }
    for (i in 1:nv) {
      e <- outliers(rid, cid2$Species, cid2$dups, env[, i])
      ee <- cbind(ee, e)
    }
    nv2 <- nv * 2
    jj <- as.data.frame(ee[, seq(1, nv2, 2)])
    ex <- as.data.frame(ee[, seq(2, nv2, 2)])
    nj <- paste(vars, "_j", sep = "")
    ne <- paste(vars, "_e", sep = "")
    names(jj) <- nj
    names(ex) <- ne
    
    ##########
    if (nchar(elevc) == 0) {
      elevMismatch <- rep(0, nr)
      demElevation<-rep(NA,nr)
    } else {
    fieldsmissing(dat, fields = elevc) # does the elevc column exist?
    elev<-elevcheck(dat, dem, elevc=elevc,diff=diff)
    elevMismatch<-elev$elevMismatch
    demElevation<-elev$demElevation
    }
    
    error <- (rowSums(cbind(CountryMismatch, wrongEnv, lowprec, 
                            jj, ex,elevMismatch)) > 1) * 1
    spperr <- error * 0
    f <- which(error == 1)
    sp <- unique(dat$Species[f])
    for (i in 1:length(sp)) {
      f <- which(dat$Species == sp[i])
      spperr[f] <- 1
    }
    dat1 <- data.frame(dat, cellid = cid, dups, country_ext, 
                       CountryMismatch, wrongEnv, lowprec = lowprec, ex, jj,elevMismatch,demElevation, error, 
                       spperr)
    return(dat1)
  }
