testComposition <- function(Z, Y, ML, AFE) {
# test for the existence of a given combination of
# chemical composition (Z, Y, AFE) and mixing length

    data(compositions)

    afe <- compositions$afe
    ml <- compositions$ml
    z <- compositions$z
    y <- compositions$y
    
    testAFE <- AFE %in% afe
    if(!testAFE) 
        warning("[alpha/Fe] value not in the database")
    
    testML <- ML %in% ml
    if(!testML) 
        warning("mixing-length value not in the database")
    
    testZ <- Z %in% z
    if(!testZ) 
        warning("z value not in the database")
    
    if( testZ ) {
        testY <- Y %in% y[z == Z,]
    } else {
        testY <- FALSE
    }
    if(!testY) 
        warning("y value not in the database")
       
    return( testAFE & testML & testZ & testY )
}

showComposition <- function() {
# Show the possible combinations of
# chemical composition (Z, Y, AFE) and mixing length

    data(compositions)

    afe <- compositions$afe
    ml <- compositions$ml
    z <- compositions$z
    y <- compositions$y

    cat("Mixing-length values:\n")
    cat("\t", paste(ml,collapse=", "), "\n\n")
    
    cat("alpha-enhancement values:\n")
    cat("\t", paste(afe,collapse=", "), " (i.e. [alpha/Fe] = 0.0 [alpha/Fe] = 0.3)", "\n\n")
    
    cat("Chemical compositions:\n")
    df <- as.data.frame(cbind(z,y))
    print(df, print.gap=2, row.names=rep("    ", dim(df)[1]))
}

getZahb <- function(z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve ZAHB for given parameters

    if(baseURL == "ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
      msg <- "CDS is unavailable; please try later"
    } else {
      msg <- "data file not found; please check path"
    }
    
    specificURL <- "zahb/ZAHB_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")
    
    if( !testComposition(z, y, ml, afe))
        stop("required data not present in the database")
    
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    URL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")
  
    DATA <- tryCatch(
              read.table(URL, skip=6),
              error=function(e) NULL,
              warning=function(e) NULL
              )
    
    if(is.null(DATA)) {
      warning(msg)
      return(NA)
    }
             
    colnames(DATA) <- c("mass", "logTe", "logL")
    L <- list(z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("zahb", "stellar")
    return(L)
}

getTrk <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from PMS to He flash) for given parameters

    if(baseURL == "ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
      msg <- "CDS is unavailable; please try later"
    } else {
      msg <- "data file not found; please check path"
    }
    specificURL <- "trk/TRK_Z"

    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")

    M <- seq(0.3, 1.1, by=0.05)
    if( !testComposition(z, y, ml, afe) | ! (as.integer(100*m) %in% as.integer(100*M)) )
        stop("required data not present in the database")
    
    M <- format(m, nsmall=2)
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "/", sep="")
    
    URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")

    DATA <- tryCatch(
              read.table(URL, skip=5),
              error=function(e) NULL,
              warning=function(e) NULL
              )
    
    if(is.null(DATA)) {
      warning(msg)
      return(NA)
    }

    colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
    L <- list(mass=m, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("trk", "stellar")
    return(L)
}

getIso <- function(age, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve iso (from PMS to He flash) for given parameters
  
    if(baseURL == "ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
      msg <- "CDS is unavailable; please try later"
    } else {
      msg <- "data file not found; please check path"
    }
    specificURL <- "iso/ISO_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")

    AGE <- seq(8,15,by=0.5)
    if( !testComposition(z, y, ml, afe) | ! (as.integer(age*100) %in% as.integer(AGE*100)))
        stop("required data not present in the database")
    
    AGE <- age*1000
    zero <- ifelse(AGE < 10000, "0", "")
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "/", sep="")
    
    URL <- paste(dirURL, "AGE", zero, AGE, "_Z", Z, "_He", Y, "_ML", ML, AFE, ".DAT", sep="")

    DATA <- tryCatch(
              read.table(URL, skip=6),
              error=function(e) NULL,
              warning=function(e) NULL
              )    
    
    if(is.null(DATA)) {
      warning(msg)
      return(NA)
    }
    
    colnames(DATA) <- c("logL" ,"logTe", "mass", "radius", "logg")
    L <- list(age=age, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("iso", "stellar")
    return(L)
}


getHb <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from ZAHB to thermal pulses) for given parameters
  
    if(baseURL == "ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
      msg <- "CDS is unavailable; please try later"
    } else {
      msg <- "data file not found; please check path"
    }
    specificURL <- "hb/TRK_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")

    M <- seq(0.3, 1.1, by=0.05)
    if( !testComposition(z, y, ml, afe) | ! (as.integer(m*100) %in% as.integer(M*100)) )
        stop("required data not present in the database")
    
    M <- format(m, nsmall=2)
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "_HB/", sep="")
  
  # search the mass of RGB progenitor...
    data(masshb)
    T <- c(m, z, y, ml)
    idx <- apply(masshb[,1:4], 1, function(x) all(as.numeric(x) == as.numeric(T)))
    masshb.ext <- masshb[idx,]
    sel <- masshb.ext[, 5] == substr(AFE, 2, nchar(AFE))
    massRGB <- format(masshb.ext[sel, 6], nsmall=4)
    
    URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, "_ZAHB", massRGB, ".DAT", sep="")

    DATA <- tryCatch(
              read.table(URL, skip=5),
              error=function(e) NULL,
              warning=function(e) NULL
          )

    if(is.null(DATA)) {
          warning(msg)
          return(NA)
        }
    
    colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
    L <- list(mass=m, massRGB=m, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
    class(L) <- c("hb", "stellar")
    return(L)
}


getHbgrid <- function(z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
# retrieve track (from ZAHB to thermal pulses) for given parameters

    if(baseURL == "ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {
      msg <- "CDS is unavailable; please try later"
    } else {
      msg <- "data file not found; please check path"
    }
    specificURL <- "hb/TRK_Z"
    
    if(substr(baseURL, nchar(baseURL), nchar(baseURL)) != "/")
        baseURL <- paste(baseURL, "/", sep="")

    if( !testComposition(z, y, ml, afe) )
        stop("required data not present in the database")
    
    Z <- format(z, nsmall=5, scientific=FALSE)	   
    Y <- format(y, nsmall=4)	   
    ML <- format(ml, nsmall=2)	   
    AFE <- ifelse(afe, "_AS09a3", "_AS09a0")
    
    dirURL <- paste(baseURL, specificURL, Z, "_He", Y, "_ML", ML, AFE, "_HB/grid/", sep="")
  
  # search the mass of RGB progenitor...
    data(masshbgrid)
    T <- c(z, y, ml)
    idx <- apply(masshbgrid[,2:4], 1, function(x) all(as.numeric(x) == as.numeric(T)))
    masshb.ext <- masshbgrid[idx,]
    sel <- masshb.ext[, 5] == substr(AFE, 2, nchar(AFE))
    massRGB <- format(masshb.ext[sel, 6], nsmall=4)
    M <- format(masshb.ext[1, 1], nsmall=2)
    
    n.trk <- length(massRGB)
    L <- list()

  # txt progress bar
    cat("Download in progress...\n")
    pb <- txtProgressBar()

    for(i in 1:n.trk) {
        URL <- paste(dirURL, "OUT_M", M, "_Z", Z, "_He", Y, "_ML", ML, AFE, "_ZAHB", massRGB[i], ".DAT", sep="")

        DATA <- tryCatch(
                  read.table(URL, skip=5),
                  error=function(e) NULL,
                  warning=function(e) NULL
                  )
        if(is.null(DATA)) {
          warning(msg)
          return(NA)

        }
        setTxtProgressBar(pb, i/n.trk)
        
        colnames(DATA) <- c("mod", "time", "logL" ,"logTe", "mass", "Hc", "logTc", "logRHOc", "MHEc", "Lpp", "LCNO", "L3a", "Lg", "radius", "logg")
        L[[i]] <- list(mass=round(as.numeric(massRGB[i]),2), massRGB=M, z=z, y=y, ml=ml, alpha.enh=ifelse(afe,0.3,0), data=DATA)
        class(L[[i]]) <- c("hb", "stellar")
    }
    close(pb)
    
    class(L) <- c("hbset", "stellar")
    return(L)
}


###############################


getTrkSet <- function(m, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {

    pars <- expand.grid(m=m, z=z, y=y, ml=ml, afe=afe)
    N <- nrow(pars)

  # test the compositions of the grid
    ok <- sapply(1:N, function(i) with(pars[i, ], testComposition(z, y, ml, afe)))
    
 
    
    if(sum(ok) > 1) {
         # txt progress bar
        cat("Download in progress...\n")
        pb <- txtProgressBar()
      
        trk <- lapply((1:N)[ok], function(i) {
          # update the progress bar
          setTxtProgressBar(pb, i/sum(ok))
          return(with(pars[i, ], getTrk(m, z, y, ml, afe, baseURL)))
        } )
        close(pb)
      }
    else
    # only 1 track ...
        return(with(pars[1, ], getTrk(m, z, y, ml, afe, baseURL)))
  
    if(sum(ok) > 1)
        class(trk) <- c("trkset", "stellar")
    else
    # no results
        return(NA)
  
    return(trk)
}


getIsoSet <- function(age, z, y, ml, afe, baseURL="ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/540/A26/") {

    pars <- expand.grid(age=age, z=z, y=y, ml=ml, afe=afe)
    N <- nrow(pars)

  # test the compositions of the grid
    ok <- sapply(1:N, function(i) with(pars[i, ], testComposition(z, y, ml, afe)))

    if(sum(ok) > 1) {
       # txt progress bar
        cat("Download in progress...\n")
        pb <- txtProgressBar()
        
        iso <- lapply((1:N)[ok], function(i) {
          # update the progress bar
          setTxtProgressBar(pb, i/sum(ok))
          
          return(with(pars[i, ], getIso(age, z, y, ml, afe, baseURL)))
        } )
        close(pb)
      }
    else
    # only 1 iso ...
        return(with(pars[1, ], getIso(age, z, y, ml, afe, baseURL)))
  
    if(sum(ok) > 1)
        class(iso) <- c("isoset", "stellar")
    else
    # no results
        return(NA)
  
    return(iso)
}


