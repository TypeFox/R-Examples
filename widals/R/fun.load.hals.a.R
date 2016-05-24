
if( getRversion() >= "2.15.1" ) { utils::globalVariables(
    c("GP", "run.parallel", "widals", "Ht", "Hst.ls", "Z", "Z.na", "Hs", "b.lag", "train.rng", "test.rng", "cv",
    "FUN.GP", "FUN.MH", "FUN.I", "FUN.EXIT", "rm.ndx", "locs", "lags", "xgeodesic", "stnd.d", "ltco",
    "rho.upper.limit", "rgr.lower.limit", "d.alpha.lower.limit"
    )) }

fun.load.hals.a <-
function() {
    
    xenvr <- as.environment(1)
    
    if( run.parallel ) {
        
        sfExport("Z", "Hs", "Ht", "Hst.ls", "b.lag", "train.rng", "test.rng")
        
        suppressWarnings(sfLibrary(widals))
        
    }
    
    p.ndx.ls <- list( c(1,2) )
    assign( "p.ndx.ls", p.ndx.ls, pos=xenvr )
    
    f.d <- list( dlog.norm, dlog.norm, dlog.norm, dlog.norm, dlog.norm )
    assign( "f.d", f.d, pos=xenvr )
    
    
    FUN.MH <- function(jj, GP.mx, X) {
        
        Z.als <- Hals.snow( jj, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, b.lag=b.lag, GP.mx=GP.mx )
        
        resids <- Z - Z.als
        our.cost <- sqrt( mean( resids[ train.rng, ]^2 ) )
        
        return( our.cost )
    }

    assign( "FUN.MH", FUN.MH, pos=xenvr )
    
    #FUN.GP <- NULL
    
    
    FUN.GP <- function(GP.mx) {
        GP.mx[ GP.mx[ , 1] > rho.upper.limit, 1 ] <- rho.upper.limit
        GP.mx[ GP.mx[ , 2] < rgr.lower.limit, 2 ] <- rgr.lower.limit
        return(GP.mx)
        
    }
    assign( "FUN.GP", FUN.GP, pos=xenvr )
    
    
    FUN.I <- function(envmh, X) {
        cat( "Improvement ---> ", envmh$current.best, " ---- " , envmh$GP, "\n" )
        
    }
    assign( "FUN.I", FUN.I, pos=xenvr )
    
    
    
    FUN.EXIT <- function(envmh, X) {
        
        GP.mx <- matrix(envmh$GP, 1, length(envmh$GP))
        
        Z.als <- Hals.snow( 1, Z=Z, Hs=Hs, Ht=Ht, Hst.ls=Hst.ls, b.lag=b.lag, GP.mx=GP.mx )
        
        resids <- Z - Z.als
        our.cost <- sqrt( mean( resids[ test.rng, ]^2 ) )
        
        cat( envmh$GP, " -- ", our.cost,   "\n" )
        
        assign( "Z.als", Z.als, pos=xenvr )
        assign( "our.cost", our.cost, pos=xenvr )
        assign( "GP", envmh$GP, pos=xenvr )
        cat( paste( "GP <- c(", paste(format(GP,digits=5), collapse=", "), ") ### ", format(our.cost, width=6), "\n", sep="" ) )
        
    }
    
    assign( "FUN.EXIT", FUN.EXIT, pos=xenvr )
    
}
