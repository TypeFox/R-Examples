#' Bootstrap data for the Elliptic DB (ECDB)
#' 
#' Main interface to generate data for ECDB based on the configuration.
#'
#' @method bootstrap ecdb
#'
#' @param object an object of ecdb class.
#' @param action the action operating on the ecdb. 
#' @param skip.existing logical, if \code{TRUE} (default), skip if action already done in history.
#' 
#' @return Row count.
#'
#' @keywords ecdb
#'
#' @export
#'
#' @importFrom utils str
#'
#' 
#' @include ecdb-class.R
#'
### <======================================================================>
"bootstrap.ecdb" <- function(object, action="all", skip.existing=TRUE)
{
    conn <- object@conn
    if (object@is.internal) {
        if (action != "all") {
            stop("Internal DB: Only default action=all is allowed")
        }
        # this is all we can do to fit a 1MB sqlite db
        conf <- object@conf$rect$internal
        rs <- .ecdb.bootstrap_rect(object, "internal", conf)    
        history(object) <- "bootstrap-rect-internal"
        return(rs)
    }
    
    if (action %in% names(object@conf$line)) {
        hkey <- paste("bootstrap", action, sep="-")
        if (skip.existing & hkey %in% history(object)) {
            print(paste("skip existing:", hkey))
            return(0)
        }
        conf <- object@conf$line[[action]]
        rs <- .ecdb.bootstrap_line(object, action, conf)       
        history(object) <- hkey
        return(rs)
    }
 
    act <- strsplit(action, "-")[[1]]
    if (act[1] == "rect" & act[2] %in% names(object@conf$rect)) {
        hkey <- paste("bootstrap", action, sep="-")
        if (skip.existing & hkey %in% history(object)) {
            print(paste("skip existing:", hkey))
            return(0)
        }
        conf <- object@conf$rect[[ act[2] ]]
        rs <- .ecdb.bootstrap_rect(object, action, conf)       
        history(object) <- hkey
        return(rs)
    }
    if (act[1] == "polar" & act[2] %in% names(object@conf$polar)) {
        hkey <- paste("bootstrap", action, sep="-")
        if (skip.existing & hkey %in% history(object)) {
            print(paste("skip existing:", hkey))
            return(0)
        }
        conf <- object@conf$polar[[ act[2] ]]
        rs <- .ecdb.bootstrap_polar(object, action, conf)       
        history(object) <- hkey
        return(rs)
    }
    
    if (action == "all") {
        rs <- 0
        lines <- names(object@conf$line)
        for (action in lines) {
            rs <- rs + bootstrap(object, action)
        }
        history(object) <- paste("bootstrap", "all", "lines", sep="-")
        
        rects <- paste("rect", names(object@conf$rect), sep="-")
        for (action in rects) {
            if (action == "rect-internal") next
            rs <- rs + bootstrap(object, action)
        }
        history(object) <- paste("bootstrap", "all", "rects", sep="-")
 
        polars <- paste("polar", names(object@conf$polar), sep="-")
        for (action in polars) {
            rs <- rs + bootstrap(object, action)
        }
        history(object) <- paste("bootstrap", "all", "polars", sep="-")
        
        history(object) <- paste("bootstrap", "all", "actions", sep="-")
        return(rs)
    }
    
    stop(paste("Unknown ecdb action:", action))
}
### <---------------------------------------------------------------------->
#' @rdname bootstrap.ecdb
setGeneric("bootstrap", function(object, action="all", skip.existing=TRUE) standardGeneric("bootstrap"))
#' @rdname bootstrap.ecdb
setMethod("bootstrap", signature("ecdb"), bootstrap.ecdb)
### <---------------------------------------------------------------------->
.ecdb.bootstrap_rect <- function(object, action, conf) {
    
    use.mpfr <- .ecdb.conf_use_mpfr(conf)
    if (use.mpfr) print("Flag use.mpfr is set to TRUE")

    rs <- 0
    a_arr <- seq(conf$from, conf$to, by=conf$step)
    print(paste("bootstrap_rect: a=seq(", conf$from, conf$to, "by", conf$step,")"))
    print(paste("bootstrap_rect: a-length=", length(a_arr), "; now is", Sys.time()))

    timer <- .ecdb.timer()    
    for (a in a_arr) {
        pairs <- ecdattr.pairs(c(a), seq(conf$from, conf$to, by=conf$step), 
                               use.mpfr=use.mpfr)
        for (pg in .ecdb.pair_group(pairs, use.mpfr)) {        
            rs <- rs + write(pg, object)
            message <- function() .ecdb.message(object, action, conf, a) 
            timer <- .ecdb.pace(timer, message)
        }
    }
    return(rs) # row count
}
### <---------------------------------------------------------------------->
.ecdb.bootstrap_polar <- function(object, action, conf) {
    
    rs <- 0
    conf1 <- conf$log2R
    conf2 <- conf$angle
    log2R_arr <- seq(conf1$from, conf1$to, by=conf1$step)
    angle_arr <- seq(conf2$from, conf2$to, by=conf2$step)
    print(paste("bootstrap_polar: log2R=seq(", conf1$from, conf1$to, "by", conf1$step,")"))
    print(paste("bootstrap_polar: angle=seq(", conf2$from, conf2$to, "by", conf2$step,")"))
    print(paste("bootstrap_polar: log2R-length=", length(log2R_arr), "; now is", Sys.time()))
    print(paste("bootstrap_polar: angle-length=", length(angle_arr), "; now is", Sys.time()))
    
    timer <- .ecdb.timer()
    use.mpfr <- .ecdb.conf_use_mpfr(conf)
    if (use.mpfr) print("Flag use.mpfr is set to TRUE")
    for (a in angle_arr) {
        pairs <- ecdattr.pairs_polar(R=2^log2R_arr, theta=c(a/180*pi), 
                                     use.mpfr=use.mpfr)
        for (pg in .ecdb.pair_group(pairs, use.mpfr)) {        
            rs <- rs + write(pg, object)
            message <- function() .ecdb.message(object, action, conf2, a) 
            timer <- .ecdb.pace(timer, message)
        }
    }
    return(rs) # row count
}
### <---------------------------------------------------------------------->
.ecdb.bootstrap_line <- function(object, action, conf) {

    use.mpfr <- .ecdb.conf_use_mpfr(conf)
    if (use.mpfr) print("Flag use.mpfr is set to TRUE")
    
    arr <- seq(conf$from, conf$to, by=conf$step)
    var_name <- if("var" %in% names(conf)) conf$var else "alpha"
    var_val  <- if("val" %in% names(conf)) conf$val else NaN
    cusp  <- if("cusp" %in% names(conf)) conf$cusp else 0

    print(paste("bootstrap:", action, "seq(", conf$from, conf$to, "by", conf$step,")"))
    print(paste("var_name=", var_name, "var_val=", var_val, "cusp=", cusp))

    pairs <- NULL
    
    if (var_name == "critical") {
        pairs <- ecdattr.pairs(arr, NaN, cusp, use.mpfr)
        print("special pairing: critical")
    } else if (var_name == "gamma") {
        pairs <- ecdattr.pairs(arr, c(var_val), cusp, use.mpfr)
    } else if (var_name == "alpha") {
        pairs <- ecdattr.pairs(c(var_val), arr, cusp, use.mpfr)
    } else if (var_name == "4x") {
        pairs1 <- ecdattr.pairs(c(0),  arr, cusp, use.mpfr)
        pairs2 <- ecdattr.pairs(c(0), -arr, cusp, use.mpfr)
        pairs3 <- ecdattr.pairs( arr, c(0), cusp, use.mpfr)
        pairs4 <- ecdattr.pairs(-arr, c(0), cusp, use.mpfr)
        pairs <- rbind(pairs, pairs1)
        pairs <- rbind(pairs, pairs2)
        pairs <- rbind(pairs, pairs3)
        pairs <- rbind(pairs, pairs4)
    }
    
    N <- length(pairs)
    pairs_group <- .ecdb.pair_group(pairs, use.mpfr)
    
    rs <- 0
    timer <- .ecdb.timer()
    for (pg in pairs_group) {
        if (! is.list(pg)) {
            stop(paste("The pair group is not a list:", str(pg)))
        }
        rs <- rs + write(pg, object)

        x1 <- pg[[1]]
        progress <- if (var_name == "gamma") x1@alpha  else x1@gamma
        message <- function() .ecdb.message(object, action, conf, 
                                            paste(progress,",", rs,"/",N))
        timer <- .ecdb.pace(timer, message)
    }
    return(rs) # row count
}
### <---------------------------------------------------------------------->
