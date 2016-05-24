### <======================================================================>
#
# This is an internal helper library for ecdb, primarily for bootstrap
#
### <---------------------------------------------------------------------->
# insert df to ECDATTR
setGeneric("ecdattr<-", function(object, value) standardGeneric("ecdattr<-"))

setMethod("ecdattr<-", "ecdb", function(object, value){

    # handle list
    if (is.list(value) | is.vector(value)) {
        for (p in value) {
            ecdattr(object) <- p
        }
        return(invisible(object))
    }

    # handle ecdattr object
    if (class(value) != "ecdattr") {
        stop("value must be a list or vector of ecdattr objects")
    }

    conn <- object@conn
    sql <- "REPLACE INTO ECDATTR ( 
                alpha_m, gamma_m, cusp, stdev, kurtosis,
                discr, jinv, ellipticity, const, time_stamp
            ) VALUES (
                :alpha_m, :gamma_m, :cusp, :stdev, :kurtosis,
                :discr, :jinv, :ellipticity, :const, :time_stamp
            )"
            
    a <- value@attr
    df <- data.frame(
        alpha_m = c(value@alpha_m),
        gamma_m = c(value@gamma_m),
        cusp    = c(value@cusp),
        stdev    = c(a$stdev),
        kurtosis = c(a$kurtosis),
        discr    = c(a$discr),
        jinv     = c(a$jinv),
        ellipticity = c(a$ellipticity),
        const = c(a$const),
        time_stamp  = c(a$time_stamp)
    )
    rs <- RSQLite::dbGetPreparedQuery(conn, sql, bind.data=df)
    ecdb.protectiveCommit(object)
    invisible(object)

})
### <---------------------------------------------------------------------->
