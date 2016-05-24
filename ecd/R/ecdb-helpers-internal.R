### <======================================================================>
#
# This is an internal helper library for ecdb, primarily for bootstrap
#
### <---------------------------------------------------------------------->
.ecdb.epoch <- function() as.integer(Sys.time())

# a list to keep track of the long and short timers and pace of sleeps
.ecdb.timer <- function() list(short= .ecdb.epoch(), 
                               short.duration= 60,
                               short.sleep= 5,
                               long= .ecdb.epoch(),
                               long.duration= 60*12,
                               long.sleep= 60)

# pace the operation by sleeping after duration
# message is printed before the short sleep
# if type is null, invoke the standard two-level pace
.ecdb.pace <- function(timer, message, type=NULL) {
    tm <- .ecdb.epoch()
    if (is.null(type)) {
        timer <- .ecdb.pace(timer, message, "short")
        timer <- .ecdb.pace(timer, message, "long")
    } else if (type == "short") {
        duration <- timer$short.duration
        sleep    <- timer$short.sleep
        if (tm > timer$short + duration) {
            message()
            Sys.sleep(sleep*(tm-timer$short)/duration) # prorated sleep
            timer$short <- .ecdb.epoch()
        }
    } else if (type == "long") {
        # avoid overheating the CPU for your little laptop
        duration <- timer$long.duration
        sleep    <- timer$long.sleep
        if (tm > timer$long + duration) {
            print("Extended sleep...")
            Sys.sleep(sleep)
            timer$long <- .ecdb.epoch()
        }
    }
    timer
}
.ecdb.message <- function(object, action, conf, var) {
    rows <- .ecdb.row_count(object)
    name <- if("var" %in% names(conf)) conf$var else "alpha"
    print(paste(action, ":", name, "/", var, "step=", conf$step, 
                "rows=", rows, "; now is", Sys.time()))
}        

### <---------------------------------------------------------------------->
.ecdb.row_count <- function(object) {
    sql <- "SELECT COUNT(*) as cnt FROM ECDATTR a"
    cnt <- RSQLite::dbGetQuery(object@conn, sql)$cnt
    unname(cnt)
}
### <---------------------------------------------------------------------->
# parse the mpfr flag from conf
.ecdb.conf_use_mpfr <- function(conf) {
    mpfr <- FALSE
    if ("mpfr" %in% names(conf)) {
        if(conf$mpfr > 0) mpfr <- TRUE
    }
    return(mpfr)
}
### <---------------------------------------------------------------------->
# split pairs (list) to pair group (list of lists)
# so that mclapply can take a break once in a while
.ecdb.pair_group <- function(pairs, use.mpfr) {
    chunk <- if (use.mpfr) 10 else 1000
    N <- length(pairs)
    print(paste("total-length-of-pairs=", N, "chunk-size=", chunk,
                "; now is", Sys.time()
                ))
    idx <- ceiling(seq(1,N)/chunk)
    split(pairs, idx)
}
### <---------------------------------------------------------------------->

