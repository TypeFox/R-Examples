## 2014-10-21 moved from mcompanion

                                                          # 2015-10-15 - default value for sign
sim_real <- function(abs, sign = rep(NA_real_, length(abs)), signprob = 0.5,
                              absgen = "runif", absarg = list(0,1), ...){
  abs[is.na(abs)] <- do.call(absgen, c(sum(is.na(abs)), absarg))

  sign[is.na(sign)] <- rbinom(sum(is.na(sign)),1,signprob)
  sign[sign==0] <- -1

  abs * sign                                    # note: result has same shape as abs and sign.
}                                               #  maybe it should be always a vector! ???

                                                           # 2015-10-15 - default value for arg
sim_complex <- function(abs, arg  = rep(NA_real_, length(abs)),
                        absgen = "runif", absarg = list(0,1),
                        arggen = runif,   argarg = list(-pi,pi), ... ){
  abs[is.na(abs)] <- do.call(absgen, c(sum(is.na(abs)), absarg))
  arg[is.na(arg)] <- do.call(arggen, c(sum(is.na(arg)), argarg))

  complex(modulus=abs, argument=arg)            # note: result is a vector
}

                                        # 2014-10-21: changed argument "eigval" and component
                                        #          "eigval" of the returned value to "values"
sim_numbers <- function(type = rep(as.character(NA), length(abs)),
                        abs  = rep(as.numeric(NA), length(type)),
                        sign = rep(as.numeric(NA), length(type)),
                        values = NULL,
                        ... ){
    stopifnot(is.character(type))    # new 2015-10-15

    if(is.null(values)){
        res <- rep(as.numeric(NA), length(type))   # 2014-06-01 was: as.complex(NA)
    }else{
        stopifnot(length(values) == length(type)) # new 2015-10-15

        flags <- !is.na(values)
        if(any(flags)){                                           # 2015-10-15 - modified
                                       # replace NA's in type with inferred types (from values)
            type[is.na(type) & flags & Im(values) == 0] <- "r"
            type[is.na(type) & flags & Im(values) != 0] <- "cp" # complex pair, not "c"

                                                              # stop if values contradicts type
            if(!all(type[flags & Im(values) == 0] == "r")  ||
               !all(type[flags & Im(values) != 0] == "c" |
                    type[flags & Im(values) != 0] == "cp" )
               ## all(!is.na(type))--ne proveryavam, mozhe tova da e polezno ponyakoga.
               ## (note on 2015-10-15: this is an old comment, below it is ignored!)
               ){
                stop("Argument 'values' is inconsistent with 'type' or 'abs'")
            }
        }

        res <- values
    }

    if(any(is.na(type)))    # new 2015-10-15 ; todo: check that they are in c("r","c","cp") ?
        stop(paste("Some NA elements in type could not be resolved:\n",
                   "\ttype = [", paste0(type, collapse=", "), "]"))

    sel <- type=="r" & is.na(res)
    res[sel] <- sim_real(abs[sel], sign[sel],...)

    sel <- (type=="cp" | type=="c") & is.na(res)
                   # 2014-06-01 conditional, to keep it real if no complex values are required
                   #            otherwise res will become complex even if length(sel)==0
                   # 2015-11-13 was:  if(length(sel) > 0) ...
                   #                  (doesn't do what is intended!)
    if(any(sel))
        res[sel] <- sim_complex(abs[sel], sign[sel], ...)

    list(values=res, type=type)
}
