check_section <- function(fn,
                          tpldat,
                          tplsec,
                          R_list,
                          check,
                          bounds=NULL,
                          data_type,
                          phase=NULL,
                          re,
                          mcmcpars,
                          profpars,
                          secname,
                          objfunname="f",
                          intcheck=c("strict","sloppy","none")) {
    intcheck <- match.arg(intcheck)
    Rnames  <- names(R_list)
    msg <- ""
    if (check!="write") {
        info <- tpldat$info[[tplsec]]
        tplnames <- info$vname
        if (length(setdiff(tplnames,Rnames))>0) {
            msg <- paste("missing values in list ",
                         "(",secname," section): ",
                         paste(setdiff(tplnames,Rnames),sep=","),sep="")
        } else if (length(setdiff(Rnames,tplnames))==0 && !all(tplnames==Rnames)) {
            msg <- "all values present, but order doesn't match"
        } else {
            msg <- ""
            if (length(grep("\\.",tplnames)>1)) {
                msg <- paste(msg,"dots in parameter/variable names")
            }
            ## FIXME: attach() throws a NOTE in R CMD check
            ##  no longer remember what it's doing:
            ##  figure it out and document, and try to
            ##   find another way around it ...
            ## I believe the issue is that we need to try
            ## to evaluate the various parameter dimensions in
            ## the context of R_list (see eval(parse(...)) calls
            ## below), but otherwise allow lookup to fall back
            ## to the regular search path; there is probably a better
            ## way to handle this
            attach(R_list,name="R_list",warn.conflicts=FALSE)
            on.exit(detach(R_list),add=TRUE)
            ## now need to check dimensions etc...
            for (i in 1:nrow(info)) {
                v <- info[i,]
                ## x <- get(v$vname)
                ##  i.e. search only in data, not everywhere ...
                x <- R_list[[v$vname]]
                v$type <- gsub("init_","",v$type)
                if (v$type %in% c("int","ivector","imatrix")) {
                    if (any(trunc(x)!=x)) msg <- paste(msg,v$vname,
                                 "non-integer;")
                }
                if (v$type %in% c("int","number")) {
                    if (length(x)>1) msg <- paste(msg,"length(",v$vname,
                                                  ")>1;")
                }
                if (v$type %in% c("ivector","vector")) {
                    if (any(is.na(c(v$X1,v$X2)))) {
                        msg <- paste(msg,"NAs in dimensions in ",v$vname)
                    } else {
                        tpllen <- eval(parse(text=paste(v$X2,"-",v$X1)))+1
                        if (length(x)!=tpllen)
                            msg <- paste(msg,"length mismatch in ",v$vname,
                                         ": ",length(x)," (r) != ",tpllen," (tpl)",
                                         sep="")
                    }
                }
                if (v$type %in% c("imatrix","matrix")) { 
                    tpldim <- with(v,c(
                                       eval(parse(text=paste(v$X2,"-",v$X1)))+1,
                                       eval(parse(text=paste(v$X4,"-",v$X3)))+1))
                    rdim <- dim(x)
                    if (is.null(rdim)) {
                        msg <- paste(msg,v$vname,"not a matrix;")
                    } else {
                        if (any(rdim!=tpldim))
                            msg <- paste(msg,"dimension mismatch in ",v$vname,
                                         ": (",rdim[1],",",rdim[2],"), (r) != ",
                                         " (",tpldim[1],",",tpldim[2],") (tpl)",
                                         sep="")
                    } ## dimensions
                } ## if matrix
                if (length(grep(v$type,"array"))>0) {
                    arraydim <- as.numeric(substr(v$type,1,1))
                    tpldim <- numeric(arraydim)
                    for (j in 1:arraydim) {
                        tpldim[j] <-
                            eval(parse(text=paste(v[2*j+2],"-",v[2*j+1])))+1
                    }
                    rdim <- dim(x)
                    if (is.null(rdim)) {
                        msg <- paste(msg,v$vname,"not an array;")
                    } else {
                        if (any(rdim!=tpldim))
                            msg <- paste(msg,"dimension mismatch in ",v$vname,
                                         ": (",
                                         paste(rdim,sep=","),"), (r) != ",
                                         " (",
                                         paste(tpldim,sep=","),") (tpl)",
                                         sep="")
                    } ## !is.null(rdim)
                } ## array
            } ## loop over variables
        } ## checking
        return(msg)
    } else { ## check=="write"
        ## WRITING THE TPL SECTION
        ## FIXME: handle bounds, phases, and 1D matrices??
        ## bounds: list OR (named) 2-col matrix or 2xn matrix
        ## phases: list OR (named) vector or n-vector
        ## arrays must be 1-based: if you want zero-based arrays
        ##     then write the TPL section yourself
        ## if (!missing(phase))
        ##    stop("phase support not yet implemented")
        ## no existing section: need title line
        sectitle <- paste(secname,"_SECTION",sep="")
        nvals <- length(R_list)
        if (length(grep("\\.",names(R_list)>1))) {
            warning("dots changed to underscores in parameter/variable names")
            names(R_list) <- gsub("\\.","_",names(R_list))
        }
        objstr <- NULL
        pars <- (secname=="PARAMETER")  ## only check bounds for parameters
        if (pars) {
            if (is.null(objfunname)) stop("must specify a name for the objective function")
            objstr <- indent(paste("objective_function_value",objfunname))
        }
        parstr <- character(nvals)
        ## FIXME: parameter tables should be more standardized --
        ## room for all possible dimensions, bounds, phase?
        partab <- data.frame(type=character(nvals),
                             vname=names(R_list),
                             dim1=rep(NA,nvals),
                             dim2=rep(NA,nvals),
                             dim3=rep(NA,nvals),
                             dim4=rep(NA,nvals),
                             lower=rep(NA,nvals),
                             upper=rep(NA,nvals),
                             phase=rep(NA,nvals),
                             stringsAsFactors=FALSE)
        for (i in 1:length(R_list)) {
            x <- R_list[[i]]
            n <- names(R_list)[i]
            ## attempt to coerce (if not all numeric, will end up as character and stop ...)
            errmsg <- paste("(param #",
                            i,": ",n,")",sep="")
            if (is.data.frame(x)) x <- as.matrix(x)
            ## attempt to detect whether a variable is an integer or not
            ##  (too complicated!)
            is.int <- function(x,n=NULL,debug=FALSE) {
                opt_in <- (!is.null(data_type) && !is.null(n) &&
                           n %in% names(data_type)[data_type=="integer"])
                opt_out <- (!is.null(data_type) && !is.null(n) &&
                            n %in% names(data_type)[data_type!="integer"])
                trunc_test <- all(trunc(x)==x)
                mode_test <- storage.mode(x)=="integer"
                if (debug) cat("opt_in",opt_in,"opt_out",opt_out,
                               "trunc_test",trunc_test,"mode_test",mode_test,"\n")
                opt_in || (!opt_out &&
                           !intcheck=="none" &&
                           ((intcheck=="strict" && mode_test) ||
                            (intcheck=="trunc" && trunc_test)))
            }
            if (is.int(x,n) && !is.null(dim(x)) &&
                length(dim(x))>2) {
                stop("can't handle integer arrays of dimension>2",
                     errmsg)
            }
            ## FIXME: bounds and phase support not implemented for integers ... ?
            ##  (does it make sense?)
            if (is.int(x,n)) {
                if (length(x)==1 && is.null(dim(x))) {
                    parstr[i] <- paste("init_int",n)
                    partab$type[i] <- "int"
                } else if (length(x)>1 && is.null(dim(x))) {
                    parstr[i] <- paste("init_ivector ",n,"(1,",length(x),")",sep="")
                    partab$type[i] <- "ivector"
                    partab$dim1 <- length(x)
                } else if (!is.null(dim(x)) && length(dim(x))==2) {
                    parstr[i] <- paste("init_imatrix ",n,
                                       "(1,",dim(x)[1],",1,",dim(x)[2],")",sep="")
                    partab$dim1[i] <- dim(x)[1]
                    partab$dim2[i] <- dim(x)[2]
                }
            } else if ((!is.null(data_type) &&
                        n %in% names(data_type)[data_type=="numeric"]) ||
                       storage.mode(x) %in% c("numeric","double") ||
                       (storage.mode(x)=="integer" && intcheck=="none"))
            {
                bpstr <- function(pars,bounds=NULL,phase=NULL,n) {
                    pstr <- NULL
                    boundstr <- ""
                    if (pars && !is.null(bounds) && n %in% names(bounds)) {
                        boundstr <- "bounded_"
                        pstr <- paste0(bounds[[n]][1],",",bounds[[n]][2])
                    }
                    if (pars && !is.null(phase) && n %in% names(phase)) {
                        pstr <- paste(c(pstr,phase[[n]]),collapse=",")
                    }
                    list(boundstr=boundstr,pstr=pstr)
                }
                pfun <- function(x) if (is.null(x)) "" else paste0("(",x,")")
                bp <- bpstr(pars,bounds,phase,n)
                if (length(x)==1 && is.null(dim(x))) {
                    parstr[i] <- paste0("init_",bp$boundstr, "number ",n, pfun(bp$pstr))
                    partab$type[i] <- "number"
                } else if (length(x)>1 && is.null(dim(x))) {
                    parstr[i] <- paste0("init_",bp$boundstr,"vector ",n,
                                        pfun(paste(c("1",length(x),bp$pstr),collapse=",")))
                    partab$type[i] <- "vector"
                    partab$dim1[i] <- length(x)
                } else if (!is.null(dim(x)) && length(dim(x))==2) {
                    parstr[i] <- paste0("init_",bp$boundstr,"matrix ",n,
                                        pfun(paste(c("1",dim(x)[1],"1",dim(x)[2],bp$pstr),collapse=",")))
                    partab$type[i] <- "matrix"
                    partab$dim1 <- dim(x)[1]
                    partab$dim2 <- dim(x)[2]
                } else if (!is.null(dim(x)) && length(dim(x))>2) {
                    ndim <- length(dim(x))
                    if (ndim>7) stop("can't handle arrays of dim>7",errmsg)
                    parstr[i] <- paste0("init_",bp$boundstr,ndim,"darray",n,
                                        pfun(paste(c(c(rbind(rep(1,ndim),dim(x))),bp$pstr),collapse=",")))
                    partab$type[i] <- "array"
                    ## FIXME: store array dimensions?
                    if (!is.null(mcmcpars))
                        stop("arrays currently incompatible with MCMC",errmsg)
                } ## multi-dim array
            } else {
                stop("can only handle numeric values",errmsg)
            }
        } ## loop over R list
        parstr <- indent(parstr)
        cursec <- tpldat$secs[[secname]]
        if (!is.null(cursec)) {
            cursec <- cursec[-1] ## drop title
            cursec <- grep("^ *$",cursec,invert=TRUE,value=TRUE)  ## drop blank lines
        }
        restr <- mcmcstr <- profstr <- NULL
        if (pars) {
            ## deal with random effects vectors
            if (!is.null(re)) {
                redim <- sapply(re,length)
                re_vectors <- re[redim==1]
                re_mats <- re[redim==2]
                restr <- ""
                if (length(re_vectors)>0)
                    restr <- c(restr,indent(paste("random_effects_vector ",names(re_vectors),
                                                  "(1,",re_vectors,")",sep="")))
                if (length(re_mats)>0)
                    restr <- c(restr,indent(paste("random_effects_matrix ",names(re_mats),
                                                  "(1,",
                                                  sapply(re_mats,"[",1),",1,",
                                                  sapply(re_mats,"[",2),")",sep="")))
            }
            ## FIXME: uuuuuugly! need a better, more consistent way
            ## of handling parameter attributes ...
            make_names <- function(z,pref1="sdreport_",pref2="r_") {
                if (z %in% partab$vname) { ## in newly specified parameters
                    i <- match(z,partab$vname)
                    type <- partab$type[i]
                    name <- partab$vname[i]
                    dimvals <- na.omit(unlist(partab[i,c("dim1","dim2","dim3","dim4")]))
                    dimvals <- c(rbind(rep(1,length(dimvals)),dimvals))
                } else if (z %in% tt$vname) { ## in existing parameters
                    i <- match(z,tt$vname)
                    type <- tt$type[i]
                    name <- tt$vname[i]
                    dimvals <- na.omit(unlist(tt[i,3:7]))
                }
                indent(paste(pref1,type," ",pref2,name,
                             if (length(dimvals)==0) "" else
                             paste("(",paste(dimvals,collapse=","),")",sep=""),
                             sep=""))
            }
            if(!is.null(mcmcpars) || !is.null(profpars)) {
                tt <- tpldat$info$other
                allnames <- c(partab$vname,tt$vname)
                if (!is.null(mcmcpars)) {
                    bad <- which(!mcmcpars %in% allnames)
                    if (length(bad)>0) {
                        stop("some mcmcpars not found in parameter table:",
                             paste(mcmcpars[bad],collapse=", "))
                    }
                    mcmcstr <- sapply(mcmcpars,make_names,pref1="sdreport_",pref2="r_")
                }
                if (!is.null(profpars)) {
                    bad <- which(!profpars %in% allnames)
                    if (length(bad)>0) {
                        stop("some profpars not found in parameter table:",
                             paste(profpars[bad],collapse=", "))
                    }
                    profstr <- sapply(profpars,make_names,pref1="likeprof_",pref2="p_")
                }
            }
        }
        return(c(sectitle,"",objstr,parstr,cursec,restr,mcmcstr,profstr))
    }
}
