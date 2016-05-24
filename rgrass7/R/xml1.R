parseGRASS <- function(cmd, legacyExec=NULL) {
    cmdCACHE <- get("cmdCACHE", envir=.GRASS_CACHE)
    res <- cmdCACHE[[cmd]]
    if (is.null(legacyExec))
        legacyExec <- get.legacyExecOption()
    stopifnot(is.logical(legacyExec))
    if (get("SYS", envir=.GRASS_CACHE) != "unix" && !legacyExec) {
        warning("legacyExec set TRUE on non-unix platform")
        legacyExec <- TRUE
    }
    if (is.null(res)) {
        ext <- get("addEXE", envir=.GRASS_CACHE)
        WN_bat <- get("WN_bat", envir=.GRASS_CACHE)
        if (get("SYS", envir=.GRASS_CACHE) == "WinNat" && nchar(WN_bat) == 0) {
            WN_bat <- sub(".bat", "", 
                list.files(paste(Sys.getenv("GISBASE"), "bin", sep="/"),
                pattern=".bat$"))
            if (nchar(Sys.getenv("GRASS_ADDON_BASE")) > 0) {
                t0 <- try(sub(".bat", "", 
                   list.files(paste(Sys.getenv("GRASS_ADDON_BASE"),
                       "bin", sep="/"), pattern=".bat$")), silent=TRUE)
                if (class(t0) != "try-error" && is.character(t0) &&
                   nchar(t0) > 0)
                   WN_bat <- c(WN_bat, t0)
            }
            assign("WN_bat", WN_bat, envir=.GRASS_CACHE)
        } 
        prep <- ""
        if ((get("SYS", envir=.GRASS_CACHE) == "WinNat") &&
            cmd %in% get("WN_bat", envir=.GRASS_CACHE))
            ext <- ".bat"

        cmd0 <- paste(paste(prep, cmd, ext, sep=""), "--interface-description")
        if (legacyExec) {
            tr <- try(system(cmd0, intern=TRUE))
	    if (class(tr) == "try-error") stop(paste(cmd, "not found"))
        } else {
            errFile <- tempfile()
            outFile <- tempfile()
            command <- paste(prep, cmd, ext, sep="")
            arguments <- "--interface-description"
            res <- system2(command=command, args=arguments, stdout=outFile,
                stderr=errFile)
            resErr <- readLines(errFile)
            if (res == 127L) {
               stop("The command\n", "   ", command, " ", arguments,
                   "\ncould not be run (", res,
                   "), and produced the error message:\n", "   ", resErr)
            } else if (res == 1L) {
                stop("The command\n", "   ", command, " ", arguments,
                    "\nproduced an error (", res,
                    ") during execution:\n", "   ", resErr)
            }
            tr <- readLines(outFile)
        }
        enc <- get("override_encoding", envir=.GRASS_CACHE)
        if (nchar(enc) > 0) {
          if (length(grep("UTF-8", tr[1])) > 0) {
            tr[1] <- sub("UTF-8", enc, tr[1])
          }        
        }
        tr <- try(xmlTreeParse(tr))
	if (inherits(tr, "try-error")) stop(paste(cmd, "not parsed"))
        tr1 <- xmlChildren(xmlRoot(tr))
        res <- vector(mode="list", length=9)
        names(res) <- c("cmd", "description", "keywords", "parameters", "flags",
            "pnames", "fnames", "ext", "prep")
        res$cmd <- cmd
        res$ext <- ext
        res$prep <- prep
        res$description <- xmlValue(tr1[[1]])
        res$keywords <- xmlValue(tr1[[2]])
        o0 <- names(tr1)
        pseq <- which(o0 == "parameter")
        np <- length(pseq)
        ps <- vector(mode="list", length=np)
        for (i in seq(along=pseq)) {
            obj <- tr1[[pseq[i]]]
            pv <- xmlAttrs(obj)
            pvCh <- xmlChildren(obj)
            pv <- c(pv, xmlValue(pvCh[["description"]]))
            default <- pvCh[["default"]]
            if (is.null(default)) {
                strdef <- as.character(NA)
            } else {
                strdef <- xmlValue(pvCh[["default"]])
                if (length(strdef) == 0) strdef <- ""
            }
            pv <- c(pv, strdef)
            kd <- pvCh[["keydesc"]]
            if (is.null(kd)) {
                nkd <- as.integer(NA)
                strkd <- as.character(NA)
            } else {
                kda <- xmlApply(kd, xmlValue)
                nkd <- length(kda)
                strkd <- paste(sapply(kda, c), collapse=",")
            }
            pv <- c(pv, nkd, strkd)
            names(pv) <- c(names(xmlAttrs(obj)), "desc", "default",
                "keydesc_count", "keydesc")
            ps[[i]] <- pv
        }
        res$parameters <- ps
        fseq <- which(o0 == "flag")
        nf <- length(fseq)
        fs <- vector(mode="list", length=nf)
        for (i in seq(along=fseq)) {
            obj <- tr1[[fseq[i]]]
            fv <- xmlAttrs(obj)
# find description, don't assume first place 101206
            nobj <- sapply(xmlChildren(obj), xmlName)
            di <- match("description", nobj)
            fv <- c(fv, xmlValue(xmlChildren(obj)[[di]]))
            suppr_req <- as.character("suppress_required" %in% nobj)
            fv <- c(fv, suppr_req)
            names(fv) <- c(names(xmlAttrs(obj)), "desc", "suppr_req")
            fs[[i]] <- fv
        }
        res$flags <- fs
        res$pnames <- sapply(res$parameters, function(x) x["name"])
        names(res$pnames) <- NULL
        res$fnames <- sapply(res$flags, function(x) x["name"])
        names(res$fnames) <- NULL
        class(res) <- "GRASS_interface_desc"
        cmdCACHE[[cmd]] <- res
        assign("cmdCACHE", cmdCACHE, envir=.GRASS_CACHE)
    } 
    res
} 

setXMLencoding <- function(enc) {
  if (!is.character(enc) || length(enc) > 1)
    stop("enc must be a character string")
  invisible(assign("override_encoding", enc, envir=.GRASS_CACHE))
}

getXMLencoding <- function() {
 get("override_encoding", envir=.GRASS_CACHE)
}
print.GRASS_interface_desc <- function(x, ...) {
    cat("Command:", x$cmd, "\n")
    if (nchar(x$ext) > 0) cat("Extension:", x$ext, "\n")
    if (nchar(x$prep) > 0) cat("Shell prefix:", x$prep, "\n")
    cat("Description:", x$description, "\n")
    cat("Keywords:", x$keywords, "\n")
    cat("Parameters:\n")
    for (i in x$parameters) {
        cat("  name: ", i["name"], ", type: ",
            i["type"], ", required: ", i["required"], ", multiple: ",
            i["multiple"], "\n", sep="")
        if (!is.na(i["default"])) cat("  default: ", i["default"],
            "\n", sep="")
        if (!is.na(i["keydesc"])) cat("  keydesc: ", i["keydesc"],
            ", keydesc_count: ", i["keydesc_count"], "\n", sep="")
        cat("[", i["desc"], "]\n", sep="")
    }
    cat("Flags:\n")
    for (i in x$flags) cat("  name: ", i["name"], " [",
        i["desc"], "] {", i["suppr_req"], "}\n", sep="")
    invisible(x)
}

doGRASS <- function(cmd, flags=NULL, ..., parameters=NULL, echoCmd=NULL, legacyExec=NULL) {
    defFlags <- get.defaultFlagsOption()
    if (!is.null(defFlags)) flags <- unique(c(flags, defFlags))
    if (all(c("quiet", "verbose") %in% flags)) {
        flags <- flags[flags != "quiet"]
    }
    if (!is.null(flags)) stopifnot(is.character(flags))
    if (!is.null(parameters)) stopifnot(is.list(parameters))
    if (is.null(echoCmd))
        echoCmd <- get("echoCmd", envir = .GRASS_CACHE)
    stopifnot(is.logical(echoCmd))
    if (is.null(legacyExec))
        legacyExec <- get.legacyExecOption()
    stopifnot(is.logical(legacyExec))
    if (get("SYS", envir=.GRASS_CACHE) != "unix" && !legacyExec) {
        warning("legacyExec set TRUE on non-unix platform")
        legacyExec <- TRUE
    }

#    G6 <- get("GV", envir=.GRASS_CACHE) < "GRASS 7"

    dlist <- list(...)
    if (!is.null(parameters) && (length(dlist) > 0))
        stop(paste("Use either GRASS parameters as R arguments,",
        "or as a parameter argument list object, but not both", sep="\n"))
    if (is.null(parameters) && (length(dlist) > 0)) parameters <- dlist
    pcmd <- parseGRASS(cmd, legacyExec=legacyExec)
    cmd <- paste(pcmd$prep, cmd, pcmd$ext, sep="")
    res <- cmd
    suppress_required <- FALSE
    if (!is.null(flags)) {
        fm <- match(flags, pcmd$fnames)
        if (any(is.na(fm))) {
            print(pcmd)
            stop(paste("Invalid flag value:", flags[is.na(fm)]))
        }
        suppr_req <- as.logical(sapply(pcmd$flags, "[", "suppr_req"))
        if (!suppress_required && any(suppr_req[fm])) suppress_required <- TRUE
        dble_dash <- c("verbose", "overwrite", "quiet")
        dbm <- na.omit(match(dble_dash, flags))
        if (length(dbm) > 0) flags[dbm] <- paste("-", flags[dbm], sep="")
        res <- paste(res, paste("-", flags, collapse=" ", sep=""))
    }
    pt <- do.call("rbind", pcmd$parameters)
    if (!is.null(pt)) {
# g.version no parameters exception caught by Rainer M Krug 090923
      req <- pt[pt[, "required"] != "no", "name"]
# patch for multiple values Patrick Caldon 090524
#      mult <- pt[pt[, "multiple"] != "no", "name"]
# patch to accept no or multiple keydesc_count 090902
      mults <- pt[, "multiple"] != "no" | (!is.na(pt[, "keydesc_count"]) &
        pt[, "keydesc_count"] > 1)
      mult <- pt[mults, "name"]
      parameters <- insert_required(pcmd=pcmd, parameters=parameters,
          pt=pt, req=req, suppress_required=suppress_required)
      if (!suppress_required && length(req) > 0 && is.null(parameters)) {
        print(pcmd)
        stop("No parameters given where some are required with defaults declared")
      }
      if (!is.null(parameters)) {
        parnms <- names(parameters)
        pm <- match(parnms, pcmd$pnames)
        if (any(is.na(pm))) {
            print(pcmd)
            stop(paste("Invalid parameter name:", parnms[is.na(pm)]))
        }
        if (!suppress_required && length(req) > 0) {
            pmr <- match(req, parnms)
            if (any(is.na(pmr))) {
                print(pcmd)
                stop(paste("Missing required parameter:", req[is.na(pmr)]))
            }
        }
        pmv <- pt[pm, "type"]
# patch for multiple values Patrick Caldon 090524
        pmmult <- match(mult, parnms)
        for (i in seq(along=parameters)) {
# patch for multiple values Patrick Caldon 090524
             if (length(parameters[[i]]) > 1 && !(i %in% pmmult))
                stop(paste("Parameter <", names(parameters)[i],
                    "> has multiple values", sep=""))
            if (pmv[i] == "string") {
                if (!is.character(parameters[[i]]))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> does not have string value", sep=""))
# added any() 140516
                if (any(is.na(parameters[[i]])))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> is NA", sep=""))
# string space protection 091108 Martin Mainzer
                Space <- length(grep(" ", parameters[[i]])) > 0
# Rainer Krug 110128
                Paran <- length(grep("\\(", parameters[[i]])) > 0 ||
                    length(grep(")", parameters[[i]])) > 0
                if (Space || Paran) {
# extra protection against existing escaping of quotes 100422
                  if (length(grep("\"", parameters[[i]])) == 0) {
                    parameters[[i]] <- paste("\"", parameters[[i]], "\"",
                        sep="")
                  }
                }
            } else if (pmv[i] == "float") {
                if (!is.numeric(parameters[[i]]))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> does not have numeric value", sep=""))
                if (any(!is.finite(parameters[[i]])))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> is not finite", sep=""))
            } else if (pmv[i] == "integer") {
                if (!is.numeric(parameters[[i]]))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> does not have numeric value", sep=""))
                if (any(!is.finite(parameters[[i]])))
                    stop(paste("Parameter <", names(parameters)[i],
                    "> is not finite", sep=""))
                if (!is.integer(parameters[[i]])) {
                    opi <- parameters[[i]]
                    if (all(as.integer(opi) == opi)) {
                        parameters[[i]] <- as.integer(opi)
                    } else {
                        stop(paste("Parameter <", names(parameters)[i],
                        "> is not integer", sep=""))
                    }
                }
            } else warning("unknown parameter type")
# patch for multiple values Patrick Caldon 090524
            param <- paste(parameters[[i]], collapse=",")
            res <- paste(res, paste(names(parameters)[i], param,
                sep="="))
        }
      }
    }
    if ((!is.null(pt) && is.null(parameters)) && is.null(flags))
        if (get.stop_on_no_flags_parasOption())
            stop("No flags or parameters provided")
        else warning("No flags or parameters provided")
    if (echoCmd) cat("GRASS command:", res, "\n")
    attr(res, "cmd") <- cmd
    res
}

insert_required <- function(pcmd, parameters, pt, req, suppress_required) {
    defs <- pt[match(req, pt[, "name"]), "default"]
    nadefs <- which(is.na(defs))
    nms <- pt[match(req, pt[, "name"]), "name"]
    nadefnms <- nms[nadefs]
    pnms <- names(parameters)
    nadefnms1 <- nadefnms[is.na(match(nadefnms, pnms))]
    if (!suppress_required && length(nadefnms1) > 0) {
        print(pcmd)
        stop(paste("required parameters with no defaults missing:",
        paste(nadefnms1, collapse=" ")))
    }
    types <- pt[match(req, pt[, "name"]), "type"]
    defnms <- nms[!is.na(defs)]
    deftypes <- types[!is.na(defs)]
    defdefs <- defs[!is.na(defs)]
    if (is.null(parameters)) parameters <- list()
    for (i in seq(along=defnms)) {
        if (!(defnms[i] %in% pnms)) {
            parameters[[defnms[i]]] <- switch(deftypes[i],
                integer = as.integer(defdefs[i]),
                float = as.numeric(defdefs[i]),
                string = defdefs[i])
        }
    }
    if (length(parameters) == 0) parameters <- NULL
    parameters
}

execGRASS <- function(cmd, flags=NULL, ..., parameters=NULL, intern=NULL,
    ignore.stderr=NULL, Sys_ignore.stdout=FALSE, Sys_wait=TRUE,
    Sys_input=NULL, Sys_show.output.on.console=TRUE, Sys_minimized=FALSE,
    Sys_invisible=TRUE, echoCmd=NULL, redirect=FALSE, legacyExec=NULL) {
    if (is.null(ignore.stderr))
        ignore.stderr <- get.ignore.stderrOption()
    stopifnot(is.logical(ignore.stderr))
    if (is.null(intern))
        intern <- get.useInternOption()
    stopifnot(is.logical(intern))
    if (is.null(legacyExec))
        legacyExec <- get.legacyExecOption()
    stopifnot(is.logical(legacyExec))
    if (get("SYS", envir=.GRASS_CACHE) != "unix" && !legacyExec) {
        warning("legacyExec set TRUE on non-unix platform")
       legacyExec <- TRUE
    }

    syscmd <- doGRASS(cmd, flags=flags, ..., parameters=parameters,
        echoCmd=echoCmd, legacyExec=legacyExec)
    if (legacyExec) {
        if (redirect) {
            syscmd <- paste(syscmd, "2>&1")
            intern=TRUE
        }
        if (get("SYS", envir=.GRASS_CACHE) == "unix") {
            res <- system(syscmd, intern=intern, ignore.stderr=ignore.stderr,
                ignore.stdout=Sys_ignore.stdout, wait=Sys_wait, input=Sys_input)
        } else {
            res <- system(syscmd, intern=intern, ignore.stderr=ignore.stderr,
                ignore.stdout=Sys_ignore.stdout, wait=Sys_wait, input=Sys_input,
                show.output.on.console=Sys_show.output.on.console, 
                minimized=Sys_minimized, invisible=Sys_invisible)
        }
        if (intern) return(res)
    } else {
        command <- attr(syscmd, "cmd")
        arguments <- substring(syscmd, (nchar(command)+2), nchar(syscmd))

        errFile <- tempfile(fileext=".err")
        outFile <- tempfile(fileext=".out")

        res <- system2(command, arguments, stderr = errFile,
            stdout = outFile, wait=Sys_wait, input=Sys_input)

        resErr <- readLines(errFile)
        if (res == 127L) {
             stop("The command:\n", command, " ", arguments,
                 "\ncould not be run (", res,
                 "), and produced the error message:\n",
                 paste(resErr, collapse="\n"))
        } else if (res == 1L) {
            stop("The command:\n", command, " ", arguments,
                "\nproduced an error (", res,
                ") during execution:\n", paste(resErr, collapse="\n"))
        } else if (res == 0L & !ignore.stderr) {
            if (length(grep("ERROR:", resErr)) > 0) {
                stop("The command:\n", command, " ", arguments,
                    "\nproduced an error (", res,
                    ") during execution:\n",
                    paste(resErr, collapse="\n"))
            } else if (length(grep("WARNING:", resErr)) > 0) {
                warning("The command:\n", command, " ", arguments,
                    "\nproduced at least one warning during execution:\n",
                    paste(resErr, collapse="\n"))
            }
        }

        resOut <- readLines(outFile)
        if (intern) return(resOut)

        if (length(resOut) > 0) cat(resOut, sep="\n")
        if (length(resErr) > 0) cat(resErr, sep="\n")
       
        attr(res, "resOut") <- resOut
        attr(res, "resErr") <- resErr
    }
    invisible(res)
}


