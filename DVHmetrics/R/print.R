## print method for one DVH
print.DVHs <- function(x, ...) {
    dots <- list(...)

    structVol <- if(!is.na(x$structVol)) {
        paste(signif(x$structVol, 2), x$volumeUnit)
    } else {
        "? CC"
    }

    if(("noID" %in% names(dots)) && (dots$noID == TRUE)) {
        cat("DVH: Structure '", x$structure,
            "' (", structVol, "),",
            " Dose: ", paste(signif(range(x$dvh[ , "dose"], na.rm=TRUE), 2), collapse="-"),
            x$doseUnit, "\n", sep="")
    } else {
        cat("DVH: Patient '", x$patName,
            "' (ID ", x$patID,
            "), structure '", x$structure,
            "' (", structVol, "),",
            " Dose: ", paste(signif(range(x$dvh[ , "dose"], na.rm=TRUE), 2), collapse="-"),
            x$doseUnit, "\n", sep="")
    }

    return(invisible(NULL))
}

## dvhInfo method for DVH patient/structure list
print.DVHLst <- function(x, ...) {
    dots <- list(...)
    doseRx <- if(!is.na(x[[1]]$doseRx)) {
        c(", prescription dose ", x[[1]]$doseRx, x[[1]]$doseUnit)
    } else {
        NULL
    }

    if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == TRUE)) {
        cat("DVH list:\nPatient '",  x[[1]]$patName,
            "' (ID ", x[[1]]$patID, doseRx, ") with ", length(x), sep="")
    } else if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == FALSE)) {
        cat("DVH list:\nStructure '", x[[1]]$structure, "' with ", length(x), sep="")
    } else {
        cat("DVH list with", length(x))
    }

    if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == TRUE)) {
        cat(" Structures:\n")
    } else if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == FALSE)) {
        cat(" Patient IDs:\n")
    } else {
        cat(" DVHs:\n")
    }

    if(("verbose" %in% names(dots)) && (dots$verbose == TRUE)) {
        if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == TRUE)) {
            Map(print, x, noID=TRUE, ...)
        } else {
            Map(print, x, noID=FALSE, ...)
        }
    } else {
        cat(names(x), sep=", ", fill=TRUE)
    }

    cat("\n")
    return(invisible(NULL))
}

## dvhInfo method for DVH patient list of lists
print.DVHLstLst <- function(x, ...) {
    if(!is.null(attributes(x)$comment)) {
        cat(attributes(x)$comment, "\n\n")
    }

    if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == TRUE)) {
        cat("DVH list of", length(x), "lists - 1 for each patient:\n\n")
    } else if(!is.null(attributes(x)$byPat) && (attributes(x)$byPat == FALSE)) {
        cat("DVH list of", length(x), "lists - 1 for each structure:\n\n")
    } else {
        cat("DVH list of", length(x), "lists:\n\n")
    }

    Map(print, x, ...)

    return(invisible(NULL))
}
