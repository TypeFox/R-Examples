checkType <- function(dsList, id, dat=NULL)
{
    stopifnot(is.data.frame(dsList))
    stopifnot(length(id) == 1)
    if (is.character(id)) {
        ind <- which(dsList$identification == id)
        stopifnot(length(ind) == 1)
        id <- ind
    }
    cat(dsList$identification[id], ": Check of data types ")
    if (is.null(dat)) {
        dat <- dsRead(dsList, id, keepContents=TRUE)
    }
    OK <- TRUE
    expected <- split.comma(dsList$originalColsType[id])
    obtainedClass <- rep("", times=ncol(dat))
    lengthUnique <- rep(NA, times=ncol(dat))
    columnOK <- rep(FALSE, times=ncol(dat))
    if (length(expected) == ncol(dat)) {
        for (i in seq.int(length.out=length(expected))) {
            x <- dat[[i]]
            x <- x[!is.na(x)]
            if (is.double(x) && all(x == round(x)))
                x <- as.integer(x)
            x.class <- class(x)
            stopifnot(length(x.class) == 1)
            obtainedClass[i] <- x.class
            lengthUnique[i] <- length(unique(x))
            if (length(x) == 0) {
                obtainedClass[i] <- "null"
                columnOK[i] <- expected[i] == "1"
            } else if (obtainedClass[i] == "numeric") {
                columnOK[i] <- expected[i] == "n"
            } else if (obtainedClass[i] == "integer") {
                columnOK[i] <- expected[i] == "n" || as.integer(expected[i]) >= length(unique(x))
            } else if (obtainedClass[i] == "character") {
                columnOK[i] <- expected[i] != "n" && as.integer(expected[i]) >= length(unique(x))
            }
        }
        if (any(!columnOK)) {
            cat("\n")
            cat("\n")
            cat("The types of columns do not match.\n")
            out <- data.frame(expected, obtainedClass, lengthUnique)
            print(out[!columnOK, ])
            OK <- FALSE
        }
    } else {
        cat("\n")
        cat("\n")
        cat("The number of columns does not match.\n")
        cat("stored = ", length(expected), ", obtained = ", ncol(dat), "\n", sep="")
        OK <- FALSE
    }
    if (OK) {
        cat("OK\n")
    } else {
        cat("Check of data types FAILED\n")
        cat("\n")
    }
    invisible(OK)
}

checkConsistency <- function(dsList, outputInd=FALSE)
{
    stopifnot(is.data.frame(dsList))
    OK <- TRUE
    tab <- table(dsList$identification)
    if (any(tab > 1))
    {
        cat("The following databases have several occurrences in the list\n")
        print(cbind(occurrences=tab[tab > 1]))
        cat("\n")
        OK <- FALSE
    }
    if (ncol(dsList) != nrow(fieldInfo)) {
        cat("dsList has", ncol(dsList), "columns instead of", nrow(fieldInfo), "\n")
        OK <- FALSE
    }
    if (OK) {
        obtained <- rep("", times=ncol(dsList))
        for (i in seq.int(length.out=ncol(dsList))) {
            obtained[i] <- class(dsList[, fieldInfo[i, 1]])[1]
        }
        if (any(obtained != fieldInfo[, 2])) {
            cat("incorrect class of columns\n")
            print(cbind(column=fieldInfo[, 1], expected=fieldInfo[, 2], obtained))
            OK <- FALSE
        }
    }
    ind <- c()
    for (i in 1:nrow(dsList))
    {
        orig.num <- dsList$originalColsNumber[i]
        calc.orig.type.num <- length(split.comma(dsList$originalColsType[i]))
        calc.orig.names.num <- length(split.comma(dsList$originalColsNames[i]))
        delete <- as.integer(split.comma(dsList$delete[i]))

        ok1 <- !(dsList$responsePos[i] %in% delete)
        ok2 <- orig.num == calc.orig.type.num
        ok3 <- calc.orig.names.num == 0 || orig.num == calc.orig.names.num

        if (!all(ok1, ok2, ok3))
        {
            cat("identification: ", dsList$identification[i], "\n")
            cat("\n")
            ind <- c(ind, i)
            OK <- FALSE
        }

        if (!ok1)
        {
            cat("delete                    :", delete, "\n")
            cat("responsePos               :", dsList$responsePos[i], "\n")
            cat("\n")
            OK <- FALSE
        }

        if (!ok2)
        {
            cat("originalColsNumber        :", orig.num, "\n")
            cat("length(originalColsType)  :", calc.orig.type.num, "\n")
            cat("\n")
            OK <- FALSE
        }

        if (!ok3)
        {
            cat("originalColsNumber        :", orig.num, "\n")
            cat("length(originalColsNames) :", calc.orig.names.num, "\n")
            cat("\n")
            OK <- FALSE
        }
    }
    if (!OK) {
        cat("Check of dsList consistency FAILED\n")
    }
    if (outputInd)
        return(ind)
    else
        return(invisible(OK))
}

