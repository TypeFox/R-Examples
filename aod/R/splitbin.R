splitbin <- function(formula, data, id = "idbin"){
    #------------------------------------------------------------------------------------------------------------------------
    splitbin.w <- function(formula, data, id){
    f <- formula
    listvar <- all.vars(f)
    # build the data set
    dfr <- model.frame(formula = f, data = data)
    rownam <- rownames(data)
    resp <- model.response(dfr)
    #n <- rowSums(resp)
    #y <- resp[ , 1]
    n <- round(rowSums(resp))
    y <- round(resp[ , 1])
    # check the response
    if(any(n == 0)){
        err <- rownam[n == 0]
        rownam <- rownam[n > 0]
        y <- y[n > 0]
        dfr <- dfr[n > 0, ]
        n <- n[n > 0]   ### beware: condition must be applied on "y" and "dfr" BEFORE "n"
        warning("Lines with row names ",
            paste(sapply(err, function(x) dQuote(x)), collapse = ", "),
            " were discarded (0 observation).\n")
        }
    if(any(y > n))
        stop("Some values of ", dQuote(listvar[1]), " were > ", dQuote(listvar[2]), ".")

    nc <- length(listvar) - 2
    if(nc > 0)
        dfr <- dfr[ , -1, drop = FALSE]
    List <- vector(mode = "list", length = length(n))
    for(i in seq(length(n))){
        dat <- data.frame(id = rep(rownam[i], each = n[i]),
            resp = rep(c(0, 1), times = c(n[i] - y[i], y[i])))
        names(dat) <- c(id, listvar[1])
        List[[i]] <- if(nc == 0) dat else merge(dat, dfr[i, , drop = FALSE], all = TRUE)
        }

    tab <- do.call("rbind", List)
    rownames(tab) <- seq(nrow(tab))
    tab
    }
    #------------------------------------------------------------------------------------------------------------------------
    dfname <- as.character(substitute(data))
    if(!is.data.frame(data))
        stop(paste("\n\nThe object ", dfname, " is not data frame.\n\n", sep = ""))
    CALL <- match.call()
    # in case the formula was provided as an object
    CALL[[2]] <- formula(deparse(formula))
    # convert any variable of mode character into a factor
    data <- as.data.frame(lapply(data, function(x) if(mode(x) == "character") factor(x) else x))
    # Grouped data of the form cbind(success, failure)
    if(substring(deparse(formula)[1], 1, 5) == "cbind")
        tab <- splitbin.w(formula = formula, data = data, id = id)
    # Grouped data with weights
    else {
        listvar <- all.vars(formula)
        names(data)[match(listvar[1], names(data))] <- "nxxx"
        data$nbsuccessxxx <- data$nxxx
        data <- data[data$nxxx > 0, ] # remove row(s) of null weight
        newf <- paste("cbind(nbsuccessxxx, nxxx - nbsuccessxxx) ~ ", paste(listvar[-1], collapse = " + "), sep = "")
        newf <- formula(newf)
        tab <- splitbin.w(formula = newf, data = data, id = id)
        tab <- tab[, -match("nbsuccessxxx", names(tab))]
        }
    ## Output
    list(CALL = CALL, tab = tab)
    }
