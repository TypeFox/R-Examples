

### assumes that all chars are single-byte
write.dbf.SSN <- function(dataframe, file, factor2char = TRUE, max_nchar = 254)
{

    suppressWarnings(warning("write.dbf.SSN"))
### need to check precision
    allowed_classes <- c("logical", "integer", "numeric", "character",
                         "factor", "Date")

    if (!is.data.frame(dataframe)) dataframe <- as.data.frame(dataframe)
    if (any(sapply(dataframe, function(x) !is.null(dim(x)))))
        stop("cannot handle matrix/array columns")
    cl <- sapply(dataframe, function(x) class(x[1L]))
    asis <- cl == "AsIs"
    cl[asis & sapply(dataframe, mode) == "character"] <- "character"
    if(length(cl0 <- setdiff(cl, allowed_classes)))
        stop(sprintf(ngettext(length(cl0),
                              "data frame contains columns of unsupported class %s",
                              "data frame contains columns of unsupported classes %s"),

                     paste(dQuote(cl0), collapse = ",")), domain = NA)
    m <- ncol(dataframe)
    DataTypes <- c(logical="L", integer="N", numeric="F", character="C",
                   factor=if(factor2char) "C" else "N", Date="D")[cl]
    for(i in seq_len(m)) {
        x <- dataframe[[i]]
        if(is.factor(x))
            dataframe[[i]] <-
                if(factor2char) as.character(x) else as.integer(x)
        else if (inherits(x, "Date"))
            dataframe[[i]] <- format(x, "%Y%m%d")
    }
    precision <- integer(m)
    scale <- integer(m)
    dfnames <- names(dataframe)
    for (i in seq_len(m)) {
        nlen <- nchar(dfnames[i], "b")
        x <- dataframe[, i]
        if (is.logical(x)) {
            precision[i] <- 1L
            scale[i] <- 0L
        } else if (is.integer(x)) {
            rx <- range(x, na.rm = TRUE)
            rx[!is.finite(rx)] <- 0 # added RSB 2005-04-17
	    if (any(rx == 0)) rx <- rx + 1 # added RSB 2005-03-10
            mrx <- as.integer(max(ceiling(log10(abs(rx))))+3L)
            precision[i] <- min(max(nlen, mrx), 19L)
            scale[i] <- 0L
        } else if (is.double(x)) {
            precision[i] <- 19L
            rx <- range(x, na.rm = TRUE)
            rx[!is.finite(rx)] <- 0 # added RSB 2005-04-17
            mrx <- max(ceiling(log10(abs(rx))))
            scale[i] <- min(precision[i] - ifelse(mrx > 0L, mrx+3L, 3L), 15L)
                    # modified RSB 2005-03-10 and 2005-04-17
        } else if (is.character(x)) {
            if(length(x) == sum(is.na(x))){
                mf = max_nchar
            } else {
                mf <- max(nchar(x[!is.na(x)], "b"))
            }
            p <- max(nlen, mf)
            if(p > max_nchar)
                warning(gettext("character column %d will be truncated to %d bytes", i, max_nchar), domain = NA)
            precision[i] <- min(p, max_nchar)
            scale[i] <- 0L
        } else stop("unknown column type in data frame")
    }
    if (any(is.na(precision))) stop("NA in precision") # added RSB 2005-04-17
    if (any(is.na(scale))) stop("NA in scale") # added RSB 2005-04-17
#################################################################################
    invisible( .Call("DoWritedbfSSN", as.character(file),
                     dataframe, as.integer(precision), as.integer(scale),
                     as.character(DataTypes)))

################################################################################
}
