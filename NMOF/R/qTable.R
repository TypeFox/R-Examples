qTable  <- function(X, xmin = NULL, xmax = NULL, labels = NULL,
                    at = NULL, unitlength = "5cm", linethickness = NULL,
                    cnames = colnames(X), circlesize = 0.01,
                    xoffset = 0, yoffset = 0, dec = 2L, filename = NULL,
                    funs = list(median=median, min=min, max=max),
                    tabular.format, skip = TRUE) {

    if (missing(tabular.format))
        tabular.format <- paste("r",
                                paste(rep("r", length(funs)), collapse = ""),
                                "r", sep ="")
    
    X <- as.matrix(X)
    if (nrow(X) < 10L)
        warning(sQuote("X"), " has less than 10 rows")
    if (!is.null(at) && is.null(labels))
        stop(sQuote("labels"), " must be provided")
    if (!is.null(labels) && is.null(at))
        stop(sQuote("at"), " must be provided")

    ## compute quantiles
    A <- apply(X, 2L, quantile, c(0.25, 0.5, 0.75))
    iqr <- abs(A[3L, ] - A[1L, ])

    ## compute whiskers
    B <- array(0, dim = c(5L, dim(A)[2L]))
    B[1L, ] <- pmax(A[1L,] - 1.5 * iqr, apply(X, 2L, min))
    B[2L:4L, ] <- A
    B[5L, ] <- pmin(A[3L,] + 1.5 * iqr, apply(X, 2L, max))

    ## ranges of plot
    if (is.null(xmin))
        xmin <- min(B)
    if (is.null(xmax))
        xmax <- max(B)
    joli <- pretty(c(xmin, xmax), n = 3L)
    a <- min(joli)
    b <- max(joli)
    if (is.null(labels))
        labels <- joli
    if (is.null(at))
        at <- joli

    ## ranges of picture (map to [0,1])
    B  <-  (B - a)/(b - a)
    at <- (at - a)/(b - a)

    ## helper functions
    fff <- function(x)
        formatC(x, digits = dec, format = "f")
    ff <- function(x)
        formatC(x, digits = 5L,  format = "f")

    STR <- paste("\n & ", paste(names(funs), collapse=" & "), "\\\\")
    for (cc in seq_len(dim(X)[2L])) {

        if (length(funs)>0L) {
            str0 <- paste(cnames[cc], " & ", 
                          paste(fff(sapply(funs, do.call, list(X[ ,cc]))),
                                collapse = " & "),
                          sep = "")
        } else {
            str0 <- cnames[cc]
        }
        
        str1 <- paste("& \\begin{picture}(1,0)(",
                      xoffset, ",", yoffset, ")", sep = "")

        str2 <- paste("\\put(", ff(B[1,cc]), ",0){\\line(1,0){",
                      ff(B[2,cc] - B[1,cc]), "}}", sep = "")

        str3 <- paste("\\put(", ff(B[3,cc]),
                      ",0){\\circle*{", circlesize, "}}", sep = "")

        str4 <- paste("\\put(", ff(B[4,cc]), ",0){\\line(1,0){",
                      ff(B[5,cc] - B[4,cc]), "}}", sep = "")

        str5 <- "\\end{picture}\\\\"

        temp <- paste(str0, str1, str2, str3, str4, str5, sep = "")
        STR  <- rbind(STR, "\n", temp)
    }

    ## ... axis line
    amp <- paste(rep("&", length(funs)+1L), collapse = "")
    strScale <- paste(amp,"\\begin{picture}(1,0)(",
                      xoffset, ",", yoffset,
                      ")\\put(0,0){\\line(1,0){1}}", sep = "")

    for (i in seq(along.with = labels)) {
        strScale <- paste(strScale, "\\put(", ff(at[i]),
                          ",0) {\\line(0,-1){0.01}}\n", sep = "")

        strScale <- paste(strScale, "\\put(", ff(at[i]), ",-0.1){",
                          labels[i], "}\n", sep = "")
    }

    ## end line
    strScale <- paste(strScale,"\\end{picture}\\\\", sep = "")

    STR <- rbind(paste("\\begin{tabular}{", tabular.format, "}", sep =""),
                 STR,
                 "\n", strScale,
                 ifelse(skip,"\n\\\\\\end{tabular}","\n\\end{tabular}"))
    
    ## add line thickness
    if (!is.null(linethickness))
        STR <- rbind(paste("\\linethickness{",
                           linethickness, "}\n", sep = ""), STR)

    ## add unit length
    temp <- paste("\\setlength{\\unitlength}{",
                  unitlength, "}\n", sep = "")

    STR  <- rbind(temp, STR)

    ## add braces
    STR <- rbind("{\n", STR, "\n}\n")
    if (!is.null(filename)) {
        cat(STR, "\n", sep = "", file = filename)
    }
    STR
    
}

