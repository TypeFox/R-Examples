`XYplot` <- function(x, y, data, relation = "nec", size = 15, mguides = TRUE,
                     jitter = FALSE, ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    other.args <- list(...)
    funargs <- unlist(lapply(match.call(), deparse)[-1])
    
    if (missing(x)) {
        cat("\n")
        stop(simpleError("Argument x is mandatory.\n\n"))
    }
    
    negated <- logical(2)
    
    xname <- yname <- ""
    
    testit <- capture.output(tryCatch(eval(x), error = function(e) e))
    if (length(testit) == 1 & is.character(testit)) {
        if (grepl("Error", testit)) {
            x <- as.vector(funargs["x"])
        }
    }
    
    if (!is.character(x)) {
        
        testit <- capture.output(tryCatch(eval(x), error = function(e) e))
        if (length(testit) == 1 & is.character(testit)) {
            if (grepl("Error", testit)) {
                x <- as.vector(deparse(funargs["x"]))
            }
            else if (grepl("~", testit)) {
                negated[1] <- TRUE
                # try check if it's an object
                if (eval.parent(parse(text=paste0("\"", gsub("~", "", testit), "\" %in% ls()")), n = 1)) {
                    x <- 1 - eval.parent(parse(text=paste("get(\"", gsub("~", "", testit), "\")", sep="")), n = 1)
                }
                else {
                    x <- testit
                }
            }
        }
    }
    else {
        if (x == tolower(x) & x != toupper(x)) {
            # check if it's an object negated with lower case letters
            if (eval.parent(parse(text=paste0("\"", toupper(x), "\" %in% ls()")), n = 1)) {
                conds <- toupper(x)
                x <- 1 - eval.parent(parse(text=paste("get(\"", toupper(x), "\")", sep="")), n = 1)
                negated[1] <- TRUE
            }
        }
    }
    
    
    if (!is.data.frame(x) & !is.matrix(x) & !missing(y)) {
        testit <- capture.output(tryCatch(eval(y), error = function(e) e))
        if (length(testit) == 1 & is.character(testit)) {
            if (grepl("Error", testit)) {
                y <- as.vector(funargs["y"])
            }
        }
        
        if (!is.character(y)) {
            testit <- capture.output(tryCatch(eval(y), error = function(e) e))
            if (length(testit) == 1 & is.character(testit)) {
                if (grepl("Error", testit)) {
                    y <- deparse(funargs["y"])
                }
                else if (grepl("~", testit)) {
                    negated[2] <- TRUE
                    # try check if it's an object
                    if (eval.parent(parse(text=paste0("\"", gsub("~", "", testit), "\" %in% ls()")), n = 1)) {
                        y <- 1 - eval.parent(parse(text=paste("get(\"", gsub("~", "", testit), "\")", sep="")), n = 1)
                    }
                    else {
                        y <- testit
                    }
                }
            }
        }
        else {
            if (y == tolower(y) & y != toupper(y)) {
                # check if it's an object negated with lower case letters
                if (eval.parent(parse(text=paste0("\"", toupper(y), "\" %in% ls()")), n = 1)) {
                    conds <- toupper(y)
                    y <- 1 - eval.parent(parse(text=paste("get(\"", toupper(y), "\")", sep="")), n = 1)
                    negated[2] <- TRUE
                }
            }
        }
    }
    
    
    if (is.character(x)) {
        if (length(x) == 1) {
            x <- QCA::splitstr(x)
        }
        
        if (length(x) == 1) {
            if (missing(y)) {
                cat("\n")
                stop(simpleError("The outcome's name is missing.\n\n"))
            }
            else if (!is.character(y)) {
                cat("\n")
                stop(simpleError("x and y should be column names from the data.\n\n"))
            }
        }
        else {
            if (!missing(y)) {
                if (is.data.frame(y)) {
                    data <- y
                }
            }
            
            y <- x[2]
            x <- x[1]
        }
        
        if (missing(data)) {
            cat("\n")
            stop(simpleError("Data is missing.\n\n"))
        }
        else {
            verify.qca(data)
        }
        
        
        x <- gsub("[[:space:]]", "", x)
        y <- gsub("[[:space:]]", "", y)
        negated <- grepl("1-|~", c(x, y))
        x <- gsub("1-|~", "", x)           
        y <- gsub("1-|~", "", y)
        
        if (!all(c(x, y) %in% colnames(data))) {
            cat("\n")
            stop(simpleError("x and y should be column names from the data.\n\n"))
        }
        
        xname <- x
        yname <- y
        x <- data[, x]
        y <- data[, y]
        if (negated[1]) {
            x <- 1 - x
        }
        if (negated[2]) {
            y <- 1 - y
        }
        
    }
    else if (is.data.frame(x) | is.matrix(x)) {
        verify.qca(as.data.frame(x))
        if (ncol(x) < 2) {
            cat("\n")
            stop(simpleError("At least two columns are needed.\n\n"))
        }
        xname <- colnames(x)[1]
        yname <- colnames(x)[2]
        y <- x[, 2]
        x <- x[, 1]
    }
    else if (!missing(y)){
        
        if (length(x) > 1 & is.numeric(x)) { # X is a numeric vector
            negated[1] <- grepl("1-|~", gsub("[[:space:]]", "", funargs[1]))
            xname <- "X"
            tc <- capture.output(tryCatch(QCA::getName(funargs[1]), error = function(e) e, warning = function(w) w))
            if (!grepl("simpleError", tc)) {
                xname <- gsub("~", "", QCA::getName(funargs[1]))
            }
        }
        
        if (length(y) > 1 & is.numeric(y)) { # Y is a numeric vector
            negated[2] <- grepl("1-|~", gsub("[[:space:]]", "", funargs[2]))
            yname <- "Y"
            tc <- capture.output(tryCatch(QCA::getName(funargs[2]), error = function(e) e, warning = function(w) w))
            if (!grepl("simpleError", tc)) {
                yname <- QCA::getName(funargs[2])
            }
        }
        
        if (length(y) == 1 & is.character(y)) {
            if (missing(data)) {
                cat("\n")
                stop(simpleError("Data is missing.\n\n"))
            }
            else {
                verify.qca(data)
            }
            
            y <- gsub("[[:space:]]", "", y)
            negated[2] <- grepl("1-|~", y)
            y <- gsub("1-|~", "", y)
            
            if (!is.element(y, colnames(data))) {
                cat("\n")
                stop(simpleError("x and y should be column names from the data.\n\n"))
            }
            else {
                yname <- as.vector(y)
                y <- data[, y]
                if (negated[2]) {
                    y <- 1 - y
                }
            }
        }
    }
    else {
        cat("\n")
        stop(simpleError("Either a dataframe with two columns or two vectors are needed.\n\n"))
    }
    
    
    # x and y should be bound between 0 and 1
    if (any(x > 1) | any(y > 1)) {
        cat("\n")
        stop(simpleError("Values should be bound between 0 and 1.\n\n"))
    }
    
    xcopy <- x
    ycopy <- y
    
    if (dev.cur() == 1) {
        dev.new(width = (size + 1)/2.54, height = (size + 1)/2.54)
    }
    
    jitfactor <- 0.01
    jitamount <- 0.01
    cexpoints <- 0.8
    cexaxis <- 0.8
    pch <- 21
    bgpoints <- "#707070"
    
    if (length(testarg <- which(names(other.args) == "factor")) > 0) {
        jitfactor <- other.args$factor
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "amount")) > 0) {
        jitamount <- other.args$amount
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "cex")) > 0) {
        cexpoints <- other.args$cex
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "bg")) > 0) {
        bgpoints <- other.args$bg
        other.args <- other.args[-testarg]
    }
    
    if (jitter) {
        x <- jitter(x, jitfactor, jitamount)
        y <- jitter(y, jitfactor, jitamount)
    }
    
    toplot <- list(as.name("plot"), x, y)
    
    xlabel <- paste0(ifelse(negated[1], "~", ""), xname)
    ylabel <- paste0(ifelse(negated[2], "~", ""), yname)
    
    if (length(testarg <- which(names(other.args) == "xlab")) > 0) {
        xlabel <- other.args$xlab
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "ylab")) > 0) {
        ylabel <- other.args$ylab
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "cex.axis")) > 0) {
        cexaxis <- other.args$cex.axis
        other.args <- other.args[-testarg]
    }
    
    if (length(testarg <- which(names(other.args) == "pch")) > 0) {
        pch <- other.args$pch
        other.args <- other.args[-testarg]
    }
    
    toplot$type <- "n"
    toplot$xlim <- c(0, 1)
    toplot$ylim <- c(0, 1)
    toplot$axes <- FALSE
    toplot$xlab <- ""
    toplot$ylab <- ""
    
    if (length(other.args) > 0) {
        toplot <- c(toplot, other.args)
    }
    
    if (dev.cur() == 1) {
        dev.new(width = (size + 1)/2.54, height = (size + 1)/2.54)
    }
    
    par(mar = c(3, 3.5, 2.5, 0.5))
    suppressWarnings(eval(as.call(toplot)))
    
    mtext(xlabel, 1, 2, font = 2, cex = cexaxis + 0.1)
    mtext(ylabel, 2, 2.5, font = 2, cex = cexaxis + 0.1, las = 3)
    
    box()
    axis(1, at = seq(0, 1, .1), labels = rep("", 11), cex.axis = cexaxis, tck = -.015)
	axis(2, at = seq(0, 1, .1), labels = rep("", 11), cex.axis = cexaxis, tck = -.015, las = 2)
	
	tickpos <- seq(0, 1, .1)
	ticklabels <- format(tickpos, decimals = 2)
	
	for (i in seq(length(tickpos))) {
	    mtext(ticklabels[i], 1, 0.4*size/15, at = tickpos[i], cex = cexaxis)
	    mtext(ticklabels[i], 2, 0.6*size/15, at = tickpos[i], cex = cexaxis, las = 2)
	}
	
	if (mguides) {
        abline(v = .5, lty = 2, col="gray")
        abline(h = .5, lty = 2, col="gray")
    }
    
    abline(0, 1, col="gray")
    
    plotpoints <- list(as.name("points"), x, y, pch = pch, cex = cexpoints, bg = bgpoints)
    suppressWarnings(eval(as.call(c(plotpoints, other.args))))
    
    inclcov <- pof(xcopy, ycopy, relation = relation)$incl.cov
    
    mtext(paste(ifelse(QCA::nec(relation), "Necessity", "Sufficiency"), "relation"),
          line = 1.3, cex = cexaxis/0.8, font = 2)
    mtext(paste("Inclusion:", format(inclcov[, 1], digits = 3)), at = 0, adj = 0, cex = cexaxis)
    mtext(paste("Coverage:", format(inclcov[, 3], digits = 3)), at = 0.3, adj = 0, , cex = cexaxis)
    mtext(paste(ifelse(QCA::nec(relation), "Relevance:", "PRI:"),
          format(inclcov[, 2], digits = 3)), at = 0.6, adj = 0, , cex = cexaxis)
    
    return(invisible(list(x = x, y = y)))
    
}

