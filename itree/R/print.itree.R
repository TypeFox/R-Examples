#ALG: itree copies this direclty from rpart.

#SCCS  @(#)print.rpart.s	1.15 06/06/01
print.itree <- function(x, minlength=0, spaces=2, cp,
               digits=getOption("digits"), ...) {
    if(!inherits(x, "itree")) stop("Not legitimate itree object")
    if (!is.null(x$frame$splits)) x <- rpconvert(x)  #help for old objects

	method.int <- pmatch(x$method, c("anova", "poisson", "class", "exp",
					"regression_extremes","regression_purity",
					"class_extremes","class_purity"))
	is_classification <- (method.int %in% c(3,7,8))
	
    if (!missing(cp)) x <- prune.itree(x, cp=cp)
    frame <- x$frame
    ylevel <- attr(x, "ylevels")
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32L), collapse = "")
    #32 is the maximal depth
    if(length(node) > 1L) {
        indent <- substring(indent, 1L, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")

    tfun <- (x$functions)$print
    if (!is.null(tfun)) {
		
		#itree adds this next 'if'
		if(is_classification){
			newyval2 <- cbind(frame$yval,frame[,grep("wt.",colnames(frame))])
			colnames(newyval2)[1] <- "yval"
			newyval2 <- as.matrix(newyval2) #make everything numeric
		}
		if (!is_classification) #changed in from rpart's check for yval2
					yval <- tfun(frame$yval,  ylevel, digits)
			else    yval <- tfun(newyval2,  ylevel, digits) #newyval2 instead of frame$yval2
		}
    else yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, digits=digits, minlength=minlength, ...)
    n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
               yval, term)

    omit <- x$na.action
    if (length(omit))
    cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep="")
    else cat("n=", n[1L], "\n\n")

    #This is stolen, unabashedly, from print.tree
    if (x$method=="class")
         cat("node), split, n, loss, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval\n")
    cat("      * denotes terminal node\n\n")

    cat(z, sep = "\n")
    return(invisible(x))
    #end of the theft
    }
