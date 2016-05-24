#ALG: itree takes this from rpart but with stop/warnings
# when cp is not defined.

#SCCS  @(#)printcp.s	1.6 01/20/00
# print out the cptable, along with some summary of the tree
printcp <- function(x, digits=getOption("digits")-2)
{
    if (!inherits(x, 'itree')) stop ("x must be an itree object")
    
    if(x$method %in% c("class_purity","class_extremes","regression_purity","regression_extremes"))
  		stop("cp not defined for this method!")
    
    if(!is.null(x$penalty)){
    	warning("cp is impacted by using penalties and is NOT comparable to unpenalized cp's.") 
    }
    
    cat(switch(x$method,anova = "\nRegression tree:\n" ,
			class = "\nClassification tree:\n" ,
			poisson="\nRates regression tree:\n",
			exp = "\nSurvival regression tree:\n")
        )

    if(!is.null(cl <- x$call)) {
	dput(cl, control=NULL)
	cat("\n")
    }
    frame <- x$frame
    leaves <- frame$var == "<leaf>"
    used <- unique(frame$var[!leaves])

    if(!is.null(used)) {
        cat("Variables actually used in tree construction:\n")
        print(sort(as.character(used)), quote=FALSE)
        cat("\n")
    }


    cat("Root node error: ", format(frame$dev[1L], digits=digits), '/',
        frame$n[1], ' = ',
        format(frame$dev[1L]/frame$n[1L], digits=digits),
        '\n\n', sep='')


    n <- x$frame$n
    omit <- x$na.action
    if (length(omit))
    cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep="")
    else cat("n=", n[1L], "\n\n")

    print (x$cptable, digits=digits)
    invisible(x$cptable)
}


