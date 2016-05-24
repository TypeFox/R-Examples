
################### TRANS - AXIS FUNCTION ####################
transAxis <- function (formula, subset = NULL, var.subset = NULL)
{

    # get the correct device for opening a plot window, note this is platform specific.
    dev.name <- getOption("device")


    if (deparse(terms(formula)[[2]], width.cutoff = 500) == eval(all.vars(eval(formula), unique = FALSE)[1])) {
        trans <- FALSE
        yaxlabs <- NULL
        yaxvalues <- NULL
    }
    else {
        y <- get(all.vars(eval(formula), unique = FALSE)[1])
        if (is.matrix(y))
            p <- dim(y)[2]
        else {
            p <- 1
            y <- as.matrix(y)
        }
        if (is.matrix(eval(terms(formula)[[2]]))) {
            q <- dim(eval(terms(formula)[[2]]))[2]
        }
        else q <- 1
        if (p != q) {
            trans <- FALSE
            yaxlabs <- NULL
            yaxvalues <- NULL
        }
        else {
            if (missing(var.subset)) {
                var.subset <- 1:p
            }
            else if (is.logical(var.subset) & any(!is.na(var.subset)))
                var.subset <- which(var.subset[!is.na(var.subset)])
            if (missing(subset))
                subset <- 1:nrow(y)
            trans <- TRUE
            y <- get(all.vars(eval(formula), unique = FALSE)[1])
            if (is.matrix(y))
                p <- dim(y)[2]
            else p <- 1
            origayticks <- list()
            yaxlabs <- list()
            yaxvalues <- list()
            miny <- list()
            maxy <- list()
            formula.term <- deparse(terms(formula)[[2]], width.cutoff = 500)
            while (regexpr(" ", formula.term) != -1) {
                formula.term <- sub(" ", "", formula.term)
            }

	    do.call(dev.name, args=list())

	    #Determine the transformation
	    if (regexpr("log\\(",formula.term)[1] != -1) transform <- "log" 
	    else if (regexpr("sqrt\\(",formula.term)[1] != -1 || regexpr("\\^0.5",formula.term)[1] != -1) transform <- "sqrt"
	    else if (regexpr("\\^0.25",formula.term)[1] != -1) transform <- "sqrt4"
	    else transform="no"

	    ymax <- max(y)
	    ymin <- min(y[y>0])
	    y.tick <- axisTicks(transform=transform, max=ymax, min=ymin, tran.lab="o")

            for (i in var.subset) {
                plot(y[subset, i], type = "n", axes = FALSE, xlab = "", ylab = "")
		
		seq.tick <- seq(1,length(axTicks(2)),by=2)
                origayticks[[i]] <- axTicks(2)#[seq.tick]
                
		miny[[i]] <- min(y[, i])
                maxy[[i]] <- max(y[, i])
                yaxlabs[[i]] <- sapply(origayticks[[i]], function(x) sub(eval(all.vars(eval(formula), unique = FALSE)[1]), x, formula.term))
                yaxvalues[[i]] <- sapply(yaxlabs[[i]], function(x) eval(formula(paste("f~",x, sep = ""))[[3]]))
                yaxlabs[[i]] <- c(y.tick$x.ticlab, "")
                yaxvalues[[i]] <- c(y.tick$x.tic,maxy[[i]])
            }
            dev.off()

        }
    }
    list(trans = trans, yaxvalues = yaxvalues, yaxlabs = yaxlabs)
}


############# CREATE TICK-MARKS FOR AXIS SCALING ############## 
axisTicks <- function (transform, max, min, tran.lab) 
{
        #Transform
        if (transform == "log") {
                c.log <- 5/log(max/min + 1)
                seqTicLab <-round(min*(exp((0:5)/c.log)-1))
                seqTic <- log(seqTicLab/min + 1) 
		if (tran.lab == "t") seqTicLab <- round(seqTic,1)

        } else if (transform == "sqrt") {
                c.sqr <- max/25
                seqTicLab <- round((0:5)^2 * c.sqr)
                seqTic <- sqrt(seqTicLab)
		if (tran.lab == "t") seqTicLab <- round(seqTic,1)

        } else if (transform == "sqrt4") {
                c.4th <- max/5^4
                seqTicLab <- round((0:5)^4 * c.4th)
                seqTic <- sqrt(sqrt(seqTicLab))
		if (tran.lab == "t") seqTicLab <- round(seqTic,1)

        } else if (transform == "no") {
                c <- 5/max
                seqTicLab <- round((0:5)/c)
                seqTic <- seqTicLab
		if (tran.lab == "t") seqTicLab <- round(seqTic,1)
        } else {
                
        }

        df <- data.frame(x.tic= seqTic, x.ticlab= seqTicLab)
        return(df)
}