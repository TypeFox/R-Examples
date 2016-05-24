craer <- function(x, type = c("regular", "fast"), rotate = FALSE, loss, loss.args = NULL, interp = "bicubic",
    method = "BFGS", stages = TRUE, verbose = FALSE, ...) { 

    if(verbose) begin.tiid <- Sys.time()

    if(class(x) != "matched") stop("craer: invalid x argument.  Must be of class matched.")

    type <- tolower(type)
    type <- match.arg(type)

    if(missing(loss) && type == "regular") loss <- "QlossRigid"
    else if(missing(loss)) loss <- NULL

    n <- dim(x$matches)[ 1 ]

    # If no matches, quietly return NULL.
    if(n == 0) return(NULL)

    a <- attributes(x)
    a$names <- NULL

    Xfeats <- x$X.feats
    Yfeats <- x$Y.feats

    # obtain rigid transformations either by optimizing a loss function (regular), or
    # by only calculating the centroid and possibly the orientation angle differences (fast).

    # loc <- a$loc
    loc <- cbind(rep(1:a$xdim[ 1 ], a$xdim[ 2 ]), rep(1:a$xdim[ 2 ], each = a$xdim[ 1 ]))

    rfun <- function(ij, Xfeats, Yfeats, loc, typ, rot, leas, leas.args, itp, omet, stg, vb, ...) {

	i <- ij[ 1 ]
	j <- ij[ 2 ]

	if(vb) cat("\n", "Comparing observed feature ", j, " to forecast feature ", i, "\n")

	obj <- rigider(x1 = as.matrix(Yfeats[[ i ]]), x0 = as.matrix(Xfeats[[ j ]]), p0 = loc, type = typ,
			rotate = rot, loss = leas, loss.args = leas.args, interp = itp, method = omet,
			stages = stg, verbose = vb, ...)

	return(obj)

    } # end of internal 'rfun' function.

    if(verbose) cat("\n", "Applying rigid transformations.  This may take time even if type is fast.\n")

    transforms <- apply(x$matches, 1, rfun, Xfeats = Xfeats, Yfeats = Yfeats, loc = loc, typ = type, rot = rotate,
		leas = loss, leas.args = loss.args, itp = interp, omet = method, stg = stages, vb = verbose, ...)


    sfun <- function(x, stages, rotate) {

	par <- x$par
	parnom <- names(par)

	x0 <- x$x0
	x1 <- x$x1

	d2 <- (x1 - x0)^2
	N <- sum(colSums(!is.na(d2)), na.rm = TRUE)
	MSE.total <- sum(colSums(d2, na.rm = TRUE), na.rm = TRUE) / N

	if(stages && rotate) {

	    x1.shift <- x$x1.translated
	    x1.shift.rot <- x$x1.transformed

	    d2.shift <- (x1.shift - x0)^2
	    N.shift <- sum(colSums(!is.na(d2.shift)))
	    MSE.shift <- sum(colSums(d2.shift, na.rm = TRUE), na.rm = TRUE) / N.shift

	    MSE.translation <- MSE.total - MSE.shift

	    d2.shift.rot <- (x1.shift.rot - x0)^2
	    N.shift.rot <- sum(colSums(!is.na(d2.shift.rot)))
	    MSE.shift.rot <- sum(colSums(d2.shift.rot, na.rm = TRUE), na.rm = TRUE) / N.shift.rot

	    fm <- sum(colSums(x1.shift.rot, na.rm = TRUE), na.rm = TRUE) / sum(colSums(!is.na(x1.shift.rot), na.rm = TRUE), na.rm = TRUE)
            om <- sum(colSums(x0, na.rm = TRUE), na.rm = TRUE) / sum( colSums( !is.na( x0 ), na.rm = TRUE), na.rm = TRUE)
            MSE.volume <- (fm - om)^2


	    res <- c(par, MSE.total, MSE.shift, MSE.translation, MSE.shift.rot, MSE.translation - MSE.shift.rot, MSE.volume,
			MSE.shift.rot - MSE.volume)
	    names(res) <- c(parnom, "MSE.total", "MSE.shift", "MSE.translation", "MSE.displacement", "MSE.rotation",
				"MSE.volume", "MSE.pattern")

	} else if(rotate) {

	    x1.shift.rot <- x$x1.transformed
	    d2.shift.rot <- (x1.shift.rot - x0)^2
            N.shift.rot <- sum(colSums(!is.na(d2.shift.rot)))
            MSE.shift.rot <- sum(colSums(d2.shift.rot, na.rm = TRUE), na.rm = TRUE) / N.shift.rot

	    fm <- sum(colSums(x1.shift.rot, na.rm = TRUE), na.rm = TRUE) / sum(colSums(!is.na(x1.shift.rot), na.rm = TRUE), na.rm = TRUE)
            om <- sum(colSums(x0, na.rm = TRUE), na.rm = TRUE) / sum( colSums( !is.na( x0 ), na.rm = TRUE), na.rm = TRUE)
            MSE.volume <- (fm - om)^2

	    res <- c(par, MSE.total, MSE.shift.rot, MSE.total - MSE.shift.rot, MSE.volume, MSE.shift.rot - MSE.volume)
	    names(res) <- c(parnom, "MSE.total", "MSE.shift.rot", "MSE.displacement", "MSE.volume", "MSE.pattern")

	} else {

	    # Note: Does not allow option to not translate.

	    x1.shift <- x$x1.transformed

	    d2.shift <- (x1.shift - x0)^2
            N.shift <- sum(colSums(!is.na(d2.shift)))
            MSE.shift <- sum(colSums(d2.shift, na.rm = TRUE), na.rm = TRUE) / N.shift

            MSE.translation <- MSE.total - MSE.shift
	    fm <- sum(colSums(x1.shift, na.rm = TRUE), na.rm = TRUE) / sum( colSums( !is.na( x1.shift ), na.rm = TRUE), na.rm = TRUE)
	    om <- sum(colSums(x0, na.rm = TRUE), na.rm = TRUE) / sum( colSums( !is.na( x0 ), na.rm = TRUE), na.rm = TRUE)
	    MSE.volume <- (fm - om)^2

	    res <- c(par, MSE.total, MSE.shift, MSE.translation, MSE.volume, MSE.shift - MSE.volume)
	    names(res) <- c(parnom, "MSE.total", "MSE.shift", "MSE.displacement", "MSE.volume", "MSE.pattern")

	}

	return(res)

    } # end of 'sfun' function.

    if(verbose) cat("\n", "Rigid transforms found.  Calculating MSE and its various breakdowns.\n")

    out <- lapply(transforms, sfun, stages = stages, rotate = rotate)
    # TO DO: see what happens here if there is only one matched object.
    nomen <- names(out[[ 1 ]])
    out <- matrix(unlist(out), nrow = n, byrow = TRUE)
    colnames( out ) <- nomen

    attr(out, "number.of.matches") <- n
    attr(out, "information") <- a

    class(out) <- "craered"

    if(verbose) print(Sys.time() - begin.tiid)
    return(out)

} # end of 'craer' function.

print.craered <- function(x, ...) {

    a <- attributes(x)
    x <- matrix(c( x ), nrow = a$number.of.matches)
    colnames( x ) <- a$dimnames[[ 2 ]]

    info <- a$information

    cat("\n", "CRA results for ", info$data.name, "\n")
    cat("model ", info$model, " at time ", info$time.point, "\n")
    cat(info$msg, "\n\n")

    print(x)

    invisible()

} # end of 'print.craered' function.
