## ====================================================
## Generic function for plotting state sequence objects
## ====================================================

seqplot <- function(seqdata, group=NULL, type="i", title=NULL,
	cpal=NULL, missing.color=NULL,
	ylab=NULL, yaxis=TRUE, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1,
	use.layout=(!is.null(group) | withlegend!=FALSE), legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(seqdata,"stslist"))
		stop(call.=FALSE, "seqplot: data is not a sequence object, use seqdef function to create one")

	## Storing original optional arguments list
	oolist <- list(...)

	## ==============================
	## Preparing if group is not null
	## ==============================

        if (type == "pc") { # modification of Reto Buergin 16.08.2012
          oolist <- append(oolist, list(group = group, rows = rows, cols = cols))
          group <- NULL
        }

	if (!is.null(group)) {
          group <- group(group)

          ## Check length
          if (length(group)!=nrow(seqdata))
            stop(call.=FALSE, "group must contain one value for each row in the sequence object")

          nplot <- length(levels(group))
          gindex <- vector("list",nplot)

          for (s in 1:nplot)
            gindex[[s]] <- which(group==levels(group)[s])

          ## Title of each plot
          if (!is.null(title))
            title <- paste(title,"-",levels(group))
          else
            title <- levels(group)
	} else {
          nplot <- 1
          gindex <- vector("list",1)
          gindex[[1]] <- 1:nrow(seqdata)
	}

	## ===================
	## Defining the layout
	## ===================
	if (type=="Ht" | type =="pc") { withlegend=FALSE }

	## IF xaxis argument is provided
	## it interferes with axes argument
	if ("xaxis" %in% names(oolist)) {
		tmpxaxis <- oolist[["xaxis"]]
		if (tmpxaxis==TRUE) {axes="all"}
		else if (tmpxaxis==FALSE) {axes=FALSE}
		oolist <- oolist[!names(oolist) %in% "xaxis"]
	}

	if (use.layout | !is.null(group) ) {
		## Saving graphical parameters
		savepar <- par(no.readonly = TRUE)

		lout <- TraMineR.setlayout(nplot, rows, cols, withlegend, axes, legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)

		## Axis should be plotted or not ?
		xaxis <- 1:nplot==lout$axisp

		legpos <- lout$legpos
	}
	else {
		if (axes!=FALSE) {xaxis <- TRUE}
		else {xaxis <- FALSE}
		legpos <- NULL
	}

	## =======
	## Ploting
	## =======
	for (np in 1:nplot) {
		## Storing ... arguments in a list
		olist <- oolist
		if ("sortv" %in% names(olist)) {sortv <- olist[["sortv"]]}
		if ("dist.matrix" %in% names(olist)) {dist.matrix <- olist[["dist.matrix"]]}
		if ("with.miss" %in% names(olist)) {
			missidx <- which(names(olist)=="with.miss")
			names(olist)[missidx] <- "with.missing"
			message(" [i] argument 'with.miss' is obsolete and replaced by 'with.missing'")
		}

		plist <- list(main=title[np], cpal=cpal, missing.color=missing.color,
			ylab=ylab, yaxis=yaxis, xaxis=xaxis[np],
			xtlab=xtlab, cex.plot=cex.plot)

		## Selecting sub sample for x
		## according to 'group'
		subdata <- seqdata[gindex[[np]],]

		## State distribution plot or Entropy index
		if (type=="d" || type=="Ht") {
			f <- seqstatd
			plist <- c(list(type=type), plist)

			plist <- plist[!names(plist) %in% "missing.color"]

			## Removing the 'cpal' argument which is not used
			## in Entropy index plots
			if (type=="Ht") {plist <- plist[!names(plist) %in% "cpal"]}
		}
		## Sequence frequency plot
		else if (type=="f") {
			with.missing <- TRUE
			f <- seqtab
		}
		## Sequence index plot
		else if (type=="i" || type=="I") {
			f <- function(seqdata) {return(seqdata)}
			with.missing <- TRUE

			## Selecting sub sample for sort variable
			## according to 'group'
			if ("sortv" %in% names(olist)) {
				if (!length(sortv)==1) {
					olist[["sortv"]] <- sortv[gindex[[np]]]
				}
			}

			if (type=="I") {
				if (!"tlim" %in% names(olist)) {olist <- c(olist, list(tlim=0))}
				if (!"space" %in% names(olist)) {olist <- c(olist, list(space=0))}
				if (!"border" %in% names(olist)) {olist <- c(olist, list(border=NA))}
			}
		}
		## Mean times
		else if (type=="mt") {f <- seqmeant}
		## Mean times
		else if (type=="ms") {
			f <- seqmodst
		}
		## Representative sequence
		else if (type=="r") {
			f <- seqrep
			with.missing <- TRUE

			## Removing unused arguments
			plist <- plist[!names(plist) %in% "yaxis"]
			
			## Selecting distances according to group
			if (!"dist.matrix" %in% names(olist))
				stop("You must provide a distance matrix", call.=FALSE)
			else {
				if (inherits(dist.matrix, "dist")) {
        				dist.matrix <- dist2matrix(dist.matrix)
				}

				olist[["dist.matrix"]] <- dist.matrix[gindex[[np]],gindex[[np]]]

				## Max theoretical distance for the scale
				if (!"dmax" %in% names(olist)) {
					dmax <- max(dist.matrix)
					olist <- c(olist,list(dmax=dmax))
				}
			}

                      } else if (type == "pc") { # modification of Reto Buergin 16.08.2012

                        plist$title <- title
                        olist <- c(olist, plist)
                        olist$plot <- FALSE
                        f <- seqpcplot
                        olist <- olist[names(olist) %in% names(formals(f))]
                        plist <- list()
                      }
		else
			stop("Unknown 'type' argument.")		

		## Calling appropriate function and plotting
		flist <- names(formals(f))

		if ("with.missing" %in% names(olist)) {
			with.missing <- olist[["with.missing"]]
		} else if ("with.missing" %in% flist) {
			with.missing <- formals(f)$with.missing
		}

		## Xlim when plotting individual sequences
		if (type %in% c("i", "I", "f")) {
			if (!"xlim" %in% names(olist)) {
				olist <- c(olist, list(xlim=c(0, ncol(seqdata))))
			}
		}

		match.args <- names(olist) %in% flist
		fargs <- olist[match.args]
		fargs <- c(list(seqdata=subdata), fargs)
		res <- do.call(f, args=fargs)

		olist <- olist[!match.args]
		plist <- c(list(x=res), plist, olist)
		do.call(plot, args=plist)
	}	

	## Plotting the legend
	if (!is.null(legpos)) {
		## Extracting some sequence characteristics
		nr <- attr(seqdata,"nr")

		if (is.null(ltext)) ltext <- attr(seqdata,"labels")

		if (is.null(missing.color)) missing.color <- attr(seqdata,"missing.color")

		if (is.null(cpal)) cpal <- attr(seqdata,"cpal")

		density <- if ("density" %in% names(oolist)) { oolist[["density"]] } else { NULL }
		angle <- if ("angle" %in% names(oolist)) { oolist[["angle"]] } else { NULL }

		## Adding an entry for missing in the legend
		if (with.missing & any(seqdata==nr)) {
			cpal <- c(cpal,missing.color)
			ltext <- c(ltext,"missing")
		## statl <- c(statl,nr)
		## nbstat <- nbstat+1
		}

		TraMineR.legend(legpos, ltext, cpal, cex=cex.legend, density=density, angle=angle)
	}

	## Restoring graphical parameters
	if (use.layout | !is.null(group)) {par(savepar)}
}
