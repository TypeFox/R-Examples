## ======================================================
## subtree: extracts a given group or position from a PST
## ======================================================

setMethod("subtree", "PSTf", function(object, group=NULL, position=NULL) {

	if (!is.null(group) & !object@segmented) {
		stop(" [!] this is not a segmented PST")
	}

	A <- alphabet(object)
	cpal <- cpal(object)
	labels <- stlab(object)
	data <- object@data
	cdata <- object@cdata


	if (!is.null(group)) {
		if (is.numeric(group)) {
			gid <- object@group==levels(object@group)[group]
		} else {
			gid <- object@group==group
		}
		data <- data[gid,]
		cdata <- if (has.cdata(object)) { cdata[gid,] } else { cdata }
	}

	object <- as(object, "list")

	for (i in length(object):1) {
		object[[i]] <- lapply(object[[i]], select.segment, group=group, position=position)
		remove <- unlist(lapply(object[[i]], is.null))
		object[[i]] <- object[[i]][!remove]
		if (length(object[[i]])==0) {
			object <- object[-i]
		}
	}

	object <- new("PSTf", object, data=data, cdata=cdata, 
		alphabet=A, cpal=cpal, labels=labels, segmented=FALSE, call=match.call())
	return(object)
}
)

