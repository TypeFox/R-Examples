subset.DVH.list <- function(x, structure=NULL, patient=NULL, ID=NULL, constraint=NULL, ignore.case=TRUE, select=NULL) {
	select <- match.arg(select, choices=c("all", "any", "none"))
	switch(select,
		all = results <- rep(TRUE, length(x)),
		any = results <- rep(FALSE, length(x)),
		none = results <- rep(FALSE, length(x))
	)
	if (!is.null(structure)) {
		if (ignore.case) {
			structure <- paste("(?i)", structure, sep="")
		}
		results.structure <- rep(FALSE, length(x))
		for (i in structure) {
			results.structure <- as.logical(results.structure + grepl(i, unlist(lapply(x, names))))
		}
		switch(select,
			all = results <- as.logical(results * results.structure),
			any = results <- as.logical(results + results.structure),
			none = results <- as.logical(results + results.structure)
		)
	}
	if (!is.null(patient)) {
		if (ignore.case) {
			patient <- paste("(?i)", patient, sep="")
		}
		results.patient <- rep(FALSE, length(x))
		for (i in patient) {
			results.patient <- as.logical(results.patient + grepl(i, as.character(x$patients)))
		}
		switch(select,
			all = results <- as.logical(results * results.patient),
			any = results <- as.logical(results + results.patient),
			none = results <- as.logical(results + results.patient)
		)
	}
	if (!is.null(ID)) {
		results.ID <- rep(FALSE, length(x))
		for (i in ID) {
			results.ID <- as.logical(results.ID + (as.character(x$ID) == i))
		}
		switch(select,
			all = results <- as.logical(results * results.ID),
			any = results <- as.logical(results + results.ID),
			none = results <- as.logical(results + results.ID)
		)
	}
	if (!is.null(constraint)) {
		results.constraint <- rep(FALSE, length(x))
		for (i in constraint) {
			results.constraint <- unlist(lapply(x, function(y) {y <- y[i]; y * (names(y) == "logical")}))
			if (any(is.na(results.constraint))) {
				results.constraint[which(is.na(results.constraint))] <- FALSE
			}
			switch(select,
				all = results <- as.logical(results * results.constraint),
				any = results <- as.logical(results + results.constraint),
				none = results <- as.logical(results + results.constraint)
			)
		}
	}
	switch(select,
		all = return(x[results]),
		any = return(x[results]),
		none = return(x[!results])
	)		
}
	
setGeneric("subset",
	subset
)
	
setMethod("subset", "DVH",
	function (x, ...) {
		subset.DVH.list(as(x, "DVH.list"), ...)
	}
)	

setMethod("subset", "DVH.list",
	function (x, ...) {
		subset.DVH.list(x, ...)
	}
)