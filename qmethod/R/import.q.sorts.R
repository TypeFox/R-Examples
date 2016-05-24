import.q.sorts <- function(q.sorts.dir, q.set, q.distribution, conditions=NULL, manual.lookup=NULL) {
  # Input validation (also see validation at the bottom!)
  if (!is.matrix(q.set)) {
    stop("The q.set specified is not a matrix.")
  }
  if (!is.vector(q.distribution)) {
   stop("The q.distribution specified is not a vector.")
  }
  if (!is.null(conditions) &  !is.vector(conditions)) {
    stop("The conditions specified are not a vector.")
  }
  if (!is.null(manual.lookup) & !is.matrix(manual.lookup)) {
    stop("The manual.lookup specified is not a matrix.")
  }
  if (!is.null(conditions)) {  # test conditions subdir only if there are conditions
    for (cond in conditions) {
      if (!file.exists(paste(q.sorts.dir, cond, sep="")))  # this must not have a trailing slash, file.exists does not like that on win http://r.789695.n4.nabble.com/file-exists-does-not-like-path-names-ending-in-td4683717.html
      {
        stop(
          paste(
            "Folder for condition",
            cond,
            "could not be found."
          )
        )
      }
    }
  }

  # Deal with no conditions
  if (is.null(conditions)) {
    conditions <- "only.one"
  }
  conditions <- factor(conditions) #  such as before, after as factors

  # Set up preliminary data structure ==========================================
  p.set <- c() #  p.set are the participants in q lingo, empty vector
  for (cond in conditions) {  # gather *all* participants for all conds
    p.set.cond <- list.files(  # gather people by listing files
      path = if (cond == "only.one"){  # if no conditions
        q.sorts.dir  # this is the path
      } else {  # if more conditions
        paste (  # here comes the path
          q.sorts.dir,
          cond,  # consider condition in path
          "/",
          sep = ""
        )
      },
      no.. = TRUE,  # no dotfiles
      pattern = "\\.csv$"  # only csv
    )
    p.set.cond <- file_path_sans_ext(p.set.cond) #  kill extensions
    p.set <- append(p.set, p.set.cond) # append vector
  }
  p.set <- unique(p.set)  # make participants unique, also for no-cond (just in case)

  # Set up empty array =========================================================
  q.sorts <- array(
    #  participants, conditions, items makes 3 dimensions, all of which integer
    data = , #  no such thing yet, so all NAs (that should be the baseline!)
    dim = c(
      nrow(q.set), #  number of items
      length(p.set), #  number of participants
      length(conditions) #  number of conditions
    ),
    dimnames = list( #  dims should be called accordingly ...
      rownames(q.set), #  those are the items by meaningful short names
      p.set, #  those are the participants
      conditions #  those are the conditions
    )
  )

  # Create lookup table =======================================================
	if (is.null(manual.lookup)) {  # automatic hashing, same as in make.cards
    lookup.table <- apply(  # replace every language field with its hash
      X = q.set,
      MARGIN = c(1,2),
      digest,
        algo = "crc32",
        serialize = FALSE
    )
	} else {  # manually entered lookup table
    lookup.table <- manual.lookup
	}
	if (any(duplicated(lookup.table))) {  # test lookup table
	  stop ("There are duplicate IDs in the lookup table.")
	}

	# Import Loops ==============================================================
	for (cond in conditions) {  # loop over the conditions (such as before, after)
		for (part in p.set) {  # loop over participants
			path <- paste(
        q.sorts.dir,
        if (cond!="only.one") {
          paste(
            cond,
            "/",
            sep = ""
          )
        },
        part,
        ".csv",
        sep = ""
      )  # establish path
			if (!file.exists(path)) {  # there may be missing cases for some conditions
				warning(
          paste(  # it's not a deal-breaker just a warning
					  "There is no file for",
					  part,
					  "under condition",
					  cond,
					  ". NAs remains in array.",
					  sep = " "
				  )
        )
			} else {
			  current.sort <- read.csv(path, # let's do one sort at a time
				  header = FALSE, #  colnames will do
				  stringsAsFactors = FALSE, #  would only add confusion
				  nrows = max(q.distribution), # stuff below is ignored (item feedback, scores etc.)
				  na.strings = "", #  empty cells become NAs
				  colClasses = "character"  # just to make sure R doesn't choke (mistakenly identified) irrational numbers :)
			  )
			  current.sort <- as.matrix(current.sort) #  because read.csv makes dataframe
        for (id in na.omit(as.vector(current.sort))) {  # loops over ids
          if (any(lookup.table == id)) {  # do we know the id in the current sort?
            row <- which(lookup.table == id, arr.ind=TRUE)[,1]  # where is it in the table?
            handle <- rownames(lookup.table)[row]  # what is the short handle?
          } else {
            stop(
              paste(
                "The qsort in",
                path,
                "contains id",
                id,
                "is not defined as per manual.lookup and was ignored.",
                "Check whether you defined manual.lookup argument as intended."
              )
            )
          }
			    current.sort[current.sort==id] <- handle  # reassign it, in both cases
        }
			  # Test content of current sort for consistency ==========================
			  if (!all(colSums(!is.na(current.sort)) == q.distribution)) {  # distr. ok?
  				stop(paste(
	  				"The qsort in",
		  			path,
			  		"does not conform to the set distribution.",
				  	"Offending columns:",
					  colSums(!is.na(current.sort)) == q.distribution,
					  sep = " "
				  ))
			  }
			  current.sort.handles <- c(as.vector(current.sort))  # preparation
			  current.sort.handles <- na.omit(current.sort.handles)  # NAs confuse the comparisons
			  all.handles <- rownames(q.set)  # preparation
			  in.sample.not.sort <- all.handles[!(all.handles %in% current.sort.handles)]  # all sampled items in sort?
        if (length(in.sample.not.sort) > 0) { # if there are any missing from sort
  				stop(  # must error out b/c of inconsistency
            paste(
		  			  "Some items from the sample cannot be found in",
			  		  path,
				  	  ". Missing items: "
            ),
            paste (  # give detail
              in.sample.not.sort,
              collapse = ", "
            )
				  )
        }
        in.sort.not.sample <- current.sort.handles[!(current.sort.handles %in% all.handles)]  # anything not in sample but in sort?
        if (length(in.sort.not.sample) > 0) {  # if there are items in sort, but not in sample
          stop(  # must error out b/c of inconsistency
            paste(
              "Some items from",
              path,
              "cannot be found in the sample",
              ".",
              "Orphaned items:"
            ),
            paste (  # give detail
              in.sort.not.sample,
              collapse = ", "
            )
          )
        }
			  # Transposing and aggregation of individual qsort =========================
			  for (item in rownames(q.set)) {  # loops over items
  				q.sorts[item,part,cond] <- which( #  write in qsort content of
	  				current.sort == item,  # participant cell of curr. looped item
		  			arr.ind = TRUE,  # return complete index (row, column)
			  		useNames = TRUE)[,2] # return only column
				  q.sorts[item,part,cond] <- as.integer(q.sorts[item,part,cond]-((length(q.distribution)+1)/2))
				  # make it into integers, revert to original score
			  }
		  }
		}
	}
  if (length(conditions) == 1) {
    q.sorts <- q.sorts[,,1]  # drops redundant dim for conditions
  }
	return(q.sorts)
}
