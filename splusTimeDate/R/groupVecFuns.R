"groupVec" <- 
function(names, classes, columns)
{
  # Constructor function for a groupVec
	# Arguments are a vector of names of its columns, a vector of the classes
	# corresponding to names, and the columns (optional)
	if(missing(names) && missing(classes) && missing(columns)) return(new(
			"groupVec"))
	if(missing(names))
		stop("names argument missing with no default")
	if(missing(classes))
		classes <- rep("ANY", length(names))
	if(missing(columns))
		columns = lapply(classes, "new")
	l <- length(columns)
	if(l != length(names) || l != length(classes))
		stop("columns, class, and name arguments must be same length")
	if(l > 0) {
		for(i in 1:l)
			if(!is(columns[[i]], classes[i])) stop(
					"columns does not have correct classes")
		ls <- sapply(columns, "length")
		l <- ls[1]
		if(!all(ls == l))
			stop("All columns must be same length")
	}
	# create an object with desired variables and return it
	obj = new("groupVec", columns = columns, names = names, classes = classes)
	obj
}

"groupVecClasses<-" <- 
function(object, value)
{
  ## classes replacement method function for groupVec object
  ## replace the data column classes
	# if the new classes have different length, the names and data will
	# be truncated/extended appropriately
	object@classes <- as(value, "ANY")
	lnew <- length(object@classes)
	lold <- length(object@names)
	if(lnew > lold) {
		# we have added to classes
		# add to names and data appropriately
		length(object@names) <- lnew
		newclasses <- object@classes[(lold + 1):lnew]
		datalen <- length(object)
		newdata <- lapply(lapply(newclasses, "new"), "length<-",
			datalen)
		object@columns <- c(object@columns, newdata)
	}
	else {
		# new length is shorter, so truncate
		length(object@names) <- lnew
		length(object@columns) <- lnew
	}
	# coerce the data to the new classes
	object@columns <- lapply(1:length(object@columns), function(colnum, dat,
		cls)
	as(dat[[colnum]], cls[colnum]), object@columns, object@classes)
	object
}

"groupVecClasses" <- 
function(object)
{
  # classes method function for groupVec object
	# extract the data column classes
	object@classes
}

"groupVecColumn<-" <- 
function(object, col.name, value)
{
  # Column replacement method function for groupVec object
	# Replace data from the col.name column(s) of object with value
	positions <- match(col.name, object@names, nomatch = 0)
	if(!all(positions))
		warning("non-matching column names will be ignored")
	if(!any(positions))
		return(object)
	if(length(col.name) == 1)
		value <- list(value)
	object@columns[positions] <- lapply(1:length(value), function(colnum,
		dat, cls)
	as(dat[[colnum]], cls[colnum]), value, object@classes)
	# check lengths
	ls <- sapply(object@columns, "length")
	l <- ls[1]
	if(!all(ls == l))
		stop("All columns must be same length")
	object
}

"groupVecColumn" <- 
function(object, col.name)
{
  # Column extraction function for groupVec object
	# Extract data from the col.name column(s) of object
	positions <- match(col.name, object@names, nomatch = 0)
	if(!all(positions))
		warning("non-matching column names will be ignored")
	if(length(positions) == 1)
		object@columns[[match(col.name, object@names, nomatch = 0)
			]]
	else object@columns[match(col.name, object@names, nomatch = 0)]
}

"groupVecData<-" <- 
function(object, value)
{
  # data method function for groupVec object
	# replace the data
	# if the new data list has different length, the names and classes
	# will be truncated/extended appropriately
  lold <- length(object@names)
	object@columns <- value
	lnew <- length(object@columns)
	if(lnew > lold) {
		# we have added to data
		# add to names and classes appropriately
		length(object@names) <- lnew
		newclasses <- sapply(object@columns, "class")[(lold + 1):lnew]
		object@classes <- c(object@classes, newclasses)
	}
	else {
		# new length is shorter, so truncate
		length(object@names) <- lnew
		length(object@classes) <- lnew
	}
	l <- length(object@columns)
	if(l > 0) {
		for(i in 1:l)
			if(!is(object@columns[[i]], object@classes[i])) 
					stop(
					"data does not have correct classes")
		ls <- sapply(object@columns, "length")
		l <- ls[1]
		if(!all(ls == l))
			stop("All columns must be same length")
	}
	object
}

"groupVecData" <- 
function(object)
{
  # data method function for groupVec object
	# extract the data
	object@columns
}

"groupVecExtValid" <- 
function(object, names, classes, checkrest = FALSE)
{
  ## validity function especially for classes that extend groupVec.
  ## It checks that the object is a valid groupVec first, and then
  ## compares the column names and classes to the passed in arguments, 
  ## and if checkrest is TRUE, also checks to see that the slots of object
  ## that did not come from the groupVec class are of length <= 1
  GVret <- groupVecValid(object)
  if(GVret != TRUE)
    return(GVret)
  ## check names
  tmpnames <- groupVecNames(object)
  if((length(tmpnames) != length(names)) || (all.equal(tmpnames, names) !=
              TRUE))
    return("groupVec object column names are incorrect")
  ## check classes
  tmpnames <- as(groupVecClasses(object), "character")
  if((length(tmpnames) != length(classes)) || (all.equal(tmpnames, 
              classes) != TRUE))
    return("groupVec object column classes are incorrect")
  ## check that rest of slots are length <= 1, if desired
  if(checkrest) return(groupVecNonVec(object, names(getSlots("groupVec"))
                                      ))
  TRUE
}

"groupVecNames<-" <- 
function(object, value)
{
  # names method replacement function for groupVec object
	# replace the data column names
	# if the new names have different length, the classes and data will
	# be truncated/extended appropriately (extended with numerics)
	object@names <- as(value, "character")
	lnew <- length(object@names)
	lold <- length(object@classes)
	if(lnew > lold) {
		# we have added to names.  For now, add numerics to end
		# of classes and data; make sure the numeric data vectors
		# are the right length too
		newclasses <- rep("numeric", lnew - lold)
		datalen <- length(object)
		newdata <- rep(list(numeric(datalen)), lnew - lold)
		object@classes <- c(object@classes, newclasses)
		object@columns <- c(object@columns, newdata)
	}
	else {
		# new lengths are shorter, so truncate classes and data
		length(object@classes) <- lnew
		length(object@columns) <- lnew
	}
	object
}

"groupVecNames" <- 
function(object)
{
  # names method function for groupVec object
	# extract the data column names
	object@names
}

"groupVecNonVec" <- 
function(object, exceptSlots)
{
  ## Function that checks that all slots in Object except those
  ## in exceptSlots have length <= 1 if they are vectors
  ## return TRUE if all slots (except those in optional argument
  ## exceptSlots) are of length 1 or less, and
  ## an error string if that's not true
  objSlots <- names(getSlots(class(object)))
  ## get required subset of slots
  if(!missing(exceptSlots))
    objSlots <- objSlots[is.na(match(objSlots, exceptSlots))]
  ## check these slots -- must either be non-vector or have length <= 1
  for(slname in objSlots) {
    sldata <- slot(object, slname)
    if(is(sldata, "vector") && (length(sldata) > 1))
      return(paste("Slot", slname, "has length > 1"))
  }
  TRUE
}

"groupVecValid" <- 
function(object)
{
  ## validity function for groupVec objects
  ## check that lengths of components match
  ## check that classes of columns components are correct
  ## check that all vector lengths are the same
  if(!is(object, "groupVec")) return("Object is not a groupVec")
  ## check lengths of components
  l <- length(object@columns)
  if(l != length(object@names) || l != length(object@classes)
     )
    return(paste("groupVec name, class, and data component", 
                 "lengths do not match"))
  if(l == 0)
    return(TRUE)
  ## check classes of columns elements
  for(i in 1:l)
    if(!is(object@columns[[i]], object@classes[i]))
      return( paste("groupVec data class ", i, " should be ",
                    object@classes[i], " not ", class(object@columns[[i]])))
  ## check that all vector lengths are the same
  ls <- sapply(object@columns, "length")
  l <- ls[1]
  if(!all(ls == l))
    return("groupVec data lengths not all the same")
  TRUE
}
