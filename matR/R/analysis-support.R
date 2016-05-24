
#-----------------------------------------------------------------------------------------
#  functions in support of analysis methods.
#-----------------------------------------------------------------------------------------


resolve <- function (first, second) {
#-----------------------------------------------------------------------------------------
#  merge two named lists, eliminating duplicates and giving priority to the first
#-----------------------------------------------------------------------------------------
	append (first, second) [ !duplicated (c (names(first), names(second))) ]
	}


parMap <- function (object, name.map, value.map) {
#-----------------------------------------------------------------------------------------
#  workhorse function enabling metadata mapping in graphics.  input is:
#    object			biom
#    map			mapping of par variable names to metadata field names (character)
#    values			optional specific value mappings for par variable (list)
#  return value is specific values for all mapped par variables.
#-----------------------------------------------------------------------------------------
	if (is.null (name.map)) return (list())

	pars <- names (name.map)
	metadata.values <- lapply (name.map, getMetColumns, object)

####  make NA any metadata that is "NA".

	metadata.values <- lapply (metadata.values,
		function (xx) {
			levels (xx) [levels (xx) == "NA"] <- NA
			xx
			})

####  create mappings where they are not specified at all.
####  in the mapping, the missing level must be identified by name "NA",
####  not by a missing name.

	automatic <- setdiff (pars, names (value.map))
	value.map [automatic] <- mapply(
		function (par, met) {
			met <- addNA (met, ifany=TRUE)
			yy <- parAuto (par, nlevels (met))
			names (yy) <- levels (met)
			names (yy) [is.na (names (yy))] <- "NA"
			yy
			},
		as.list (automatic),
		metadata.values [automatic],
		SIMPLIFY=FALSE)

	value.map <- value.map [pars]

####  make NA any unspecified levels of partially-specified metadata

	metadata.values <- mapply(
		function (xx, map) {
			levels (xx) [! (levels (xx) %in% names (map))] <- NA
			xx
			},
		metadata.values,
		value.map,
		SIMPLIFY=FALSE)

####  complete partial mappings by adding a value for NA, if necessary

	value.map <- mapply(
		function (map, par) {
			if (! ("NA" %in% names (map))) {
				c(map, "NA" = parAuto (par))
			} else map },
		value.map,
		names (value.map),
		SIMPLIFY=FALSE)

####  apply all mappings to create par specifications.
####  the missing level must be renamed to "NA",
####  in order to get the value intended for it.

	yy <- mapply(
		function (values, map) {
			values <- addNA(values,ifany=TRUE)
			yy <- levels (values)
			yy [is.na (yy)] <- "NA"
			levels (values) <- map [yy]
			if (is.numeric (map)) {								# ensure correct class
				as.numeric (as.character (values))
			} else as.character (values)
			},
		metadata.values,
		value.map,
		SIMPLIFY=FALSE)
	}


#-----------------------------------------------------------------------------------------
#  utilities to facilitate metadata mapping
#-----------------------------------------------------------------------------------------

getMetColumns <- function (name, xx) {
	yy <- columns (xx, name)
	if (ncol (yy) != 1)
		stop ("\'", name, "\' does not identify a unique column metadata field")
	yy [[1]]
	}

getMetRows <- function (name, xx) {
	yy <- rows (xx, name)
	if (ncol (yy) != 1)
		stop ("\'", name, "\' does not identify a unique row metadata field")
	yy [[1]]
	}


#-----------------------------------------------------------------------------------------
#  utilities to facilitate literal metadata substitutions
#-----------------------------------------------------------------------------------------

subColumn <- function (name, xx) {
	if (length (name) != 1 ||
		!is.character (name) ||
		substr (name, 1, 2) != "$$") return (name)
	yy <- columns (xx, substring (name, 3))
	if (ncol (yy) != 1)
		stop ("\'", name, "\' does not identify a unique column metadata field")
	as.character (yy [[1]])
	}

subRow <- function (name, xx) {
	if (length (name) != 1 ||
		!is.character (name) ||
		substr (name, 1, 2) != "$$") return (name)
	yy <- rows (xx, substring (name, 3))
	if (ncol (yy) != 1)
		stop ("\'", name, "\' does not identify a unique row metadata field")
	as.character (yy [[1]])
	}
