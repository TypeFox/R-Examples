get.wpp.e0.data <- function(sex='M', start.year=1950, present.year=2015, wpp.year=2015, my.e0.file=NULL, 
							my.locations.file=NULL, verbose=FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNe0(sex=sex, wpp.year=wpp.year, my.e0.file=my.e0.file, 
								present.year=present.year, verbose=verbose)
	data <- un.object$data.object$data
	# get region and area data
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=wpp.year, my.locations.file=my.locations.file,
											package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	data_incl <- data[include,]
	nr_countries_estimation <- nrow(data_incl)
	if(any(!is.na(prediction.only))) { # move prediction countries at the end of data
		data_prediction <- data[prediction.only,]
		data_incl <- rbind(data_incl, data_prediction)
	}
	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=start.year, 
							present.year=present.year)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')

	LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									start.year, present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
									
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nr_countries_estimation,
				suppl.data=bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
				))
}

read.UNe0 <- function(sex, wpp.year, my.e0.file=NULL, ...) {
	un.dataset <- paste('e0', sex, sep='')
	un.suppl.dataset <- paste('e0', sex, '_supplemental', sep='')
	data <- bayesTFR:::do.read.un.file(un.dataset, wpp.year, my.file=my.e0.file, ...)
	suppl.data <- bayesTFR:::do.read.un.file(un.suppl.dataset, wpp.year, my.file=my.e0.file, ...)
	if(is.null(suppl.data$data)) suppl.data <- NULL
	return(list(data.object=data, suppl.data.object=suppl.data))
}

set.e0.wpp.extra <- function(meta, countries=NULL, my.e0.file=NULL, my.locations.file=NULL, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	un.object <- read.UNe0(sex=meta$sex, wpp.year=meta$wpp.year, my.e0.file=my.e0.file, 
							present.year=meta$present.year, verbose=verbose)
	data <- un.object$data.object
	extra.wpp <- bayesTFR:::.extra.matrix.regions(data=data, countries=countries, meta=meta, 
							package="bayesLife", verbose=verbose)
	if(!is.null(extra.wpp)) {
		extra.wpp <- list(e0.matrix=extra.wpp$tfr_matrix, 
						  e0.matrix.all=extra.wpp$tfr_matrix_all, 
						  regions=extra.wpp$regions, 
						  nr.countries.estimation=extra.wpp$nr_countries_estimation,
						  is_processed = extra.wpp$is_processed)
		locations <- bayesTFR:::read.UNlocations(data$data, wpp.year=meta$wpp.year, 
									my.locations.file=my.locations.file, package='bayesLife', verbose=verbose)
		suppl.wpp <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, extra.wpp, locations$loc_data, 
									meta$start.year, meta$present.year)
		extra.wpp$suppl.data <- bayesTFR:::.get.suppl.data.list(suppl.wpp, matrix.name='e0.matrix')
	}
	return(extra.wpp)
}

get.wpp.e0.data.for.countries <- function(meta, sex='M', my.e0.file=NULL, my.locations.file=NULL, verbose=FALSE) {
	sex <- toupper(sex)
	if(sex != 'M' && sex != 'F')
		stop('Allowed values for argument "sex" are "M" and "F".')
	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNe0(sex=sex, wpp.year=meta$wpp.year, present.year=meta$present.year, 
						my.e0.file=my.e0.file, verbose=verbose)
	data <- un.object$data.object$data
	# get region and area data
	locations <- bayesTFR:::read.UNlocations(data, wpp.year=meta$wpp.year, 
							my.locations.file=my.locations.file, package='bayesLife', verbose=verbose)
	loc_data <- locations$loc_data
	include <- c()
	for (i in 1:length(meta$regions$country_code)) { # put countries into the same order as in meta
		loc_index <- which(data$country_code == meta$regions$country_code[i])
		if(length(loc_index) <= 0) 
			stop('Country ', meta$regions$country_code[i], ' not found.')
		include <- c(include, loc_index)
	}
	data_incl <- data[include,]	
	LEXmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
							data_incl, loc_data, 
							start.year=meta$start.year, 
							present.year=meta$present.year)
	if (verbose) 
		cat('Dimension of the e0 matrix:', dim(LEXmatrix.regions$obs_matrix), '\n')
	LEXmatrixsuppl.regions <- bayesTFR:::.get.suppl.matrix.and.regions(un.object, LEXmatrix.regions, loc_data, 
									meta$start.year, meta$present.year)
	if(!is.null(un.object$suppl.data.object) && verbose) 
		cat('Dimension of the supplemental e0 matrix:', dim(LEXmatrixsuppl.regions$obs_matrix), '\n')
											
	return(list(e0.matrix=LEXmatrix.regions$obs_matrix, 
				e0.matrix.all=LEXmatrix.regions$obs_matrix_all, 
				regions=LEXmatrix.regions$regions, 
				nr.countries.estimation=nrow(data_incl),
				suppl.data=bayesTFR:::.get.suppl.data.list(LEXmatrixsuppl.regions, matrix.name='e0.matrix')
				))
}
