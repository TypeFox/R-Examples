utils::globalVariables("wpp.data.env")

wpp.explore <- function(wpp.year=NULL, host=NULL, ...) {
	if(!is.null(wpp.year)) set.wpp.year(wpp.year)
	if(missing(host)) host <- getOption("shiny.host", "0.0.0.0")
	shiny::runApp(system.file('explore', package='wppExplorer'), host = host, ...)
}

wpp.explore3d <- function(wpp.year=NULL) {
	if(!is.null(wpp.year)) set.wpp.year(wpp.year)
	shiny::runApp(system.file('bubbles', package='wppExplorer'))
}

get.available.wpps <- function() c(2008, 2010, 2012, 2015)
check.wpp.revision <- function(wpp.year) {
	if (!wpp.year %in% get.available.wpps())
		stop('wpp.year must be one of ', paste(get.available.wpps(), collapse=', '))
}

wpp.year.from.package.name <- function(package)
	return(as.integer(substr(package, 4, nchar(package))))

wpp.indicator <- function(what, ...) {
	data <- do.call(what, list(...))
	if(is.null(data)) return(NULL)
	merge.with.un.and.melt(data, what=what)
}

	
wpp.by.year <- function(data, year) {
  data <- data[data$Year == year,]
  data$Year <- NULL
  data
}

wpp.by.country <- function(data, country) {
  data <- data[data$charcode == country,]
  data$charcode <- NULL
  data
}

wpp.by.countries <- function(data, countries) {
  data <- data[data$charcode %in% countries,]
  data
}

set.wpp.year <- function(wpp.year) {
	check.wpp.revision(wpp.year)
	# cleanup the environment
	for (item in ls(wpp.data.env)) {
		if(!(item %in% c('indicators'))) rm(list=item, envir=wpp.data.env)
	}
	data('iso3166', envir=wpp.data.env)
	wpp.data.env$package <- paste('wpp', wpp.year, sep='')
	# Filter out non-used countries
	do.call('data', list("popM", package=wpp.data.env$package, envir=wpp.data.env))
	wpp.data.env$iso3166 <- wpp.data.env$iso3166[is.element(wpp.data.env$iso3166$uncode, wpp.data.env$popM$country_code),]
	cat('\nDefault WPP package set to', wpp.data.env$package,'.\n')
}

get.wpp.year <- function() as.integer(substr(wpp.data.env$package, 4,8))
 
tpop <- function(...) {
	# Create a dataset of total population
	if.not.exists.load('popM')
	if.not.exists.load('popF')
	tpop <- sumMFbycountry(wpp.data.env$popM, wpp.data.env$popF)
	if(wpp.year.from.package.name(wpp.data.env$package) > 2010) { #projection stored separately from observations
		if.not.exists.load('popMprojMed')
		if.not.exists.load('popFprojMed')
		tpopp <- sumMFbycountry(wpp.data.env$popMprojMed, wpp.data.env$popFprojMed)
		tpop <- merge(tpop, tpopp, by='country_code')
	}
	tpop
}

tpopF <- function(...) return(tpop.sex('F'))
tpopM <- function(...) return(tpop.sex('M'))

tpop.sex <- function(sex) {
	# Create a dataset of total population
	dataset <- paste('pop', sex, sep='')
	pop <- load.dataset.and.sum.by.country(dataset)
	if(wpp.year.from.package.name(wpp.data.env$package) > 2010) { #projection stored separately from observations
		dataset <- paste('pop', sex, 'projMed', sep='')
		popp <- load.dataset.and.sum.by.country(dataset)
		pop <- merge(pop, popp, by='country_code')
	}
	pop
}

mig <- function(...) {
	# Create a dataset of net migration
	if(wpp.year.from.package.name(wpp.data.env$package) <2015) { # sex- and age-specific migration available
		if.not.exists.load('migrationM')
		if.not.exists.load('migrationF')
		return(sumMFbycountry(wpp.data.env$migrationM, wpp.data.env$migrationF))
	}
	load.and.merge.datasets('migration', NULL) # total migration available
}

migrate <- function(...) {
	migcounts <- mig()
	pop <- tpop()
	mergepop <- merge(migcounts[,'country_code', drop=FALSE], pop, sort=FALSE)
	ncols <- ncol(mergepop)
	#browser()
	cbind(country_code=mergepop$country_code, (migcounts[,2:ncol(migcounts)]*200.)/((mergepop[,3:ncols]+mergepop[,2:(ncols-1)])/2.))
}
	
popagesex <- function(sexm, agem, ...){
	age <- agem
	sex <- sexm
	if(is.null(age)) age <- '0-4'
	if(is.null(sex)) sex <- 'F'
	if(length(sex)==0 || length(age)==0) return(NULL)
	tpop <- tpopp <- NULL			
	for(s in sex) {
		dataset.name <- paste('pop',s, sep='')
		if.not.exists.load(dataset.name)
		pop <- sum.by.country.subset.age(wpp.data.env[[dataset.name]], age)
		if(!is.null(tpop)){
			tpop <- cbind(country_code=tpop[,'country_code'], tpop[,2:ncol(tpop)] + pop[,2:ncol(pop)])
		} else tpop<-pop
		if(wpp.year.from.package.name(wpp.data.env$package) > 2010) { #projection stored separately from observations
			dataset.name <- paste('pop', s, 'projMed', sep='')
			if.not.exists.load(dataset.name)
			popp <- sum.by.country.subset.age(wpp.data.env[[dataset.name]], age)
			if(!is.null(tpopp)){
				tpopp <- cbind(country_code=tpopp[,'country_code'], tpopp[,2:ncol(tpopp)] + popp[,2:ncol(popp)])
			} else tpopp<-popp
		}
	}
	if(!is.null(tpopp)) tpop <- merge(tpop, tpopp, by='country_code')
	tpop
}

mortagesex <- function(sex, age, ...){
	if(is.null(age)) age <- '0'
	if(is.null(sex)) sex <- 'F'
	dataset.name <- paste('mx',sex, sep='')
	if.not.exists.load(dataset.name)
	sum.by.country.subset.age(wpp.data.env[[dataset.name]], age)
}

fertage <- function(age, ...){
	if(is.null(age)) age <- '15-19'
	if.not.exists.load('percentASFR')
	tfert <- fert()
	tfert <- cbind(country_code=tfert$country_code, tfert[,.get.year.cols.idx(tfert)])
	asfr <- sum.by.country.subset.age(wpp.data.env[['percentASFR']], age)
	tfert <- tfert[tfert$country_code %in% asfr$country_code,]
	tfert <- tfert[match(asfr$country_code, tfert$country_code), ] # put rows in the same order
	#browser()
	cbind(country_code=tfert[,'country_code'], tfert[,2:ncol(tfert)] * asfr[,2:ncol(asfr)] / 100.)
}

pfertage <- function(agem, ...){
	age <- agem
	if(is.null(age)) age <- '15-19'
	if.not.exists.load('percentASFR')
	sum.by.country.subset.age(wpp.data.env[['percentASFR']], age)
}

fert <- function(...) {
	name.pred <- if(wpp.data.env$package=='wpp2008') NULL else 'tfrprojMed'
	return(load.and.merge.datasets('tfr', name.pred))
}

leF <- function(...) {
	name.pred <- if(wpp.data.env$package=='wpp2008') NULL else 'e0Fproj'
	return(load.and.merge.datasets('e0F', name.pred))
}

leM <- function(...) {
	name.pred <- if(wpp.data.env$package=='wpp2008') NULL else 'e0Mproj'
	return(load.and.merge.datasets('e0M', name.pred))
}

sexratio <- function(...) {
	return(load.and.merge.datasets('sexRatio', NULL))
}

meanagechbear <- function(...) {
	# mean age of child bearing
	data <- load.and.merge.datasets('percentASFR', NULL)
	ddply(data[,-which(colnames(data) == "age")], "country_code", .fun=colwise(function(x) sum(seq(17.5, by=5, length=7)*x)/100.))
}

.sum.popFM.keep.age <- function() {
	name.preds <- if(wpp.year.from.package.name(wpp.data.env$package) <= 2010) c(NULL, NULL) else c('popFprojMed', 'popMprojMed')
	pF <- load.and.merge.datasets('popF', name.preds[1], by=c('country_code', 'age'), remove.cols=c('country', 'name'))
	pM <- load.and.merge.datasets('popM', name.preds[2], by=c('country_code', 'age'), remove.cols=c('country', 'name'))
	cbind(country_code=pF[,1], pF[,-c(1,2)] + pM[,-c(1,2)])
}

medage <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(gmedian))
}

meanageinchbearage <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(gmean.child.bearing))
}

tdratio <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(dependency.ratio, which='total'))
}

psratio <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(function(x) 1/dependency.ratio(x, which='old')))
}

chdratio <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(dependency.ratio, which='child'))
}

oadratio <- function(...) {
	ddply(.sum.popFM.keep.age(), "country_code", .fun=colwise(dependency.ratio, which='old'))
}

popgrowth <- function(...) {
	pop <- tpop()
	ncols <- ncol(pop)
	#browser()
	cbind(country_code=pop$country_code, log(pop[,3:ncols]/pop[,2:(ncols-1)])/5)
}

.pi.suffix <- function(x) c(low='l', high='u')[x]

fert.ci <- function(which.pi, bound, ...) {
	# which.pi is for '80', '95' or 'half.child'
	# bound is 'low' or 'high'
	if(wpp.data.env$package=='wpp2008') return(NULL)
	if(wpp.data.env$package=='wpp2010' && which.pi != 'half.child') return(NULL)
	dataset.name <- if(which.pi == 'half.child') paste0('tfrproj', capitalize(bound))
					else paste0('tfrproj', which.pi, .pi.suffix(bound))
	load.and.merge.datasets(dataset.name, NULL)
}

leF.ci <- function(which.pi, bound, ...) {
	e0.ci('F', which.pi, bound)
}

leM.ci <- function(which.pi, bound, ...) {
	e0.ci('M', which.pi, bound)
}

e0.ci <- function(sex, which.pi, bound) {
	if(wpp.year.from.package.name(wpp.data.env$package) <= 2010 || which.pi == 'half.child') return(NULL)
	load.and.merge.datasets(paste0('e0', sex, 'proj', which.pi, .pi.suffix(bound)), NULL)
}

tpop.ci <- function(which.pi, bound, ...) {
	# which.pi is for '80', '95' or 'half.child'
	# bound is 'low' or 'high'
	if(wpp.year.from.package.name(wpp.data.env$package) <= 2010) return(NULL)
	dataset.name <- if(which.pi == 'half.child') paste0('popproj', capitalize(bound))
					else paste0('popproj', which.pi, .pi.suffix(bound))
	load.and.merge.datasets(dataset.name, NULL)
}

popagesex.ci <- function(which.pi, bound, sexm, agem, ...) {
	# bound is 'low' or 'high'
	if((wpp.year.from.package.name(wpp.data.env$package) <= 2010) || (length(sexm) > 1) || (length(agem) > 1) || (which.pi != 'half.child')) 
		return(NULL)
	dataset.name <- paste('pop', sexm, 'proj', capitalize(bound), sep='')
	if.not.exists.load(dataset.name)
	sum.by.country.subset.age(wpp.data.env[[dataset.name]], agem)
}

load.dataset.and.sum.by.country<-function(dataset){
	if.not.exists.load(dataset)
	pop <- sum.by.country(wpp.data.env[[dataset]])
}

trim.spaces <- function (x) gsub("^\\s+|\\s+$", "", x)

if.not.exists.load <- function(name) {
	if(exists(name, where=wpp.data.env, inherits=FALSE)) return()
	do.call('data', list(name, package=wpp.data.env$package, envir=wpp.data.env))
	# special handling of the age column (mostly because of inconsistent labels in the various datasets)
	# trim spaces in age column if needed
	if('age' %in% colnames(wpp.data.env[[name]]) && is.factor(wpp.data.env[[name]]$age)) {
		levels(wpp.data.env[[name]]$age) <- trim.spaces(levels(wpp.data.env[[name]]$age))
		# 'age' in the mx dataset should be numeric but includes 100+, so it's factor
		# replace by 100 and make it numeric
		levs <- levels(wpp.data.env[[name]]$age)
		if("5" %in% levs && "100+" %in% levs) {
			levels(wpp.data.env[[name]]$age)[levs == "100+"] <- "100"
			wpp.data.env[[name]]$age <- as.integer(as.character(wpp.data.env[[name]]$age))
		}
	}
}

load.and.merge.datasets <- function(name.obs, name.pred=NULL, by='country_code', remove.cols=c('country', 'name')){
	if.not.exists.load(name.obs)
  	data <- wpp.data.env[[name.obs]]
  	if(length(remove.cols) > 0) data <- data[,-which(colnames(data)%in%remove.cols)]
  	if(!is.null(name.pred)){
  		# load predictions
  		if.not.exists.load(name.pred)
  		data.pred <- wpp.data.env[[name.pred]]
  		if(length(remove.cols) > 0) data.pred <- data.pred[,-which(colnames(data.pred)%in%remove.cols)]
  		data <- merge(data, data.pred, by=by, sort=FALSE)
  	}
	data
}

lookupByIndicator <- function(indicator, sex.mult=c(), sex=c(), age.mult=c(), age=c()) {
	indicator <- as.numeric(indicator)
	fun <- ind.fun(indicator)
	# load observed data
	#if(fun == 'mortagesex') browser()
	if(!is.null(wpp.data.env[[fun]])) return(wpp.data.env[[fun]])
	data <- wpp.indicator(fun, sexm=sex.mult, sex=sex, agem=age.mult, age=age)
	if(!ind.is.by.age(indicator))
		wpp.data.env[[fun]] <- data
	data
}

lookupByIndicatorInclArea <- function(indicator, ...) {
	if (as.numeric(indicator) == 0) {
		env <- new.env()
		data('UNlocations', envir=env, package=wpp.data.env$package)
		iso <- wpp.data.env$iso3166
		df <- merge(iso[iso$is.country,c('charcode', 'uncode')], env$UNlocations[,c('country_code','area_name')], 
							by.x='uncode', by.y='country_code')[,-1]
		colnames(df)[2] <- .indicator.title.incl.area(0)
		return(df)
	}
	lookupByIndicator(indicator, ...)
}

lookupByIndicator.mchart <- function(indicator, ...) {
	exdf <- wpp.data.env$mchart.data
	name <- .indicator.title.incl.area(indicator[1], ...)
	iso <- wpp.data.env$iso3166
	iso <- iso[iso$is.country,]
	if(!is.null(exdf) && name %in% colnames(exdf)) {
		exdf <- merge(iso[,c('charcode', 'name')], exdf)
		return(exdf[,-1])
	}
	df <- lookupByIndicatorInclArea(indicator[1])
	colnames(df)[which(colnames(df)=='value')] <- name
	if(!is.null(exdf)) 
		df <- merge(exdf, df)
	if(length(indicator) > 1) {
		for (ind in 2:length(indicator)) {
			name <- .indicator.title.incl.area(indicator[ind], ...)
			if(name %in% colnames(df)) next
			if (as.numeric(indicator[ind]) == 0) {
				env <- new.env()
				data('UNlocations', envir=env, package=wpp.data.env$package)
				locs <- merge(iso[,c('charcode', 'uncode')], env$UNlocations[,c('country_code','area_name')], 
							by.x='uncode', by.y='country_code')[,-1]
				#browser()
				df <- merge(df, locs, by='charcode')
				colnames(df)[which(colnames(df)=='area_name')] <- name
			} else {
				df <- merge(df, lookupByIndicator(indicator[ind], ...))
				colnames(df)[which(colnames(df)=='value')] <- name
			}
		}
	}
	wpp.data.env$mchart.data <- df
	df <- merge(iso[,c('charcode', 'name')], df)
	return(df[,-1])
}

.indicator.title.incl.area <- function(indicator, ...) {
	indicator <- as.numeric(indicator)
	if(indicator == 0) return('UN Areas')
	return(get.indicator.title(indicator, ...))
}

.get.pi.name <- function(x) c('80', '95', 'half.child')[x]
.get.pi.name.for.label <- function(x) c('80%', '95%', '1/2child')[x]

getUncertainty <- function(indicator, which.pi, bound='low', sex.mult=c(), sex=c(), age.mult=c(), age=c()) {
	indicator <- as.numeric(indicator)
	if(!ind.is.low.high(indicator)) return(NULL)
	if(length(which.pi) == 0) return(NULL)
	fun <- paste(ind.fun(indicator), 'ci', sep='.')
	all.data <- NULL
	for(i in 1:length(which.pi)) {
		pi.idx <- as.integer(which.pi[i])
		pi.name <-.get.pi.name(pi.idx)
		lookup.name <- paste(fun, pi.name, bound, sep='.')
		if(!is.null(wpp.data.env[[lookup.name]])) data <- wpp.data.env[[lookup.name]]
		else {
			data <- wpp.indicator(fun, pi.name, bound=bound, sexm=sex.mult, sex=sex, agem=age.mult, age=age)
			if(is.null(data)) next
			if(!ind.is.by.age(indicator))
  				wpp.data.env[[lookup.name]] <- data
  		}
  		colnames(data) <- sub('value', paste0('value.', pi.idx), colnames(data))
  		all.data <- if(is.null(all.data)) data 
  					else merge(all.data, data, by=c('charcode', 'Year'))
  	}
	all.data
}

.get.year.col.names <- function(col.names) {
	col.names <- gsub('.y', '', col.names, fixed=TRUE)
	l <- nchar(col.names)	
	substr(col.names, l-3, l)
}

.get.year.cols.idx <- function(data, remove.duplicate.columns=TRUE) {
	year.cols.idx <- grep('[0-9]{4}$|[0-9]{4}.y$', colnames(data))
	# if(remove.duplicate.columns) {
  		# dupl.year <- duplicated(.get.year.col.names(colnames(data)[year.cols.idx]), fromLast=TRUE)
 		# if(any(dupl.year)) year.cols.idx <- year.cols.idx[-which(dupl.year)]
 	# }
 	year.cols.idx
}

merge.with.un.and.melt <- function(data, id.vars='charcode', what=NULL) {
	year.cols.idx <- .get.year.cols.idx(data)
  	year.cols <- colnames(data)[year.cols.idx]
	data <- merge(wpp.data.env$iso3166[,c('uncode', 'name', 'charcode')], data, 
					by.x='uncode', by.y='country_code', sort=FALSE)
  	data <- data[,-which(colnames(data)=='uncode')] 
  	data <- melt(data,
               id.vars = id.vars, 
               measure.vars = year.cols,
               variable.name = 'Year',
               na.rm=TRUE)
	data$Year <- as.numeric(.get.year.col.names(as.character(data$Year)))
	#if(!is.null(what) && ind.mid.years(what))
	#	data$Year <- data$Year - 2
	#browser()
	data	
}

sum.by.country <- function(dataset) {
	year.cols.idx <- grep('^[0-9]{4}', colnames(dataset))
	ddply(dataset[,c(which(colnames(dataset)=='country_code'), year.cols.idx)], "country_code", .fun=colwise(sum))
}

sumMFbycountry <- function(datasetM, datasetF) {
	tpopM <- sum.by.country(datasetM)
	tpopF <- sum.by.country(datasetF)
	cbind(country_code=tpopM[,'country_code'], tpopM[,2:ncol(tpopM)] + tpopF[,2:ncol(tpopF)])
}

sum.by.country.subset.age <- function(dataset, ages) {
	if('100+' %in% ages) ages <- c(ages, "100")
	sum.by.country(with(dataset, dataset[gsub("^\\s+|\\s+$", "", age) %in% ages,]))
}



preserveStructure <- function(dataFrame) {
  structure(
    lapply(names(dataFrame), function(name) {I(dataFrame[[name]])}),
    names=names(dataFrame)
  )
}

ind.settings <- function() attr(wpp.data.env$indicators, 'settings')
ind.fun <- function(indicator) rownames(ind.settings())[indicator]
ind.is.by.age <- function(indicator) ind.settings()[indicator, 'by.age']
ind.is.low.high <- function(indicator) ind.settings()[indicator, 'low.high']
ind.no.age.sum <- function(indicator) ind.settings()[indicator, 'no.age.sum']
ind.sum.in.table <- function(indicator) ind.settings()[indicator, 'sum.in.table']
ind.mid.years <- function(indicator) ind.settings()[indicator, 'mid.years']
ind.digits <- function(indicator) ind.settings()[indicator, 'digits']
ind.definition <- function(indicator) attr(wpp.data.env$indicators, 'definition')[indicator]

set.data.env <- function(name, value) wpp.data.env[[name]] <- value

gmedian <- function(f, cats=NULL) {
	# group median
	if(is.null(cats)) cats <- seq(0, by=5, length=length(f)+1)
	nhalf <- sum(f)/2.
	cumsumf <- cumsum(f)
	medcat <- findInterval(nhalf, cumsumf) + 1
	med <- cats[medcat] + ((nhalf-cumsumf[medcat-1])/f[medcat])*(cats[medcat+1]-cats[medcat])
	return(med)
}

gmean <- function (f, cats = NULL) 
{
    if (all(is.na(f))) 
        return(NA)
    if (is.null(cats)) 
        cats <- seq(0, by = 5, length = length(f) + 1)
    l <- min(length(cats), length(f) + 1)
    mid.points <- cats[1:(l - 1)] + (cats[2:l] - cats[1:(l - 
        1)])/2
    counts <- f * mid.points
    return(sum(counts)/sum(f))
}



gmean.child.bearing <- function(f) {
	# group mean of child bearing age
	return(gmean(f[4:10], cats=seq(15, by=5, length=8)))
}

dependency.ratio <- function(counts, which='total'){
	nom <- 0
	if(which %in% c('total', 'child')) nom <- nom + sum(counts[1:3])
	if(which %in% c('total', 'old')) nom <- nom + sum(counts[14:21])
	nom/sum(counts[4:13])	
}

get.pyramid.data <- function(year, countries, which.pi=NULL, bound=NULL, indicators=c(F='popF', M='popM'), load.pred=TRUE) {
	name.preds <- name.obs <- c(NULL, NULL)
	if(is.null(which.pi)) {
		name.obs <- indicators
		if(wpp.year.from.package.name(wpp.data.env$package) > 2010 && load.pred) name.preds <- paste(indicators, 'projMed', sep='')
	} else { #PIs
		# only +-half.child available
		if(wpp.year.from.package.name(wpp.data.env$package) > 2010 && 'half.child' %in% .get.pi.name(as.integer(which.pi))) 
			name.obs <- paste(indicators, 'proj', capitalize(bound), sep='')
	}
	if(all(is.null(c(name.preds, name.obs)))) return(NULL)
	dataB <- list()
	for(i in 1:min(2,length(indicators))) {
		p <- load.and.merge.datasets(name.obs[i], name.preds[i], by=c('country_code', 'age'), remove.cols=c('country', 'name'))
		dataB[[i]] <- merge.with.un.and.melt(cbind(p, age.num=.get.age.num(p$age)), id.vars=c('charcode', 'age', 'age.num'),
				what=indicators[i])
		dataB[[i]] <- cbind(dataB[[i]], sex=names(indicators)[i])
	}	
	data <- wpp.by.year(if(length(indicators) > 1) rbind(dataB[[1]], dataB[[2]]) else dataB[[1]], year)
	wpp.by.countries(data, countries)
}

.get.pASFR <- function(year, countries) {
	if.not.exists.load('percentASFR')
	asfr <- wpp.data.env[['percentASFR']]
	asfr <- asfr[,-which(is.element(colnames(asfr), c('country', 'name')))]
	asfrm <- wpp.by.countries(wpp.by.year(
				merge.with.un.and.melt(cbind(asfr, age.num=.get.age.num(asfr$age)), id.vars=c('charcode', 'age', 'age.num')), year), countries)
	asfrm
}

get.age.profile.fert <- function(year, countries){
	asfrm <- .get.pASFR(year, countries)
	#browser()
	tfert <- fert()
	tfert <- cbind(country_code=tfert$country_code, tfert[,.get.year.cols.idx(tfert)])
	tfertm <- wpp.by.countries(wpp.by.year(
				merge.with.un.and.melt(tfert, id.vars='charcode'), year), countries)
	colnames(tfertm)[2] <- 'tfr'
	data <- merge(asfrm, tfertm, by='charcode')
	data <- ddply(data, 'charcode', mutate, value = get("value")/100. * get("tfr"))
	data$tfr <- NULL
	data
}

get.age.profile.pfert <- function(year, countries){
	.get.pASFR(year, countries)
}

get.indicator.title <- function(indicator, sex.mult=c(), sex=c(), age.mult=c(), age=c()) {
	indicator <- as.numeric(indicator)
	title <- names(wpp.data.env$indicators)[indicator]
	if (!ind.is.by.age(indicator)) return(title)
	if(ind.no.age.sum(indicator)){
		sex.string <- paste('sex: ', sex, sep='')
		age.string <- paste('age: ', age, sep='')
	} else { # multiple sex and age groups possible
		sex.string <- paste('sex: ', if(length(sex.mult)>1) 'Both' else sex.mult, sep='')
		age.string <- paste('age: ', paste(age.mult, collapse=', '), sep='')
	}
	return(paste(title, sex.string, age.string, sep='; '))
}

.get.age.num <- function(age) {	
	# Return numeric version of the age, either its index or its numeric value
	aorder <- .get.age.order()
	#browser()
	if(any(!(age %in% names(aorder)))) return(age)
	aorder[as.character(age)]
} 
.get.age.order <- function() {
	age <- c(paste(seq(0, by=5, length=20), seq(4, by=5, length=20), sep='-'), '100+')
	age.array <- 1:21
	names(age.array) <- age
	age.array
}
