tv.site <- function (db, tv_home, drop = TRUE, common.only = FALSE, iconv="CP437", verbose = TRUE, replace.names, ...) 
{
#  ow <- options('warn')
# if(quiet) { options(warn=-1) }
# if (is.list(db)) site <- tv.mysql(db, "tvhabita") else {
  if (missing(tv_home)) tv_home <- tv.home()
  site <- read.dbf(file.path(tv_home, "Data", db[1], "tvhabita.dbf"), as.is=TRUE)
  if(!missing(replace.names)) 
    for(i in 1:nrow(replace.names)) 
      names(site) <- sub(paste('^', replace.names[i,1], '$', sep=''), replace.names[i,2], names(site))
  if (suppressWarnings(any(site[,sapply(site, is.numeric)] < -1e+05, na.rm = TRUE))) 
    message(paste("WARNING! Values less than -100,000. \n", "WARNING! tvhabita.dbf may be corrupt. \n", "WARNING! Please correct by im- / exporting e.g. with OpenOffice."))
  if(length(db)>1) 
    for(i in 2:length(db)) {
	    site.tmp <- read.dbf(file.path(tv_home, 'Data', db[i],'tvhabita.dbf'))
  	if(any(site$RELEVE_NR %in% site.tmp$RELEVE_NR)) stop('Found duplicate releve numbers in ', db[i] , ' aborting!')
  	if(!missing(replace.names)) 
  	  for(i in 1:nrow(replace.names)) 
  	    names(site.tmp) <- sub(paste('^', replace.names[i,1], '$', sep=''), replace.names[i,2], names(site.tmp))
  	cols1 <- names(site)
  	cols2 <- names(site.tmp)
  	if (common.only){
	  	common <- intersect(cols1, cols2)
		  site <- rbind(site[, common], site.tmp[, common])
	  } else {
 		  All <- union(cols1, cols2)
  		miss1 <- setdiff(All, cols1)
  		miss2 <- setdiff(All, cols2)
  		site[, c(as.character(miss1))] <- NA
  		site.tmp[,c(as.character(miss2))] <- NA
  		site <- rbind(site, site.tmp)
  	} 
  }

### Conversion of factors
# fac <- sapply(site, is.factor)
    for(i in names(site)) if(is.character(site[,i])) site[,i] <- iconv(site[,i], iconv, '')

    ### Time
    if(any(is.na(site$DATE))) 
      message(sum(is.na(site$DATE)), ' releves without date. Not converted from factor to date format.') else {
    site$DATE <- gsub('/','',site$DATE)
#   Date <- rep('no date', nrow(site))
    index <- nchar(as.character(site$DATE))==4
    fun <- function(x) paste(x,'0101',sep='')
    site$DATE[index] <- fun(site$DATE[index])
#   Date[!index] <- as.character(site$DATE[!index])
    index <- nchar(as.character(site$DATE))==6
    fun <- function(x) paste(x,'01',sep='')
    site$DATE[index] <- fun(site$DATE[index])
    site$DATE <- as.Date(site$DATE, '%Y%m%d')
    }
  ### Survey Area
  n <- sum(site$SURF_AREA == 0 | is.na(site$SURF_AREA))
  if(n>0) message(paste(n, ' releves without survey area'))
  site$SURF_AREA[site$SURF_AREA==0] <- NA

### 
  fun <- function(x) all(is.na(x))
  na <- apply(site, 2, fun)
  if (drop & any(na)) {
      if(verbose) {
        message('Some columns contain no data and are omitted.')
        print(names(site)[na], quote = FALSE)
      }
      site <- site[, !na]
  }
  fun.2 <- function(x) all(x == 0 | is.na(x))
  leer <- apply(site, 2, fun.2)
  if (drop & any(leer)) {
      if(verbose) {
        message('Some numeric columns contain only 0 values and are omitted.')
  	    print(names(site)[leer], quote = FALSE)
      }
	    site <- site[, !leer]
  }
  fun.3 <- function(x) is.numeric(x) & any(x == 0, na.rm = TRUE)
  null <- logical()
  for (i in 1:length(site)) null[i] <- fun.3(site[, i])
  if (any(null) && getOption('warn') >= 0) {
    if(verbose) {
#        options(warn=1)
        message("Some numeric fields contain 0 values:")
        print(names(site)[null], quote = FALSE)
        message('Please check if these are really meant as 0 or if they are erroneously assigned because of DBase restrictions.')
        message(paste("If so, use something like:"))
        message("site$Column_name[site$Column_name==0] <- NA", '\n')
    }
#       paste("summary(site[,c('", paste(names(site)[null], collapse = "','"), "')]) \n", sep = "")
      }
  site <- site[order(site$RELEVE_NR),]
#  warning(ow)
  site
}

