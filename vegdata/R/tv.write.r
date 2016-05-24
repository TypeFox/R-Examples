if(getRversion() >= "2.15.1")  utils::globalVariables(c("DoWritedbf"))

tv.write <- function(x, site, db, name, cover=c('code','perc'), overwrite = FALSE, iconv="CP437", newTvAdmin = FALSE, ...) {
  if(missing(db)) warning('Name of original database is necessary to process TvAdmin.dbf, remarks.dbf and tvwin.set')
  cover <- match.arg(cover)
  if('veg' %in% class(x)) {
    X <- reShape.veg(x, ...)
    if(cover == 'perc') {
      X$COVER_CODE <- as.character(X$COVER_PERC)
      site$COVERSCALE <- '00'
      X <- X[,c('RELEVE_NR','TaxonUsageID', 'COVER_CODE', 'LAYER')]
    }
  } else {
      if(!any(c('tv.obs','vw.obs') %in% class(x))) 
        stop('Species observations must be of either \"tv.obs\" or \"vw.obs\" class.')
  if(!all(c('RELEVE_NR') %in% names(x))) 
    stop('column names of species observations must contain RELEVE_NR')
  names(x)[2] <- 'SPECIES_NR'
    X <- x
  }
#  names(X)[2] <- 'SPECIES_NR'
  if(!overwrite) 
    if(file.exists(file.path(options('tv_home'), 'Data', name, 'tvhabita.dbf'))) 
      stop('Database ', name, ' already exists.')
  site$DATE <- gsub('-','',site$DATE)
  for(i in names(site)) if(is.character(site[,i])) {
    site[is.na(site[,i]),i] <- ''
    site[,i] <- iconv(site[,i], '', iconv)
  }
### Add obligigatory fields
# oblig <- c('RELEVE_NR','COUNTRY','REFERENCE','TABLE_NR','NR_IN_TAB','COVERSCALE','PROJECT','AUTHOR','DATE','SYNTAXON','SURF_AREA','UTM','ALTITUDE','EXPOSITION','INCLINATIO','COV_TOTAL','COV_TREES','COV_SHRUBS','COV_HERBS','COV_MOSSES','COV_LICHEN','COV_ALGAE','COV_LITTER','COV_WATER','COV_ROCK','TREE_HIGH','TREE_LOW','SHRUB_HIGH','SHRUB_LOW','HERB_HIGH','HERB_LOW','HERB_MAX','CRYPT_HIGH','MOSS_IDENT','LICH_IDENT','REMARKS')
oblig <- c('RELEVE_NR','COUNTRY','REFERENCE','TABLE_NR','NR_IN_TAB','COVERSCALE','DATE','SURF_AREA','UTM','ALTITUDE','EXPOSITION','INCLINATIO','COV_TOTAL')
for(m in oblig[!oblig %in% names(site)]) {
  site[,m] <- ''
}

### Write
  dir.create(file.path(options('tv_home'), 'Data', name), showWarnings = if(overwrite) FALSE else TRUE)
  write.dbf(site, file.path(options('tv_home'), 'Data', name, 'tvhabita.dbf'))
  write.dbf(X, file.path(options('tv_home'), 'Data', name, 'tvabund.dbf'))
 # Remarks
 if(!missing(db)) {
  # print(file.path(options('tv_home'), 'Data', db[1], 'remarks.dbf'))
  remarks <- read.dbf(file.path(options('tv_home'), 'Data', db[1], 'remarks.dbf'), as.is=TRUE)
  if(length(db) > 1) 
    for(n in 2:length(db)) remarks <- rbind(remarks, read.dbf(file.path(options('tv_home'), 'Data', db[n], 'remarks.dbf'), as.is=TRUE))
  remarks <- remarks[remarks$RELEVE_NR %in% site$RELEVE_NR,]
  op <- options('warn')
  options(warn=-1)
  suppressWarnings(write.dbf(remarks, file.path(options('tv_home'), 'Data', name, 'remarks.dbf')))
  options(op)
 # TvAdmin
  TvAdmin <- read.dbf(file.path(options('tv_home'), 'Data', db[1], 'TvAdmin.dbf'), as.is = TRUE)
  if(length(db) > 1) 
    for(n in 2:length(db)) TvAdmin <- rbind(TvAdmin, read.dbf(file.path(options('tv_home'), 'Data', db[n], 'TvAdmin.dbf')))
  TvAdmin <- TvAdmin[TvAdmin$RELEVE_NR %in% site$RELEVE_NR,]
  TvAdmin$MOD_USER[is.na(TvAdmin$MOD_USER)] <- Sys.getenv('USER')
  TvAdmin$MOD_DATE[is.na(TvAdmin$MOD_DATE)] <- format(Sys.Date())
  write.dbf(TvAdmin, file.path(options('tv_home'), 'Data', name, 'TvAdmin.dbf'))

  if(!any(db == name)) file.copy(from = file.path(options('tv_home'), 'Data', db[1], 'tvwin.set'), to=file.path(options('tv_home'),'Data', name, 'tvwin.set'), overwrite = overwrite)
 }
  cat('Turboveg database', name, 'written to', file.path(options('tv_home'), 'Data', name),'\n')

if(newTvAdmin) {
  cat('If you want to create a new TvAdmin.dbf, please install a library with uuid capabilities (e.g. UUIDgenerate from package uuid, or uuid.gen from dplR):\n')
  cat("tvadmin <- data.frame(RELEVE_NR=site$RELEVE_NR, SOURCE_DB=db,  GUID=replicate(nrow(site), UUIDgenerate()),	CREAT_USER=Sys.getenv('USER'), CREAT_DATE=Sys.Date(), MOD_USER=Sys.getenv('USER'),	MOD_DATE=Sys.Date(), NDFF_QUAL=as.integer(0))\n")
  cat("write.dbf(tvadmin, file.path(options('tv_home'),'Data', db, 'TvAdmin.dbf'))")
}
}

