tv.biblio <- function(x='all', db, dict = tv.dict(db), quiet=FALSE, tv_home, iconv = "CP437", ...) {
  if(missing(tv_home)) tv_home <- tv.home()
  if(missing(db) & missing(dict)) {
    message('Using tvrefenc.dbf from default dictionary.')
    dict = ''
  }
  bibliopath <- file.path(tv_home, 'Popup', dict, 'tvrefenc.dbf')
  biblio <- read.dbf(bibliopath, as.is=TRUE)
  for(i in c('AUTHOR','TITLE','PUBLISHED')) biblio[,i] <- iconv(biblio[,i], iconv, "")  
  if(x[1] != 'all') {
    x <- as.numeric(unique(x))
      biblio <- biblio[match(x, as.numeric(biblio$REFERENCE)),]
      if(!quiet) print(biblio)
  }
  invisible(biblio)
}
