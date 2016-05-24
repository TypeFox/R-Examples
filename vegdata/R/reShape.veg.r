reShape.veg <- function (veg, crop=TRUE, refl) {
if(!'veg' %in% class(veg)) stop('Only applicable for objects of class \"veg\".')
if(is.null(attr(veg, 'taxreflist')) & missing(refl)) stop('Set option refl because attribute \"taxreflist\" is not set for object \"veg\".')
veg <- as.matrix(veg)
perf <- as.vector(veg)
plots <- as.integer(as.character(dimnames(veg)[[1]][row(veg)]))
spec <- dimnames(veg)[[2]][col(veg)]


  spcnames <- sapply(spec, function(x) strsplit(as.character(x), '.', fixed=TRUE))
  layer <- factor(unlist(lapply(spcnames, function(x) x[2])))
  levels(layer) <- c(levels(layer),'0')
  layer <- as.integer(layer)
  layer[is.na(layer)] <- 0

  df <- data.frame(RELEVE_NR=plots, LETTERCODE=spec, COVER_PERC=perf, LAYER=layer)
  spc <- unlist(lapply(spcnames, function(x) x[1]))
  df$TaxonUsageID <- integer(nrow(df))
  taxa <- tax(spc, syn=FALSE)
  df$TaxonUsageID <- taxa$TaxonUsageID[match(spc, taxa$LETTERCODE)]

  # o <- df[order(df[,1], df[,2]),]
  # df <- df[o,]
  class(df)<- c('tv.obs', 'data.frame')
  if(!is.null(attr(veg, 'taxreflist'))) attr(df, 'taxreflist') <- attr(veg, 'taxreflist')
  if(crop) df <- df[df$COVER_PERC!=0,]
  return(df)
}
