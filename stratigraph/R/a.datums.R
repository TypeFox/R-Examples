"a.datums" <-
function(x, depths = NULL, skip = 0,
         increasing.down = FALSE, ...){

if(!is.strat.column(x)){
  stop(paste('argument to a.datums is not of class strat.column'))
}

if(is.null(x$depths) && is.null(depths)){
  stop('appearance datums are meaningless if depths are not supplied; for indices, set depths = 1:ncol(x)')
}

y <- x$counts

nas <- sum(is.na(y))
if(nas > 0){
  warning(paste(nas, 'NAs replaced with zeros for the purposes of calculating stratigraphic ranges; original counts have not been modified'))
  y[is.na(y)] <- 0
}

if(is.null(depths)){
  depths <- x$depths
}

index.pos <- function(x) {return((1:length(x))[x > 0])}

presence.index <- apply(y, 2, index.pos)

if(skip > 0){
  depths[c(1:skip, (length(depths) - skip):length(depths))] <- NA
}

if(increasing.down){
  fads <- depths[as.numeric(lapply(presence.index, max))]
  lads <- depths[as.numeric(lapply(presence.index, min))]
}else{
  lads <- depths[as.numeric(lapply(presence.index, max))]
  fads <- depths[as.numeric(lapply(presence.index, min))]
}

if(is.null(x$taxa)){
  names(fads) <- colnames(x$counts)
  names(lads) <- colnames(x$counts)
}else{
  names(fads) <- x$taxa
  names(lads) <- x$taxa
}

return(data.frame(fads = fads, lads = lads))
}

"fads" <-
function(x, depths = NULL, skip = 0,
         increasing.down = FALSE, ...){

a.datums <- a.datums(x, depths = depths, skip = skip,
                     increasing.down = increasing.down, ...)
fads <- a.datums$fads
names(fads) <- rownames(a.datums)
return(fads)
}

"lads" <-
function(x, depths = NULL, skip = 0,
         increasing.down = FALSE, ...){

a.datums <- a.datums(x, depths = depths, skip = skip,
                     increasing.down = increasing.down, ...)
lads <- a.datums$lads
names(lads) <- rownames(a.datums)
return(lads)

}