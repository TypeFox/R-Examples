intersectFindFn <- function(e1, e2, sortby=NULL) {
##
## 1.  rbind
##
  e12Names <- intersect(names(e1), names(e2))
  k12 <- length(e12Names)
  if(k12 != length(names(e1)))
    warning('Columns in e1 not in e2;  using only shared columns.')
  if(k12 != length(names(e2)))
    warning('Columns in e2 not in e1;  using only shared columns.')
#
  xy <- rbind(e1[, e12Names], e2[, e12Names])
##
## 2.  Find and merge duplicates
##
  if(!('Link' %in% names(xy)))
    stop('Either x or y does not contain Link')
  pacFn <- paste(xy$Package, xy$Function, sep='\n')
  pF. <- table(pacFn)
  pF <- names(pF.[pF.>1])
  npF <- length(pF)
  if('Description' %in% e12Names)
    xy$Description <- as.character(xy$Description)
  xy$Link <- as.character(xy$Link)
#
  xy. <- xy[rep(1, npF), ]
  for(i in seq(length=npF)){
    idup <- which(pacFn == pF[i])
    xy.[i,] <- xy[idup[1], ]
    if('Count' %in% e12Names)
      xy.[i, 'Count'] <- max(xy[idup, 'Count'])
    if('MaxScore' %in% e12Names)
      xy.[i, 'MaxScore'] <- max(xy[idup, 'MaxScore'])
    if('TotalScore' %in% e12Names)
      xy.[i, 'TotalScore'] <- max(xy[idup, 'TotalScore'])
    if('Score' %in% e12Names)
      xy.[i, 'Score'] <- max(xy[idup, 'Score'])
    if('Description' %in% e12Names){
      nchD <- nchar(xy[idup, 'Description'])
      selD <- which.min(nchD)
      xy.[i, 'Description'] <- xy[idup[selD], 'Description']
    }
#
    nchL <- nchar(xy[idup, 'Link'])
    selL <- which.min(nchL)
    xy.[i, 'Link'] <- xy[idup[selL], 'Link']
  }
##
## 3.  Rebuild summary and resort
##
  xys <- sortFindFn(xy.[,
     c('Package', 'Score', 'Function', 'Date', 'Description', 'Link')],
                          sortby)
##
## 4.  Fix attributes
##
  attr(xys, 'matches') <- c(attr(e1, 'matches'), attr(e2, 'matches'))
  attr(xys, 'string') <- paste(attr(e1, 'string'), attr(e2, 'string'),
                               sep=' & ')
  attr(xys, 'call') <- call( "(", call( "&",
                    attr(e1, "call"), attr(e2, "call") ) )
##
## 5.  Done
##
  xys
}
