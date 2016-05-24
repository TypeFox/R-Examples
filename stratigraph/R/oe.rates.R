"oe.rates" <-
function(x, depths = NULL, breaks = NULL,
         per.capita = FALSE, remove.below = 1, ...){

# warn in x is not a strat.column and treat it like a count matrix
if(!is.strat.column(x)){
  warning(paste('argument to oe.rates() is not of class strat.column'))
}

if(is.null(x$counts)) counts <- x
else counts <- x$counts

#clobber x$depths with depths; if both null, use integers
if(is.null(depths)){
  depths <- x$depths
  if(is.null(depths)) depths <- 1:nrow(counts)
}
if(length(depths) != nrow(counts)){
  stop('depth not the same lengths as rows in x')
}

# check for missing values
nas <- sum(is.na(counts))
if(nas > 0){
  warning(paste(nas, 'NAs replaced with zeros for the purposes of calculating stratigraphic ranges; original counts have not been modified'))
  counts[is.na(counts)] <- 0
}

# get breaks from depths if they are not specified
depth.plus.lagged <- cbind(depths[-length(depths)], depths[-1])
if(is.null(breaks)) breaks <- apply(depth.plus.lagged, 1, mean)

# names for intervals based on breaks
paster <- function(x){paste(x, collapse = '/')}
int.names <- cbind(breaks[-length(breaks)], breaks[-1])
int.names <- apply(int.names, 1, paster)
#int.names <- c('above', int.names, 'below')
ints <- diff(breaks)

# names for boundaries based on sample locations
bound.names <- apply(depth.plus.lagged, 1, paster)

# get numeric presence/absence matrix with ranges filled
pa <- fill.ranges(strat.column(counts = counts, depths = depths),
                  out = 'pa')$counts + 0

# remove species that appear in remove.below or fewer rows
pa <- pa[,colSums(pa) > remove.below]
removed <- ncol(counts) - ncol(pa)
cat(paste(ncol(pa), 'species appear at more than',
          remove.below, 'level(s);', removed, 'removed\n'))

# calculate standing diversity and o/e rates at sampled levels
stand.div <- rowSums(pa)
oe <- a.datums(strat.column(counts, depths = depths),
               increasing.down = TRUE)
oe.merge <- merge(as.data.frame(table(oe$fads)),
            as.data.frame(table(oe$lads)),
            by = 1, all = TRUE)
oe.merge <- merge(as.data.frame(depths), oe.merge,
            by = 1, all = TRUE)
row.names(oe.merge) <- oe.merge$depths
oe.merge <- oe.merge[,-1]
names(oe.merge) <- c('orig', 'ext')
oe.merge[is.na(oe.merge)] <- 0
basic <- cbind(stand.div, oe.merge)
if(per.capita){
  basic$orig <- basic$orig/basic$stand.div
  basic$ext <- basic$ext/basic$stand.div
}

# count boundary crossers
bx <- rep(0, length(bound.names))
for(i in 1:(length(bound.names) - 1)){ #loop through boundaries
  for(j in 1:ncol(pa)){ #loop through spp.
    if(pa[i,j] > 0 && pa[i+1,j] > 0) bx[i] <- bx[i] + 1
  }
}
names(bx) <- bound.names

# Fundamental measures from Foote (2000b):
#  FL = # of taxa confined to interval (not crossing either
#       top or bottom)
#  bL = # of taxa going extinct in interval (crossing only
#       bottom boundary)
#  Ft = # of taxa originating in interval (crossing only top
#       boundary)
#  bt = # of taxa crossing both top and bottom of interval

# calculate standing diversity and o/e rates at breaks
fund <- matrix(NA, nrow = length(breaks) - 1, ncol = 4)
colnames(fund) <- c('FL', 'bL', 'Ft', 'bt')
for(i in 1:(length(breaks) - 1)){ #loop through intervals
  FL <- bL <- Ft <- bt <- 0
  for(j in 1:ncol(pa)){ #loop through species
    fad <- oe$fads[rownames(oe) == colnames(pa)[j]]
    lad <- oe$lads[rownames(oe) == colnames(pa)[j]]
    top <- breaks[i]
    bottom <- breaks[i+1]
    if(is.na(top)) top <- -Inf
    if(is.na(bottom)) bottom <- Inf
    #cat(colnames(pa)[j], '\t')
    #cat(fad, lad, bottom, top, '\t')
    if(!(fad < top || lad > bottom)){ # taxon in interval
        #cat('overlap\t')
      if(fad < bottom && lad > top){
        FL <- FL + 1 #confined to interval
        #cat('FL\n')
      }else if(fad >= bottom && lad >= top){
        bL <- bL + 1 #only bottom crossed
        #cat('bL\n')
      }else if (fad <= bottom && lad <= top){
        Ft <- Ft + 1 #only top crossed
        #cat('FT\n')
      }else if(fad >= bottom && lad <= top){
        bt <- bt + 1 #crossing both boundaries
        #cat('bt\n')
      }else{
        stop('huh?')
      }
    }#else cat ('\n')
  }
  fund[i,] <- c(FL, bL, Ft, bt)
}

rownames(fund) <- int.names

botx <- fund[,'bL'] + fund[,'bt']
topx <- fund[,'Ft'] + fund[,'bt']
stand.est <- (botx + topx) / 2

#  negative natural Log Origination Per Time (lopt)
lopt <- -log(fund[,'bt']/topx) / ints
#  negative natural Log EXtinction per Time (lext)
lext <- -log(fund[,'bt']/botx) / ints


interpol <- cbind(fund, botx, topx, stand.est, lopt, lext)

return(list(basic = basic, bx = bx,
            intervals = ints, interpol = interpol))
} # End of function