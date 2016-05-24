#### function indicate site conditions
# isc(veg, trait.db, ivnames, keyname = 'LETTERCODE', method = c('mode', 'mean'), db, ...)

isc <- function(veg,
                trait.db, 
                ivnames,
                keyname = 'LETTERCODE',
                method = c('mean', 'mode'),
                weight,
                db,
                ...
) {
  species <- tax('all', quiet = TRUE)
  if(is.character(trait.db)) iv <- tv.traits(trait.db = trait.db, refl = 'GermanSL 1.3', ...) else iv = trait.db
  if(missing(veg)) veg <- tv.veg(db, ...)
  if(!all(ivnames %in% names(iv))) stop('Not all ivnames in table of indicators.')
  if(missing(weight)) {
    iv$weight <- 1
  #  weight <- 'weight'
  } else names(iv)[names(iv) == weight] <- 'weight'
  iv <- as.data.frame(cbind(iv[, match(ivnames, names(iv))], iv[, keyname], iv[, 'weight']))
  names(iv) <- c(ivnames, keyname, 'weight')
#  iv$weight <- as.numeric(iv$weight)
# head(iv1)
  #  print(names(iv))
  # workaround
  if(length(ivnames) == 1) {
    colnames(iv)[1] <- as.character(ivnames)
    iv[,1] <- as.numeric(as.character(iv[,1]))
  }
  # ivnames <- factor(ivnames, levels = ivnames, ordered = TRUE)
#  print(names(iv))
  v <- as.matrix(iv[match(names(veg), iv[, keyname]), ivnames]) #
  rownames(v) <- names(veg)
  if(length(ivnames) == 1) {
    veg <- veg[,!is.na(v)]
    v <- as.matrix(v[!is.na(v),])
  } else  v[is.na(v)] <- 0
  # Species * indicator Matrix of available Species
  w <- as.character(iv$weight[match(names(veg), iv[, keyname])])
  w[is.na(w)] <- "1"
  w <- as.numeric(w)
  veg <- t(t(veg) * w)
  io <- matrix(0, nrow = nrow(veg), ncol = ncol(v)) # Plots * sum of WS indicators for WS
  io <- apply(v, 2, function(x) rowSums(as.matrix(veg/apply(veg, 1, sum)) %*% x, na.rm=TRUE) )
  # io <- apply(v, 2, function(x) rowSums(as.matrix(veg) %*% x, na.rm=TRUE) )
  rS <- rowSums(io, na.rm = TRUE)
  if(any(rS == 0)) {
    cat('The following plots are without a single indicator species:\n')
    print(rownames(veg)[rS == 0])
  }
  
    # Method == max
    if(method == 'mode') {
      IV <- vector('character', nrow(veg))
      for(p in 1:nrow(io)) IV[p] <- paste(ivnames[if(all(io[p,] == 0)) 0 else which(io[p,] == max(io[p,]))], collapse='/')
      # Code to sort automatically according to the order of columns (ivnames)
      IV[IV == ''] <- '/'
      IV = factor(IV, levels(factor(IV))[order(ivnames[match(sapply(strsplit(levels(factor(IV)), '/'), '[[', 1), ivnames)])])
      levels(IV)[levels(IV) == '/'] <- ''
    }
    # Method == mean
    if(method == 'mean') {
      # nis <- apply(veg, 1, function(x) sum((x * v) > 0))
      IV <- io
    }
  names(IV) <- rownames(veg)
    return(IV)
}
### end of function

showplot <- function(veg, plotids)  for(i in 1:length(plotids)) print(veg[plotids[i], veg[plotids[i],]>0])
# showplot(veg, c("361", "362", "363"))

showindiplot <- function(veg, trait.db, plotid, weight, keyname = 'LETTERCODE') {
  if(length(plotid) > 1) {
    warning('more than one plot selected. using only the first.')
    plotid <- plotid[1]
  }
  if(missing(weight)) { trait.db$weight <- 1 } else names(trait.db)[names(trait.db) == weight] <- 'weight'
  pl <- veg[plotid, veg[plotid,]>0] * trait.db[,'weight']
  indi <- trait.db[match(names(pl), trait.db[, keyname]) , 3:8]
  rownames(indi) <- colnames(pl)
  
  indsum <- matrix(unlist(sapply(indi, function(x) x * pl)), nrow=nrow(indi))
  rbind(sapply(indi, function(x) x * pl), SUM = colSums(indsum, na.rm = TRUE))
}
# showindiplot(veg, wsingo, which(ingo == '5+/4+/3+'))


