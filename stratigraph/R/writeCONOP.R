"writeCONOP" <-
function(x, tax.check = TRUE, out.dir = getwd(),
         sect.names = NULL, neg.depths = FALSE,
         scaling = rep ('yes', length(x)),
         REGION = NULL, AGE = NULL,
         REF1 = NULL, REF2 = NULL, REF3 = NULL){

if(!all(lapply(x, is.strat.column)))
  stop('x is not a list of objects of type strat.column')

if(is.null(names(x))){
  if(is.null(sect.names)){
    names(x) <- paste('Section', 1:length(x), sep = '')
    warning('section names have not been provided; using Section1, ..., SectionN')
  }else{
    names(x) <- sect.names
  }
}

# Put together a list of all the taxon codes in all sections (plus other data?)
#regions.vect <- vector(length = length(x), mode = 'character')
#ages.vect <- vector(length = length(x), mode = 'character')
for(i in 1:length(x)){
  if (is.null(x[[i]]$tax.code))
    stop(paste('section', i, 'is missing a taxon code vector (tax.code) with which to create a taxon dictionary'))
  if(i == 1) taxon.vect <- x[[1]]$tax.code
  else taxon.vect <- c(taxon.vect, x[[i]]$tax.code)
  #regions.vect[i] <- paste(x[[i]]$country, x[[i]]$state)
  #ages.vect[i] <- x[[i]]$age
}

# Put together the section dictionary table
section.vect <- names(x)
section.dic <- matrix('*', ncol = 9, nrow = length(section.vect))
colnames(section.dic) <- c('NAME', 'ABBREV', 'FILE.NAME',
                           'REGION', 'AGE', 'REF1', 'REF2',
                           'REF3', 'SCALING')
section.dic <- as.data.frame(section.dic)

if(!is.null(REGION)) region.vect <- REGION
if(!is.null(AGE)) region.vect <- AGE

section.dic$NAME <- names(x)
#section.dic$ABBREV <- substr(1, 3, names(x)) #UNIQUIZE!
section.dic$FILE.NAME <- names(x)
#section.dic$REGION <- regions.vect
#section.dic$AGE <- ages.vect
if(!is.null(REF1)) section.dic$REF1 <- REF1
if(!is.null(REF2)) section.dic$REF2 <- REF2
if(!is.null(REF3)) section.dic$REF3 <- REF3
section.dic$SCALING <- scaling

# if tax.check is true or numeric, fork to different return
if(tax.check == TRUE && is.logical(tax.check)){
  nt <- length(unique(taxon.vect))
  warning(paste('there are', nt, 'unique taxa in the dictionary; check that there are no synonymies and then rerun with tax.check != TRUE'))
}

#  tax.check.frame <- matrix('', nrow = nt, ncol = 11)
#  names(tax.check.frame) <- c('name', 'short.name',
#                              'long.name', 'family',
#                              'genus', 'species',
#                              'subspecies', 'suffix',
#                              'morphotype', 'auth', 'date')  
#  tax.check.frame$short.name <- taxon.short.vect
#  tax.check.frame$long.name <- taxon.vect
#  if(!is.null(taxon.short.vect)){
#  	tax.check.frame$name <- taxon.short.vect
#  }else{
#    tax.check.frame$name <- taxon.vect
#  }
#  # assign other taxon data to tax.check.frame
#  ### STUFF
  
#  # Remove unique rows in tax.check.frame
#  ### STUFF

#  if(is.null(out.dir)){
#    return(tax.check.frame)
#  }else{
#    write.table(tax.check.frame,
#                paste(out.dir, 'tax.check.frame', sep = ''))
#  }
#}else if (is.numeric(tax.check)){
#  #sort out the duplicate taxa
#  ### STUFF
#}

# Put together the taxon dictionary table
taxon.vect <- sort(unique(taxon.vect))
taxon.dic <- matrix('*', ncol = 11, nrow = length(taxon.vect))
colnames(taxon.dic) <- c('REF.NUM', 'TAX.ABBREV', 'TAX.CAT',
                         'GENUS', 'SP', 'SUBSP', 'MORPHOTYPE',
                         'AUTH', 'DATE', 'TAG', 'SYN.NUM')
taxon.dic <- as.data.frame(taxon.dic)

taxon.dic$REF.NUM <- paste('0', 99 + (1:length(taxon.vect)), sep = '')
taxon.dic$TAX.ABBREV <- taxon.vect
#taxon.dic$TAX.CAT <- 
taxon.dic$GENUS <- rep('Unknown', length(taxon.vect))
taxon.dic$SP <- rep('unknown', length(taxon.vect))
#taxon.dic$SUBSP <- 
#taxon.dic$MORPHOTYPE <- 
#taxon.dic$AUTH <- 
#taxon.dic$DATE <- 
#taxon.dic$TAG <- 
#taxon.dic$SYM.NUM <- tax.check

# Produce the individual section data files
sects <- vector(length = length(x), mode = 'list')
for(i in 1:length(x)){
  cts <- x[[i]]$counts
  nms <- x[[i]]$tax.code
  
  sects[[i]] <- matrix('*', ncol = 8, nrow = ncol(cts))
  colnames(sects[[i]]) <- c('SORTCODE', 'TAXCODE', 'FAD',
                           'LAD', 'SORTCODE2', 'TAXCODE2', 'REF',
                           'DATE')
  sects[[i]] <- as.data.frame(sects[[i]])
  
  sects[[i]]$TAXCODE <- nms
  sects[[i]]$TAXCODE2 <- nms
  sects[[i]]$SORTCODE <- taxon.dic$REF.NUM[match(nms, taxon.dic$TAX.ABBREV)]
  sects[[i]]$SORTCODE <- as.character(sects[[i]]$SORTCODE)
  sects[[i]]$SORTCODE2 <- sects[[i]]$SORTCODE
  if(neg.depths == TRUE){
    sects[[i]]$FAD <- format(fads(x[[i]]), nsmall = 2)
    sects[[i]]$LAD <- format(lads(x[[i]]), nsmall = 2)
  }else{
    sects[[i]]$LAD <- format(fads(x[[i]]), nsmall = 2)
    sects[[i]]$FAD <- format(lads(x[[i]]), nsmall = 2)  
  }
#  sects[[i]]$REF <- 
  sects[[i]]$DATE <- rep(toupper(format(Sys.Date(), '%d-%b-%y')), ncol(cts))

  if(!is.null(out.dir)){ 
    write.table(sects[[i]], paste(out.dir, '/', names(x)[i], sep = ''),
                quote = c(2,6), row.names = FALSE, col.names = FALSE)
  }
}

# Return or write files to disk
if(is.null(out.dir)){
  return(list(section.dic = section.dic, taxon.dic = taxon.dic,
              sects = sects))
}else{
  write.table(section.dic, paste(out.dir, '/',
                                 'section.dic', sep = ''),
              quote = TRUE, row.names = FALSE,
              col.names = FALSE)
  write.table(taxon.dic, paste(out.dir, '/',
                               'taxon.dic', sep = ''),
              quote = c(2,8), row.names = FALSE,
              col.names = FALSE)
}

} # End of function

#writeCONOP(x, tax.check = FALSE)

#data(sI)
#sec1 <- strat.column(sI[,-1], depths = sI[,1])
#data(sII)
#sec2 <- strat.column(sII[,-1], depths = sII[,1])
#x <- list(s1 = sec1, s2 = sec2)

#writeCONOP(sects)

# in strat.column(), need to check for presence of names....

# add to strat.column:
#  $age
#  $genus
#  $species
#  $family
#  $tax.code
#  $tax.prefix
#  $tax.auth
#  $tax.date