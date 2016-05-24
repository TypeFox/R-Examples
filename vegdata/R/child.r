if(getRversion() >= "2.15.1")  utils::globalVariables(c("gwindow", "gtree", "addHandlerDoubleclick", "svalue"))

##### Child taxa of a taxon
child <- function (x, refl = tv.refl(), gen=4, tree=FALSE, quiet=FALSE, syn=FALSE, ...) {
  if(missing(refl)) refl <- tv.refl()
  species <- tax("all", detailed = TRUE, refl = refl, syn = TRUE, quiet =TRUE, ...)
  if(length(x) > 1) { warning('More than one species selected, using only the first.');  x <- x[1]}
  s <- tax(x, refl = refl, strict = TRUE, quiet = TRUE, ...)
  x <- s$TaxonConceptID
# if(is.character(x) & nchar(x) != 36) stop('x must be given as Taxon ID (GUID or integer).')
  if(tree) {
    root <- child(x, gen=1, ...)
    if(!is.null(root)) {
      offspring <- function(path, ...) {
        ll <- root
        for(i in path)
          ll <- childs(i, gen=1, tree=FALSE, quiet=TRUE, syn=syn)
        off <- logical(nrow(ll))
        for(n in 1:nrow(ll)) off[n] <- !is.null(childs(ll$TaxonUsageID[n], quiet=TRUE))
        if(syn) out <- data.frame(
          Name=ll$TaxonName,
          hasOffspring=off,
          Rang=ll$TaxonRank,
          Synonym=ll$SYNONYM,
          stringsAsFactors=FALSE
        ) else
          out <- data.frame(
            Name=ll$TaxonName,
            hasOffspring=off,
            Rang=ll$TaxonRank,
            Nr=ll$TaxonUsageID,
            stringsAsFactors=FALSE
          )
        out
      }
      if (requireNamespace("gWidgets", quietly = TRUE)) {
      w <- gWidgets::gwindow(paste("Taxonomic Tree of", species$TaxonName[species$TaxonUsageID==x]))
      tr <- gWidgets::gtree(offspring=offspring, container=w)  
      gWidgets::addHandlerDoubleclick(tr, handler=function(h,...) {
        print(childs(gWidgets::svalue(h$obj), gen=1, syn=syn , quiet=TRUE)[, c('TaxonUsageID' , 'LETTERCODE' , 'TaxonName' , 'GRUPPE' , 'TaxonRank' , 'SYNONYM', 'IsChildTaxonOfID' , 'AccordingTo' , 'EDITSTATUS')], row.names=FALSE)
      })
      } else   stop('Please install gWidgetstcltk to run this function with tree = TRUE.')
  }} else {
      x <- species[match(x, species$TaxonUsageID),'TaxonConceptID']
      x <- species[match(x, species$TaxonUsageID),]
      if(syn) {
        ch <- species[which(species$IsChildTaxonOfID == x$TaxonUsageID),'TaxonUsageID']
        if(length(ch)>0) ch <- do.call(rbind, lapply(ch, function(x) syn(x, quiet=TRUE, refl=refl)))
      } else 
      	ch <- species[which(species$IsChildTaxonOfID == x$TaxonUsageID),]
      if(is.null(ch)) stop('ch is NULL')
       else
        if(nrow(ch)==0) {
          if(!quiet)  if(is.na(x$TaxonName)) message('Could not find ', s) else message(x$TaxonName, ' has no children.')
        } else {
          ch$GENERATION <- 1
          ch2 <- ch
          t <- 1
          repeat {
            t <- t+1
            if(syn) {
              ch2 <- species[which(species$IsChildTaxonOfID == x$TaxonUsageID),'TaxonUsageID']
              ch2 <- do.call(rbind, lapply(ch2, function(x) syn(x, quiet=TRUE)))
            } else 
            	ch2 <- species[which(species$IsChildTaxonOfID %in% ch2$TaxonUsageID),]
            if(nrow(ch2)== 0 ) break
            ch2$GENERATION <- t
            ch <- rbind(ch, ch2)
            if(gen <= t) break
          }
          if(!is.null(gen)) ch <- ch[ch$GENERATION <= gen,]
          if(!quiet) {
    cat('Children of ', s$TaxonName, ' (', s$TaxonUsageID, ')', if(x$TaxonUsageID != s$TaxonUsageID) {paste(" = Synonym of ", x$TaxonName, ' (', x$TaxonConceptID, ')', sep='')},':\n', sep='')
            print(ch[, names(ch)[names(ch) %in% c('TaxonUsageID','TaxonName','TaxonRank','AccordingTo','IsChildTaxonOfID','GENERATION','SYNONYM','EDITSTATUS')]], row.names=FALSE)
          }}
    invisible(ch)
    }
}
childs <- function(...) child(...)


## Parents of a taxon
parent <- function (x, refl = tv.refl(), rank, quiet = FALSE, ...) {
  taxlevels <- factor(c('FOR','VAR','ZUS','SSP','SPE','SGE','SSE','SER','SEC','AGG','GAT','FAM','ORD','UKL','KLA','UAB','ABT','AG2','ROOT'), levels= c('FOR','VAR','ZUS','SSP','SPE','SGE','SSE','SER','SEC','AGG','GAT','FAM','ORD','UKL','KLA','UAB','ABT','AG2','ROOT'), ordered=TRUE)
  species <- tax("all", detailed = TRUE, refl = refl, syn = TRUE, quiet =TRUE, ...)
  if(length(x)>1) {
  	warning('More than one match, using only first.')
  	x <- x[1]
  }
  s <- tax(x, refl = refl, strict = TRUE, quiet = TRUE, ...)
  y <- species[match(s$TaxonConceptID, species$TaxonUsageID),]
  if(y$TaxonUsageID != s$TaxonUsageID) warning('Synonym, will use valid taxon "', y$TaxonName, '" instead.')
  y$GENERATION <- 0
  p <- species[match(unique(y$IsChildTaxonOfID),species$TaxonUsageID),]
  p$GENERATION <- 1
  
  lo <- function(y, p) {
    if(nrow(p)==0) cat(y$TaxonName, 'has no parents.\n') 
    else {
      p2 <- p
      t <- 1
      repeat {
        t <- t+1
        p2 <- species[match(p2$IsChildTaxonOfID,species$TaxonUsageID),]
				if(is.na(p2$TaxonName)) break
        p2$GENERATION <- t
        p <- rbind(p, p2)
        if(p2$TaxonUsageID == 0 ) break
      }}
    return(p)
  }
  
  if(!missing(rank)) {
    if(!rank %in% taxlevels) stop(c('Rank must be one of', rank))
    if(taxlevels[match(rank, taxlevels)] <= taxlevels[match(y$TaxonRank, taxlevels)]) {
      warning('Species is of equal or higher rank than the specified parent level.')
      p <- c(TaxonName='')
    } else {
      P <- lo(y, p)
      # oblig.taxlevels <- factor(c('SPE','GAT','FAM','ORD','KLA','ABT','ROOT'), levels= c('SPE','GAT','FAM','ORD','KLA','ABT','ROOT'), ordered=TRUE)
      #  p$TAXLEVEL <- as.integer(oblig.taxlevels[match(p$TaxonRank, oblig.taxlevels)])
      P <- P[which(P$TaxonRank == rank), ]
#      if(nrow(p) == 0) p <- c(TaxonName='Incertae_sedis')
      #    tv <- oblig.taxlevels[(which(oblig.taxlevels == y$TaxonRank)+1):length(oblig.taxlevels)]
      #    if(!all(tv %in% p$TaxonRank)) 
      cat('Parent level', rank, ' of', y$TaxonName, '(', y$TaxonUsageID, '):\n')
  if(nrow(P)==0) cat('"Incertae sedis" = uncertain placement within this level.\n') 
#      else
#         print(p[,c('TaxonUsageID','TaxonName','AccordingTo','TaxonRank','GENERATION')], row.names=FALSE)
    }
  }  else P <- lo(y, p)
  if(!quiet) {
    cat('Parents of ', s$TaxonName, ' (', y$TaxonUsageID, ')', if(y$TaxonUsageID != s$TaxonUsageID) {paste(" = Synonym of", y$TaxonName, '(', y$TaxonConceptID, ')')}, ':\n', sep='')
    print(P[, names(P)[names(P) %in% c('TaxonUsageID','TaxonName','TaxonRank','IsChildTaxonOfID','GENERATION')]], row.names=FALSE)
  }
  invisible(P)
}

parents <- function(...) parent(...)


# Synonymy swarm of a taxon
syn <- function (x, refl = tv.refl(), quiet=FALSE, ...) {
  species <- tax('all', detailed = TRUE, refl = refl, syn = TRUE, strict = TRUE, quiet = TRUE, ...)
if(is.character(x)) 
	x <- tax(x, refl=refl, strict=TRUE, quiet = TRUE, ...)$TaxonUsageID
  if(length(x) > 1) {
    warning('More than one match, using only first.')
    x <- x[1]
  }
  v <- species[match(x, species$TaxonUsageID),'TaxonConceptID']
  if(length(v)==0 | is.na(v)) stop('No matching species.')
  s <- species[which(species$TaxonConceptID == v),]
  if(!quiet) {
    cat('Name swarm of', s$TaxonName[s$TaxonUsageID == x],':\n')
		if('EDITSTATUS' %in% names(s)) print(s[, c('TaxonUsageID','TaxonName','SYNONYM','EDITSTATUS')]) else
			print(s[, c('TaxonUsageID','TaxonName','SYNONYM')], row.names=FALSE)
    #    print(p[,c(1,3,8,9,12,21)], row.names=FALSE)
  }
  invisible(s)
}


agg <- function(x, refl = 'GermanSL 1.2', species, ...) {
  message('Deprecated function. Using child(x, gen=1) instead.')
  child(x, refl=refl, species=species, gen=1, ... )
}

