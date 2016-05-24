if(getRversion() >= "2.15.1")  utils::globalVariables(c("multipatt"))

syntab <- function (veg, clust, type = c('rel','abs','mean.cover'), mupa = NULL, dec = 0, refl, ...)
{
    type <- match.arg(type)
    if (missing(clust)) clust <- sample(1:2, size = nrow(veg), replace = TRUE)
    ncl <- length(unique(clust))
    cat(' Number of clusters: ', ncl, '\n')
    nb.rel.clust <- as.numeric(table(clust))
    cat(' Cluster frequency', nb.rel.clust,'\n')
    if(any(nb.rel.clust == 0)) stop("All cluster levels must be represented by plots.")
    if(is.null(levels(clust))) levels(clust) <- 1:length(table(clust))
    if(any(colSums(veg)==0)) stop('Some species without occurrence.')
    sp.veg <- split(veg, clust, drop=FALSE)
    if(type=='rel') {
    	tmp <- lapply(sp.veg, function(x) colSums(x > 0))
    	temp <- vector('list', length=ncl)
	    for(i in 1:length(nb.rel.clust))
	      temp[[i]] <- round(tmp[[i]] / nb.rel.clust[i] * 100, dec)
    }
    if(type=='mean.cover') {
      temp <- lapply(sp.veg, function(x) {x[x==0] <- NA; round(colMeans(x, na.rm=TRUE),dec)})
      is.na(temp) <- 0
    }
    if(type=='abs') 
      temp <- lapply(sp.veg, function(x) colSums(x > 0))

    st <- as.data.frame(temp)
    st[is.na(st)] <- 0
    names(st) <- levels(clust)

    if(!is.null(mupa) | class(mupa)=='multipatt' & ncl < 1) {
      if(class(mupa)!='multipatt') {
        if (requireNamespace("indicspecies", quietly = TRUE)) {
          mu <- indicspecies::multipatt(veg, clust, ...)
        } else {
          stop("Package indicspecies needed (function multipatt)")
        }} else mu <- mupa
      # st <- st[rownames(st) %in% rownames(sig),]     
      # o <- order(mu$sign[,'index'])
      df <- mu$sign
      df[, 1:ncl] <- st # mu$sign[,1:ncl] * st
      colnames(df)[1:ncl] <- levels(clust)
      st <- df
    }
      
    out <- list(syntab=st, clust=clust)
    class(out) <- c('syntab', 'list')
    return(out)
}

#--------------

if(getRversion() >= "2.15.1")  utils::globalVariables(c("st"))

#--------------
print.syntab <- function(x, zero.print='.', trait, limit = 1, minstat = 0, alpha = 0.05, ...) {
  clust <- x[[2]]
  x <- x$syntab
  if(any(c('stat','index','p.value') %in% names(x))) {
    if(any(is.na(x[1:(ncol(x)-3)]))) stop('NA values in frequency table. Do you have species without occurrences in your matrix?')
    mu = TRUE 
    } else {
      if(any(is.na(x))) stop('NA values in frequency table. Do you have species without occurrences in your matrix?')
      mu = FALSE
    }
  if(mu) {
    stat <- x[,'stat']; x <- x[,-which(names(x)=='stat')]
    index <- x[,'index']; x <- x[,-which(names(x)=='index')]
    p.value <- x[,'p.value']; x <- x[,-which(names(x)=='p.value')]
    select <- stat > minstat & !is.na(stat) & p.value < alpha & !is.na(p.value) & apply(x, 1, function(y) max(y) >= limit)
    } else  select <- apply(x, 1, function(y) max(y) >= limit)
  if(sum(select)>0) {
    x <- x[select, ] 
  if(zero.print != "0" && any(i0 <- x == 0)) {
      x[i0] <- sub("0", zero.print, x[i0])
      x[i0] <- sub("0.0", zero.print, x[i0])
      x[i0] <- sub("0.00", zero.print, x[i0]) 
  }

  if(mu) {
    if(sum(select) > 0)
      x <- cbind(x, index=index[select], stat=stat[select], p.value=p.value[select])
    x <- x[order(x$index),]
    }

  if(!missing(trait)) {
    if(is.null(names(trait))) stop('Trait vector must have names of taxa according to the vegetation matrix.')
    traitname <- names(trait) #as.character(substitute(trait))
    trait.df <- as.data.frame(trait[match(rownames(x), rownames(trait)),])
		x <- cbind(x, trait.df)
#		names(x)[names(x)=='trait'] <- traitname
	}
  } else warning('NO species exceed the chosen significance threshold.')

  ncl <- length(unique(clust))
  cat(' Number of clusters: ', ncl, '\n')
  nb.rel.clust <- as.numeric(table(clust))
  cat(' Cluster frequency', nb.rel.clust,'\n')

  if(sum(select)>0) print.data.frame(x, ...)
 }

