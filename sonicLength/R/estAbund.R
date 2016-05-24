estAbund <- function(locations, lengths, replicates=NULL, jackknife=F, kmax=0, ...)
{
  ## Purpose:estimate abundances from sonicated lengths
  ## ----------------------------------------------------------------------
  ## Arguments: locations - vector of location IDS
  ##            lengths - vector of lengths
  ##            replicates - vector of replicate IDs
  ##            jackknife - return leave one out results?
  ##            kmax - highest count to bother with (all higher values
  ##                        are globbed together in the result)
  ##            ... - other args to pass to maxEM or phi.update
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date: 27 May 2011, 13:04

  ## TODOs: minlen arg?
  
  mc <- match.call()
  noRep <- (is.null(replicates))
  stopifnot( noRep || length(replicates)==length(lengths) )
  stopifnot( noRep || length(unique(replicates))!=1 )
  stopifnot( if ( noRep )  !jackknife else TRUE )
  stopifnot( length(locations) == length(lengths) )

  dfr.call <- mc
  dfr.call[[1]] <- as.name("data.frame")
  keep.args <- c(1,pmatch(c("locations","lengths", if (!missing(replicates)) "replicates"),names(dfr.call)))
  dfr.call[ -keep.args ] <- NULL
  dfr <- eval.parent(dfr.call)
  if (noRep) {
  ## simplest case is no replicates

    emcall <- mc
    emcall$locations <- emcall$lengths <- emcall$replicates <- NULL
    emcall[[1]] <- as.name("maxEM")
    emcall$slmat <- quote(
                          xtabs( ~factor(lengths,min(lengths):max(lengths))
                                + factor(locations) ) )
    tmp <- eval(emcall)
    tmp$call <- mc
    tmp$obs <- xtabs( ~ factor(locations) )
    tmp$pr.y <-
      if (kmax==0) NULL else pr.y.given.x(eval(emcall$slmat), tmp$theta , tmp$phi, kmax)
    tmp$data <- dfr
    tmp
  } else {

    ## with replicates
    sl.start <- unique(locations) # all insertion sites
    slmat.list <- list()
    for (ix in sort(unique(replicates))) {
      slmat.list[[ ix ]] <-
        xtabs( ~factor(lengths,min(lengths):max(lengths)) +
              factor(locations,sl.start),
              subset=replicates == ix )
    }
    ## combine tables
    slmat3 <- do.call(rbind,slmat.list)
    
    ## get the lframe for phi.update
    
    tframe.list <- lapply(slmat.list,function(x) pad.tab(rowSums(x)))
    tframe3 <- do.call(rbind,tframe.list)
    ## add a strata column
    tframe3$strata <- factor( rep(1:length(tframe.list), sapply( tframe.list, nrow ) ) )
    ## slim down the rownames:
    rownames(tframe3) <- with(tframe3, paste(strata,x)) # WAS: sub(".* ","",rownames(tframe3))
    ## copy them to slmat3
    rownames(slmat3) <- rownames(tframe3)[tframe3$orig]
    if (length(levels(tframe3$strata))>20) warning( "more than 20 replicates specified...")
    ## now fit the whole thing:

    emcall <- mc
    emcall$locations <- emcall$lengths <- emcall$replicates <- NULL
    emcall[[1]] <- as.name("maxEM")
    emcall$slmat <- as.name("slmat3")
    emcall$lframe <- as.name("tframe3")
    ## be sure partial matching doesn't bite:
    if (!any(pmatch( names(emcall), c("phi.update","phi.deriv"),0)==1))  emcall$phi.update <- as.name("phi.update.lframe")
    tmp3 <- eval(emcall)
    
    ## typically this works out to: tmp3 <- maxEM(slmat3,lframe=tframe3,phi.update=phi.update.lframe)
    
##    tmp.all <- lapply(slmat.list,function(x) maxEM(x[, colSums(x)>0]))
    
    tmp3$lframe <- tframe3
    tmp3$call <- mc
    tmp3$obs <- xtabs( ~factor(locations,sl.start) + factor(replicates,sort(unique(replicates))) )
    tmp3$pr.y <-
      if (kmax==0) NULL else pr.y.given.x(slmat3, tmp3$theta , tmp3$phi, kmax)
    if (jackknife){
      ureps <- sort(unique(replicates))
      jack.res <-
        lapply( ureps, function(x) {
          lf.use <- tframe3$strata != x
          slmat.rowuse <- lf.use[ tframe3$orig ]
          slmat.coluse <- colSums(slmat3[ slmat.rowuse, ])>0
          emcall$slmat <- quote(slmat3[slmat.rowuse,slmat.coluse])
          emcall$lframe <- quote(tframe3[lf.use,])
          jres <- eval(emcall)
          jres$pr.y <-
            if (kmax==0) NULL else pr.y.given.x(eval(emcall$slmat), jres$theta , jres$phi, kmax)
          jres
        })
      tmp3$jackknife <- jack.res
    }
    tmp3$data <- dfr
    tmp3
  }    
}


