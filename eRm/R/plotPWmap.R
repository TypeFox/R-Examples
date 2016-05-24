`plotPWmap` <-
function(object, pmap=FALSE, imap=TRUE, item.subset="all", person.subset="all",
                 mainitem="Item Map", mainperson="Person Map",
                 mainboth="Item/Person Map", latdim="Latent Dimension",
                 tlab="Infit t statistic", pp=NULL, cex.gen=0.6, cex.pch=1,
                 person.pch=1, item.pch=16, personCI=NULL, itemCI=NULL, horiz=FALSE)
{
#mjm  def.par <- par(no.readonly = TRUE) ## save default, for resetting...

  ## Pathway map currently only for RM, PCM and RSM

  ## The next part of the code finds locations and standard errors for
  ## the item thresholds
  if ((object$model == "LLTM") || (object$model == "LRSM") || (object$model == "LPCM"))
    stop("Pathway Map can only be computed for RM, RSM, and PCM!")

  if (!pmap && !imap)
    stop("Pathway Map requires you to request at least one map (item or person)!")

  ## compute threshtable (from betapars for dichotomous models) and item names
  if (object$model == "RM" || max(object$X, na.rm=TRUE) < 2 ) { # dichotomous model
    dRm <- TRUE

    ## betapars are easiness parameters; only the pars need negating
    threshtable<-cbind(object$betapar * -1, object$se.beta)
    rownames(threshtable) <- colnames(object$X)

    ## shorter synonym
    tt<-threshtable
  } else { ## polytomous model
    dRm <- FALSE

    thresh <- thresholds(object)
    threshtable <- cbind(thresh$threshpar, thresh$se.thresh)
    tlevels<-apply(thresh$threshtable[[1]], 1,
                   function(x) length(na.exclude(x))) - 1
    if (!(sum(tlevels)==nrow(threshtable)))
      stop("Threshtable rows not equal to number of thresholds - oops!")

    ttl<-NULL ## threshtable labels
    for (i in rownames(as.matrix(tlevels)))
      if (tlevels[i]==1)
        ttl<-c(ttl,i)
      else
        ttl<-c(ttl,paste(i,1:tlevels[i],sep=":"))
    rownames(threshtable)<-ttl

    ## shorter synonyms
    tt<-threshtable
    tl<-tlevels
  }

  if(is.null(pp)){
    suppressWarnings(pp<-person.parameter(object))
  } else if(class(pp) != "ppar" || !identical(pp$X,object$X)){
    stop("pp is not a person.parameter object which matches the main Rasch data object!")
  }
  


  ## We will be plotting the infit data versus the parameters for
  ## both items and persons
  # item fit
  iloc  <- tt[, 1L]
  ise   <- tt[, 2L]
  ifit  <- itemfit(pp)
  ifitZ <- ifit$i.infitZ

  # person fit
  pfit       <- personfit(pp)
  pfitZ      <- pfit$p.infitZ
  if(length(pfit$excl_obs_num) > 0L){                  # mjm 2014-09-17
    temp_namevec <- pfit$excl_obs_chr                  #
    for(ex_pers in pfit$excl_obs_num){                 # workaround: add deleted persons and their names
      pfitZ <- append(pfitZ, NA, ex_pers - 1L)         #
      names(pfitZ)[ex_pers] <- temp_namevec[ex_pers]   #
    }                                                  #
    rm(temp_namevec)                                   #
  }                                                    #
  
  ploc       <- as.matrix(pp$theta.table['Person Parameter'])[,1]

  if(length(pfit$excl_obs_num) > 0L){   # mjm 2014-09-17
    ploc[pfit$excl_obs_num] <- NA       # set all parameters of excluded obs NA
  }                                     #
  
  pse        <- unlist(pp$se.theta, recursive=FALSE)
  names(pse) <- sub("^NAgroup[0-9]*\\.","",names(pse))
  pse        <- pse[names(ploc)]

  ## We can now do item and person subsetting; the item subsetting is
  ## pretty ugly as there are multiple cases.  (We dare not do it earlier
  ## as we have to take items from all of iloc, ise and ifitZ.)
  if (imap && is.character(item.subset)) {
    ## Case 1: item subsetting by item names
    if (dRm) {
      if (length(item.subset)>1 && all(item.subset %in% rownames(tt))) {
        iloc  <- iloc[item.subset]
        ise   <- ise[item.subset]
        ifitZ <- ifitZ[item.subset]
        tt    <- tt[item.subset,]
      }
      else if(length(item.subset)!=1 || !(item.subset=="all"))
        stop("item.subset misspecified. Use 'all' or vector of at least two valid item indices/names.")
    } else {
      if (length(item.subset)>1 && all(item.subset %in% rownames(as.matrix(tl)))) {
        tl    <- tl[item.subset]
        keep.subset <- c()
        for (i in rownames(as.matrix(tl)))
          if (tl[i]==1)
            keep.subset<-c(keep.subset,i)
          else
            keep.subset<-c(keep.subset,paste(i,1:tl[i],sep=":"))

        iloc  <- iloc[keep.subset]
        ise   <- ise[keep.subset]
        ifitZ <- ifitZ[item.subset]
        tt<-tt[keep.subset,]
      }
      else if(length(item.subset)!=1 || !(item.subset=="all"))
        stop("item.subset misspecified. Use 'all' or vector of at least two valid item indices/names.")
    }
  } else if (imap) {
    ## Case 2: item subsetting by item numbers
    if (dRm) {
      if (length(item.subset)>1 && all(item.subset %in% 1:nrow(tt))) {
        iloc  <- iloc[item.subset]
        ise   <- ise[item.subset]
        ifitZ <- ifitZ[item.subset]
        tt    <- tt[item.subset,]
      }
      else
        stop("item.subset misspecified. Use 'all' or vector of at least two valid item indices/names.")
    }
    else {
      if (length(item.subset)>1 && all(item.subset %in% 1:length(tl))) {
        tl    <- tl[item.subset]
        keep.subset <- c()
        for (i in rownames(as.matrix(tl)))
          if (tl[i]==1)
            keep.subset<-c(keep.subset,i)
          else
            keep.subset<-c(keep.subset,paste(i,1:tl[i],sep=":"))

        iloc  <- iloc[keep.subset]
        ise   <- ise[keep.subset]
        ifitZ <- ifitZ[item.subset]
        tt<-tt[keep.subset,]
      }
      else
        stop("item.subset misspecified. Use 'all' or vector of at least two valid item indices/names.")
    }
  }

  ## We can now do person subsetting; this is significantly easier than
  ## item subsetting, as there is no dRM/eRm distinction.
  if (pmap && is.character(person.subset)) {
    ## Case 1: person subsetting by person names
    if (length(person.subset)>1 && all(person.subset %in% names(ploc))) {
      ploc  <- ploc[person.subset]
      pse   <- pse[person.subset]
      pfitZ <- pfitZ[person.subset]
    }
    else if(length(person.subset)!=1 || !(person.subset=="all"))
      stop("person.subset misspecified. Use 'all' or vector of at least two valid person indices/names.")
  } else if (pmap) {
    ## Case 2: person subsetting by person numbers
    if (length(person.subset)>1 && all(person.subset %in% seq_along(ploc))) {
      ploc  <- ploc[person.subset]
      pse   <- pse[person.subset]
      pfitZ <- pfitZ[person.subset]
    }
    else
      stop("person.subset misspecified. Use 'all' or vector of at least two valid person indices/names.")
  }

  ## Confidence intervals for persons and items
  ##
  ## Need defaults for multiple of standard error for purpose of range
  ## calculation; these are zero as default is not to draw confidence
  ## intervals
  pci <- 0
  ici <- 0

  ## Our calculation is simplistic; we use the normal distribution to
  ## estimate our confidence interval from our standard error.  However,
  ## since this is likely to only be approximate and indicative anyway, we
  ## are not concerned by this.
  if(pmap && !is.null(personCI)) {
    if(is.null(personCI$clevel)) personCI$clevel <- 0.95
    if(is.null(personCI$col))    personCI$col    <- "orange"
    if(is.null(personCI$lty))    personCI$lty    <- "dotted"
    pci <- qnorm((1-personCI$clevel)/2, lower.tail=FALSE)
  }
  if(imap && !is.null(itemCI)) {
    if(is.null(itemCI$clevel)) itemCI$clevel <- 0.95
    if(is.null(itemCI$col))    itemCI$col    <- "red"
    if(is.null(itemCI$lty))    itemCI$lty    <- "dotted"
    ici <- qnorm((1-itemCI$clevel)/2, lower.tail=FALSE)
  }

  ## Now we can plot the Pathway Map

  if (pmap) { ## person map
    xrange.pmap <- range(pfitZ,finite=TRUE)
    xrange.pmap[1] <- min(-2.5,xrange.pmap[1])
    xrange.pmap[2] <- max(2.5,xrange.pmap[2]+1) ## need space for labels
    yrange.pmap<-range(ploc,finite=TRUE)
    yrange.pmap[1]<-yrange.pmap[1]-pci*max(pse, na.rm=TRUE)
    yrange.pmap[2]<-yrange.pmap[2]+pci*max(pse, na.rm=TRUE)
  }
  if (imap) { ## item map
    xrange.imap <- range(ifitZ,finite=TRUE)
    xrange.imap[1] <- min(-2.5,xrange.imap[1])
    xrange.imap[2] <- max(2.5,xrange.imap[2]+1) ## need space for labels
    yrange.imap<-range(iloc,finite=TRUE)
    yrange.imap[1]<-yrange.imap[1]-ici*max(ise, na.rm=TRUE)
    yrange.imap[2]<-yrange.imap[2]+ici*max(ise, na.rm=TRUE)
  }

  if (pmap && !imap) {
    xrange <- xrange.pmap
    yrange <- yrange.pmap
    maintitle <- mainperson
  } else if (!pmap && imap) {
    xrange <- xrange.imap
    yrange <- yrange.imap
    maintitle <- mainitem
  } else if(pmap && imap){
    xrange <- numeric(2)
    yrange <- numeric(2)
    xrange[1] <- min(xrange.pmap[1], xrange.imap[1], na.rm = TRUE)
    xrange[2] <- max(xrange.pmap[2], xrange.imap[2], na.rm = TRUE)
    yrange[1] <- min(yrange.pmap[1], yrange.imap[1], na.rm = TRUE)
    yrange[2] <- max(yrange.pmap[2], yrange.imap[2], na.rm = TRUE)
    maintitle <- mainboth
  } else {
    stop("error determining plot ranges.")
  }


#mjm  par(mar=c(5,4,4,2))

  if (!horiz){  # rh 2010-12-09
    plot(xrange, yrange, xlim = xrange, ylim = yrange, main = maintitle, ylab = latdim, xlab = tlab, type = "n")
    abline(v = c(-2, 2), col = "lightgreen")
  } else {
    plot(yrange, xrange, xlim = yrange, ylim = xrange, main = maintitle, ylab = tlab, xlab = latdim, type = "n")
    abline(h = c(-2, 2), col = "lightgreen")
  }

  if (pmap) { ## person map
    zt <- pfitZ
    if (!horiz){
      if (pci>0) ## draw confidence intervals
        arrows(zt,ploc+pci*pse, zt,ploc-pci*pse, angle=90, code=3, length=0.04,
               col=personCI$col, lty=personCI$lty)
      points(zt,ploc,pch=person.pch,cex=cex.pch)
      text(zt,ploc,names(ploc),cex=cex.gen,pos=4)
    } else {
      if (pci>0) ## draw confidence intervals
        arrows(ploc+pci*pse, zt,ploc-pci*pse, zt, angle=90, code=3, length=0.04,
               col=personCI$col, lty=personCI$lty)
      points(ploc, zt, pch=person.pch,cex=cex.pch)
      text(ploc, zt, names(ploc),cex=cex.gen,pos=4)
    }
  }



  if (imap) { ## item map
    if (dRm)
      zt <- ifitZ
    else
      zt <- rep(ifitZ,times=tl)

    if (!horiz){
      if (ici>0) ## draw confidence intervals
        arrows(zt,iloc+ici*ise, zt,iloc-ici*ise, angle=90, code=3, length=0.04,
               col=itemCI$col, lty=itemCI$lty)
      points(zt,iloc,pch=item.pch,cex=cex.pch)
      text(zt,iloc,rownames(tt),cex=cex.gen,pos=4)
    } else {
      if (ici>0) ## draw confidence intervals
        arrows(iloc+ici*ise, zt,iloc-ici*ise,zt, angle=90, code=3, length=0.04,
               col=itemCI$col, lty=itemCI$lty)
      points(iloc, zt,pch=item.pch,cex=cex.pch)
      text(iloc,zt, rownames(tt),cex=cex.gen,pos=4)
    }
  }

#mjm  par(def.par)

  invisible(NULL)
}
