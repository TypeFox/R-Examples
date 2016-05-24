`plot.hmmer` <-
function(x, cutoff=NULL, cut.seed=NULL, cluster=TRUE, mar=c(2, 5, 1, 1), cex=1.1, ...) {
  allowed <- c("phmmer", "hmmsearch", "jackhmmer")
  
  if(!any(inherits(x, allowed)))
    stop(paste("please provide the results of a hmmer search of type:",
               paste(allowed, collapse=", ")))

  if(is.null(x$evalue))
    stop("missing evalues")
  
  ##x$mlog.evalue=-log(x$evalue)
  x$mlog.evalue=x$score
  
  panelplot <- function(z=x$mlog.evalue, ylab="-log(Evalue)", gp=gp, ...) {
    z=as.numeric(z)
    plot(z, xlab="", ylab=ylab, col=gps, ...)
    abline(v=gp, col="gray70", lty=3)

    pos=c(rep(3, length(gp))[-length(gp)],2)
    text(  gp, z[gp], 
          labels=paste0("Nhit=",gp ,", x=", round(z[gp])), 
          col="black", pos=pos, cex=cex, ...) ##"gray50"
  }

  nrow <- 1
  if(!is.null(x$kg) & !is.null(x$species))
    nrow=nrow+2
    
  ##- Setup plot arangment
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfcol=c(nrow,1), mar=mar, cex.lab=cex)


  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)


  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  } else {

    if(cluster) {
      ## Ask USER whether to continue with clustering with many hits  
      nhit <- length(x$mlog.evalue)
      if(nhit > 1500) {
        cluster <- readline( paste0(" Note: ", nhit, 
          " hits, continue with TIME-CONSUMING clustering [y/n/q](n): ") )

        cluster <- switch(cluster, y=TRUE, yes=TRUE, q="QUIT", FALSE)
        if(cluster=="QUIT") { stop("user stop") }
      }
    }

    if(is.null(cut.seed)) {
      ## Use mid-point of largest diff pair as seed for
      ##  cluster grps (typical PDB values are ~110)
      cut.seed = mean( x$mlog.evalue[dx.cut:(dx.cut+1)] )
    }

    if(cluster){
      ##- Partition into groups via clustering 
      ##  In future could use changepoint::cpt.var
      hc <- hclust( dist(x$mlog.evalue) )
      if(!is.null(cutoff)) { cut.seed=cutoff } 
      gps <- cutree(hc, h=cut.seed)
    } 

    if(!cluster || (length(unique(gps))==1)) {
      ##- Either we don't want to run hclust or hclust/cutree 
      ##   has returned only one grp so here we will divide   
      ##   into two grps at point of largest diff
      gps = rep(2, length(x$mlog.evalue))
      gps[1:dx.cut]=1
    }
  }

  gp.inds <- na.omit(rle2(gps)$inds)
  gp.nums <- x$mlog.evalue[gp.inds]

  cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
      "           Yielding Nhits:   ", gp.inds, "\n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to cut.seed
    i <- which.min(abs(gp.nums - cut.seed))
    cutoff <- floor( gp.nums[ i ] )
  }

  inds <- x$mlog.evalue >= cutoff
  cat("  * Chosen cutoff value of:   ", cutoff, "\n",
      "           Yielding Nhits:   ", sum(inds), "\n")

  ##- Plot each alignment statistic with annotated grps
  ##panelplot(gp=gp.inds)
  panelplot(x$score, ylab="Bitscore", gp=gp.inds)
  

  ## plot kigdom / species
  if(!is.null(x$kg) & !is.null(x$species)) {
    tmpfun <- function(s) {
      s=s[1:min(2, length(s))]
      if(length(s)>1)
        paste(substr(s[1], 1,1), substr(s[2], 1,6), collapse=".")
      else
        substr(s, 1,8)
    }
    
    ## make grps for x$species
    species <- unlist(lapply(strsplit(x$species, " "), tmpfun))
    grps.sp <- rep(NA, length(species))
    unq.sp <- unique(species)
    
    for(i in 1:length(unq.sp))
      grps.sp[ which(species %in% unq.sp[i]) ] = i

    ## make grps for x$kg
    grps.kg <- rep(NA, length(x$kg))
    unq.kg <- unique(x$kg)
    for(i in 1:length(unq.kg))
      grps.kg[ which(x$kg %in% unq.kg[i]) ] = i

    ylim <- c(0, max(x$score)*1.2)
    xlim <- c(0, length(x$score))
    
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, ...)
  
    cols <- c("lightblue", "lightgreen", "lightpink", "lightsalmon",
              "lightyellow3",  "lightcoral", "lightsteelblue2", "lightgoldenrod1")
          
    for(i in 1:length(unq.kg)) {
      bs <- bounds(which(grps.kg==i))
      rect(bs[,"start"]-1, 1, bs[,"end"], max(x$score)*1.05, col=cols[i], border=NA)
    }
    
    mp <- barplot(x$score, col=grps.sp, width=1.0, space=0, border=par("fg"),
                  ylab="Annotation",  add=TRUE, ...)
    abline(v=gp.inds, col="gray70", lty=3)

    ncol <- 2
    legend("topright", unq.kg, col=cols[1:i],  pch=16,
           ncol=ncol, cex=cex*0.8, box.lwd = 0, box.col = "white",bg = "white")
    box(); axis(1);
    
    ## summary of chosen stuff
    unq.kg=unique(x$kg[which(inds)])
    unq.sp=unique(x$species[which(inds)])
    
    txt1 <- paste("N kingdoms: ", length(unq.kg), "/", length(unique(x$kg)), sep="")
    txt2 <- paste("N species: ", length(unq.sp), "/", length(unique(x$species)), sep="")
    txt <- paste(txt1, txt2, sep=", ")

    mtext(txt, side=3, line=-1.25, at=0, adj=0, cex=cex*0.8)
  }


  ## plot kigdom / species
  if(!is.null(x$kg) & !is.null(x$species)) {
    unq.kg <- unique(x$kg)
    tbl <- table(x$kg, cut(x$score, 30))
    tbl=tbl[, seq(ncol(tbl), 1), drop=FALSE]
    cols <- seq(1,nrow(tbl))
    barplot(tbl, col=cols, ylab="Frequency")
    
    legend("topright", rownames(tbl), col=cols,  pch=16,
           cex=cex*0.8, box.lwd = 0, box.col = "white",bg = "white")
    box()
  }

 
  
  ##- Return details of hits above cutoff
  out <- cbind("acc"=x$acc[inds], "group"=gps[inds])
  rownames(out) <- which(inds)
  invisible(list(hits=out, acc=x$acc[inds], "inds"=which(inds)))
}
