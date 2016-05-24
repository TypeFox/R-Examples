plot.minPtest <-
function(x, type=c("gene", "SNP", "both") ,level=0.05, lambda=1, gene.name=FALSE, sigPch=pch, nonsigPch=pch, pch=20, sigLty=lty, nonsigLty=lty, lty=1, sigCol=col, nonsigCol=col, col=NULL, xlab, ...){
  type <- match.arg(type)
  if(is.null(sigCol)) sigCol <- "red"
  if(is.null(nonsigCol)) nonsigCol <- "black"
  if (missing(xlab)) {
    xlab <- ifelse(type == "gene", "Gene", "SNP")
  }
  minp <- x$minp
  genes <- rownames(minp)
  if(length(which(minp[,"minP"]==0))>0){# to check
    minp[which(minp[,"minP"]==0),"minP"] <- 1/x$n.permute
  }
  p.adj.minp <- x$p.adj.minp
  minp.padj <- cbind(minp,p.adj.minp)
  log.minp <- -log(minp[,"minP"],10)
  ## permutation-based p-values only
  if(type=="gene"){
    plot(x=c(1,x$nrgene), y=c(0,max(log.minp)), xlab=xlab, ylab=expression(-log[10](minp)),xaxt="n",type="n", ...)
    for(i in 1:x$nrgene){
      if(minp.padj[i,"p.adjust"]<=level){
        points(i,log.minp[i], pch=sigPch, col=sigCol, ...)
      }else{
        points(i,log.minp[i], pch=nonsigPch, col=nonsigCol, ...)
      }
    }
    axis(1, at=1:x$nrgene,labels=genes,...)
    box(...)
  }else{
    psnp <- x$psnp
    SNPtoGene <- x$SNPtoGene
    p.adj.snp <- x$p.adj.psnp
    psnp.padj <- cbind(psnp,p.adj.snp)
    gensnp <- lapply(seq_along(genes), function(i){
      snpnames <- SNPtoGene[which(SNPtoGene[,2]==genes[i]),1]
      psnp.padj.temp <- psnp.padj[snpnames,]
      if(length(snpnames)==1){
        psnp.padj.temp <- matrix(psnp.padj.temp,nrow=1,ncol=2)
        rownames(psnp.padj.temp) <- snpnames
        colnames(psnp.padj.temp) <- c("p_value","p.adjust")
      }
      psnp.padj.temp
    })
    names(gensnp) <- genes
    gensnp.sort <- do.call(rbind,gensnp)
    log.psnp <- -log(gensnp.sort[,"p_value"],10)
    ## marginal p-values only
    if(type=="SNP"){
      plot(x=c(0,x$nrsnp),y=c(0,max(log.psnp)), xlab=xlab, ylab=expression(-log[10](psnp)),xaxt="n",type="n", ...)
      for(i in 1:x$nrsnp){
        if(gensnp.sort[i,"p.adjust"]<=level){
          points(i,log.psnp[i], pch=sigPch, col=sigCol, ...)
        }else{
          points(i,log.psnp[i], pch=nonsigPch, col=nonsigCol, ...)
        }
      }
      if(gene.name){
        gene.size <- unlist(lapply(gensnp, function(i){
          l <- nrow(i)
        }))
        x1 <- 1
        x2 <- gene.size[1]
        middle <- median(x1:x2)
        for(i in 1:(length(gene.size))){
          x1 <- x2+1
          x2 <- x2+gene.size[i+1]
          if(i!=length(gene.size)){
            middle <- c(middle,median(x1:x2))
          }
        }
        axis(1, at=middle,labels=genes, ...)
      }
      box(...)
    }
    ## permutation-based p-values and marginal p-values only
    if(type=="both"){
    par(mar=c(5,4,4,5))
    plot(x=c(0,x$nrsnp),y=c(0,max(log.psnp)), xlab=xlab, ylab=expression(-log[10](psnp)),xaxt="n",type="n", ...)
    for(i in 1:x$nrsnp){
      if(gensnp.sort[i,"p.adjust"]<=level){
        points(i,log.psnp[i], pch=sigPch, col=sigCol, ...)
      }else{
        points(i,log.psnp[i], pch=nonsigPch, col=nonsigCol, ...)
      }
    }
    gene.size <- unlist(lapply(gensnp, function(i){
      l <- nrow(i)
      }))
    x1 <- 1
    x2 <- gene.size[1]
    middle <- median(x1:x2)
    for(i in 1:(length(gene.size))){
      par(new=TRUE)
      plot(x=c(0,x$nrsnp),y=c(0,max(log.minp)), xlab="", ylab="",xaxt="n",type="n", yaxt="n", ...)
      if(minp.padj[i,"p.adjust"]<=level){
        lines(c(x1,x2),c((log.minp[i])*lambda,(log.minp[i])*lambda),col=sigCol, lty=sigLty, ...)
      }else{
        lines(c(x1,x2),c((log.minp[i])*lambda,(log.minp[i])*lambda), col=nonsigCol, lty=nonsigLty, ...)
      }
      x1 <- x2+1
      x2 <- x2+gene.size[i+1]
      if(i!=length(gene.size)){
        middle <- c(middle,median(x1:x2))
      }
    }
    axis(4, at=seq(0,max(log.minp)), labels=seq(0,max(log.minp)), ...)
    if(lambda==1){
      mtext(text=expression(-log[10](minp)), side=4, line=3, ...)
    }else{
      mtext(text=bquote(.(-lambda) * (log[10](minp))), side=4, line=3, ...)
    }
    if(gene.name){
      axis(1, at=middle,labels=genes, ...)
    }
    box(...)
  }}
}
