`grpPhylogeo` <-
function(gbinfo,align,seuil=0.1,method="single",model="raw",pairwise.deletion=TRUE)
{
  al <- read.dna(align,format="fasta")
  gb <- read.csv(gbinfo,header=TRUE,colClasses="character")
    
  nomseq <- names(al)
  gbal <- gb[match(nomseq,gb$nom),]

  gbal$lon <- gsub(",",".",gbal$lon)
  gbal$lat <- gsub(",",".",gbal$lat)

  gbal$lon <- as.numeric(gbal$lon)
  gbal$lat <- as.numeric(gbal$lat)
  
  gbal <- gbal[!is.na(gbal$lat),]
  gbal <- gbal[!is.na(gbal$lon),]
  
  al <- al[gbal$nom]
  
  gdist <- dist.dna(al,as.matrix=TRUE,model=model,pairwise.deletion=pairwise.deletion)
  pdist <- distGPS(gbal)
  
  hc <- hclust(as.dist(pdist),method=method)
  
  plot(hc,main="",xlab="",ylab="physical distance")
  
  mhc <- cutree(hc,h=hc$height)
  
  edgeCor <- list(NULL)
  
  for (i in 1:length(hc$height))
  {
    groupe <- mhc[,i]
    nbgroupe <- unique(groupe)
    
    corM <- NULL
    
    for (j in 1:length(nbgroupe))
    {
      sousGrp <- names(groupe[groupe==nbgroupe[j]])
      
      m1 <- as.matrix(gdist[sousGrp,sousGrp])
      m2 <- as.matrix(pdist[sousGrp,sousGrp])
      
      Z <- mantel.test(m1,m2)
      
      corM <- c(corM,Z$p)
    }
    
    edgeCor[[i]] <- corM
  }
  
  edgeRec <- list(NULL)
  
  for(i in 1:length(edgeCor))
  {
    a <- (1:length(edgeCor[[i]]))[edgeCor[[i]]<seuil]
    
    edgeRec[[i]] <- a
  }

  out <- list(NULL)
  for (i in 1:length(edgeRec))
  {
    if(length(edgeRec[[i]])>0)
    {
      nb <- unique(edgeRec[[i]])
      grpframe <- cutree(hc,h=hc$height[i])
      for (j in 1:length(nb))
      {
      grpframe2 <- grpframe[grpframe==nb[j]]
      nom <- names(grpframe2)
      nom2 <- match(nom,hc$labels)
      deb <- min(match(nom2,hc$order)) -0.5
      fin <- max(match(nom2,hc$order)) + 0.5
      #col=heat.colors(seuil*100)[seuil*100-floor(edgeCor[[i]][nb[j]]*100)]
      col="Red"
      lines(c(deb,deb),c(0,hc$height[i]),col=col)
      lines(c(fin,fin),c(0,hc$height[i]),col=col)
      lines(c(deb,fin),rep(hc$height[i],2),col=col)
      out[[length(out)+1]] <- nom
      }
    }
  }
out2 <- unique(out[2:length(out)])
out2
}

