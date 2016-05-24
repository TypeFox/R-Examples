plot.simfam <- function(x, famid=NULL, pdf=FALSE, file=NULL, ...){
  
  if(is.null(famid)) famid <- x$famID[1]
  if(length(famid)<5) {
  	ncol <- length(famid)
  	nrow <- 1
  	}
  else {
  	ncol <- 5 
  	nrow <- ceiling(length(famid)/5)
	}

  if(is.character(file)) pdf=TRUE

  if(pdf){
	  par(mfrow=c(1,1))
	  if(is.null(file)) pdf("pedigreeplot.pdf")
	  else pdf(file=file)
  }
  else par(mfrow=c(nrow, ncol))


  for(i in famid){
  obj <- x[x$famID==i, ]
  obj$gender[obj$gender==0] <- 2
  iprob <- rep("",length(obj$indID))
  iprob[obj$proband==1] <- "proband"
  if(!is.null(obj$carrp.pheno)) iprob[is.na(obj$mgene)] <- paste("p=",round(obj$carrp.pheno[is.na(obj$mgene)], 2),sep="")
  else if(!is.null(obj$carrp)) iprob[is.na(obj$mgene)] <- paste("p=",round(obj$carrp[is.na(obj$mgene)], 2),sep="")

  ped <- pedigree(id=obj$indID, dadid=obj$fatherID, momid=obj$motherID, sex=obj$gender, affected=cbind(obj$status, obj$mgene))

 pedplot<- plot.pedigree(ped, id=paste(obj$indID, iprob, sep="\n"), col=ifelse(obj$proband==1,2,1))

 title(paste("Family ID:", i))
 pedigree.legend(ped, labels=c("Affected", "Mutation\nCarrier"), location="topright", radius=pedplot$boxw*0.7)
}


if(pdf) dev.off()

}