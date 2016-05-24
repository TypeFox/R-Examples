`print.eqtl` <- function(x,which=NULL,sig=0.01,output="bed",...){
   
  output <- match.arg(output,c("bed","single"))

  if(x$type=="full")
  {
    x <- x$eqtl
    xx <- x
    if(!is.null(which)==TRUE) xx <- x[which]

    X <- data.frame(Chr=-1,Start=-1,End=-1,Name="-1", Gene="-1")
    
    if(output=="bed")
    {
    for(testRun in 1:length(xx))
    {
      tempRes <- xx[[testRun]]
      repSNP <- tempRes[[2]][((tempRes[[3]]>-1)&(tempRes[[3]]<sig)),]
      if(nrow(repSNP)>0){
        temp <- data.frame(Chr=repSNP[,1],Start=repSNP[,4],End=repSNP[,4],Name=repSNP[,2], Gene=names(xx)[testRun])
        X <- rbind(X,temp)
      }
    }
      X <- X[-1,]
    } else if(output=="single") {
      X <- list()
      
      for(gene in 1:length(xx))
      {
	  temp <- xx[[gene]]
	  Nlocs <- length(table(temp$ProbeLoc))
	  for(sub in 1:Nlocs)
	  {
	    temp <- xx[[gene]]$TestedSNP[((xx[[gene]]$p.values>-1)&(xx[[gene]]$p.values<=sig)),c(1,2,4,5,6)]
	    temp2 <- xx[[gene]]$p.values[((xx[[gene]]$p.values>-1)&(xx[[gene]]$p.values<=sig))]
	    temp <- cbind(temp,temp2)
	    X[[gene]] <- temp
	    colnames(X[[gene]]) <- c("Chr","SNP","Position","Allele1","Allele2","p.value")
	  }
      }
      if(is.null(which)) which <- 1:length(x)
      names(X)<- names(x[which])
    }
    print(X,...)
  } else if(x$type=="sig"){
    print(x$bed,...)
  }
} 

