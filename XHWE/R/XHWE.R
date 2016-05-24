XHWE <-
function(ped,loci=NULL,header=T,status_missing=9,allele_missing=0,start.rho=0.02,simuno=1000,dv=1e-7,itertime=1000,filename="results.txt") {
  
  marker.name <- NULL
  
  # read in pedigree file if it is provided
  
  if (is.character(ped)) {
    
    pedfile <- ped

    ped <- read.table(pedfile,header=header,stringsAsFactors=FALSE)

    n.loci <- (dim(ped)[2] - 6) / 2
    
    n.ind <- dim(ped)[1]

    if(header==F){names(ped)[5+2*(1:n.loci)] <- paste("SNP",1:n.loci)}

    marker.name <- names(ped)[5+2*(1:n.loci)]

  }

  else if(is.data.frame(ped)) {
    
    n.loci <- (dim(ped)[2] - 6) / 2
    
    n.ind <- dim(ped)[1]
    
    names(ped)[5+2*(1:n.loci)] <- paste("SNP",1:n.loci)
    
    marker.name <- names(ped)[5+2*(1:n.loci)] # read in marker names from ped

  }

  else if(is.matrix(ped)) {
    
    n.loci <- (dim(ped)[2] - 6) / 2
    
    n.ind <- dim(ped)[1]
    
    if(!is.null(dimnames(ped)))
      
    marker.name <- dimnames(ped)[[2]][5+2*(1:n.loci)]

  }
 
  
  if(is.na(status_missing)){
    
    ped[which(ped[,6] == 1),6] <-0
    
    ped[which(ped[,6] == 2),6] <-1  		
    
  }

  else if (status_missing!=0){
    
    ped[which(ped[,6] == 1),6] <-0
    
    ped[which(ped[,6] == 2),6] <-1
    
  }

  else if(status_missing==0){
    
    ped[which(ped[,6] == 0),6] <-9
    
    ped[which(ped[,6] == 1),6] <-0
    
    ped[which(ped[,6] == 2),6] <-1
    
    status_missing <- 9
    
  }
  

  # read in loci file if it is provided

  if (is.character(loci)) {
    
    locifile <- loci
    
    loci.raw <- readLines(locifile)
    
    loci.raw2 <- strsplit(loci.raw, ",")
   
    for(i in 1:n.loci){
        
    marker.name[i] <- loci.raw2[[1]][i]
      
  }
    
 }


  if (is.null(marker.name)) {
    
    marker.name <- paste("SNP",1:n.loci)

  }

  
  
  # Output 
  Tstat <- matrix(0,ncol=6,nrow=n.loci)
  
  dimnames(Tstat) <- list(marker.name, c("LRT0","LRT1","LRT2","Z0","Z1","Z2"))
  
  pvalue <- matrix(0,ncol=8,nrow=n.loci)
  
  dimnames(pvalue) <- list(marker.name, c("LRT0","LRT0b","LRT1","LRT2","LRT2b","Z0","Z1","Z2"))
 
  Estimates_H1 <- matrix(0,ncol=3,nrow=n.loci)

  dimnames(Estimates_H1) <- list(marker.name, c("pm","pf","rho"))
  
  Estimates_H01 <- matrix(0,ncol=2,nrow=n.loci)
  
  dimnames(Estimates_H01) <- list(marker.name, c("p","rho"))
  
  Estimates_H02 <- matrix(0,ncol=2,nrow=n.loci)
  
  dimnames(Estimates_H02) <- list(marker.name, c("pm","pf"))
  
  Estimates_H0 <- matrix(0,ncol=1,nrow=n.loci)
  
  dimnames(Estimates_H0) <- list(marker.name, c("p"))



  family.id <- unique(ped[,1])

  n.family <- length(family.id)

  select.fou <- which(ped[,3] == 0 & ped[,4] == 0 & ped[,7] != 0)

  ped <- ped[select.fou,]
 
  for(i in 1:n.loci){
    
    j <- 2*i+5

    result1 <- result(ped[,c(1:6,j:(j+1))],loci,dv,start.rho,simuno,status_missing,allele_missing,header,itertime,marker.name[i])

    Tstat[i,1] <- result1[[1]]$LRT.test
    
    Tstat[i,2] <- result1[[1]]$LRT.test2
    
    Tstat[i,3] <- result1[[1]]$LRT.test1

    Tstat[i,4] <- result1[[1]]$z0	
    
    Tstat[i,5] <- result1[[1]]$z1
    
    Tstat[i,6] <- result1[[1]]$z2


    pvalue[i,1] <- result1[[2]]$Pvalue
    
    pvalue[i,2] <- result1[[2]]$Pvalue.boot
    
    pvalue[i,3] <- result1[[2]]$Pvalue2
    
    pvalue[i,4] <- result1[[2]]$Pvalue1
    
    pvalue[i,5] <- result1[[2]]$Pvalue1.boot
    
    pvalue[i,6] <- result1[[2]]$Pvalue.z0
    
    pvalue[i,7] <- result1[[2]]$Pvalue.z1
    
    pvalue[i,8] <- result1[[2]]$Pvalue.z2
      
   
    Estimates_H1[i,1] <- result1[[3]]$pm 

    Estimates_H1[i,2] <- result1[[3]]$pf1	
    
    Estimates_H1[i,3] <- result1[[3]]$rho1  
    


    Estimates_H01[i,1] <- result1[[3]]$p.01    
    
    Estimates_H01[i,2] <- result1[[3]]$rho.01  
    


    Estimates_H02[i,1] <- result1[[3]]$pm.02  
    
    Estimates_H02[i,2] <- result1[[3]]$pf.02 


    Estimates_H0[i,1] <- result1[[3]]$p.0

  }

  pvalue1 <- noquote(format(pvalue, scientific = TRUE, digits=5))
  
  results <- list(Tstat=Tstat,Pvalue=pvalue1,Estimates_H1=Estimates_H1,Estimates_H01=Estimates_H01,Estimates_H02=Estimates_H02,Estimates_H0=Estimates_H0)
  
  if (missing(filename)!=T){
		
	sink(filename)
		
        print(results)
		
        sink()

  }

  print(results)


}
