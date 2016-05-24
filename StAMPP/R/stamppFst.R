####################################
#
# Fixation Index (Fst)
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
###################################

stamppFst <-
function (geno, nboots=100, percent=95, nclusters=1){
  
  if(class(geno)=="genlight"){  #if genotype object is a genlight object convert to a data.frame
        
    geno2 <- geno
    
    geno <- as.matrix(geno2) #extract genotype data from genlight object
    sample <- row.names(geno) #individual names
    pop.names <- pop(geno2) #population names
    ploidy <- ploidy(geno2) #ploidy level
    geno=geno*(1/ploidy) #convert genotype data (number of allele A) to precentage allele frequency
    geno[is.na(geno)]=NaN
    format <- vector(length=length(geno[,1])) 
    format[1:length(geno[,1])]="genlight"
    
    
    pops <- unique(pop.names) #population names
    
    pop.num <- vector(length=length(pops)) #create vector of population ID numbers
    
    for (i in 1:length(geno[,1])){
      pop.num[i]=which(pop.names[i]==pops) #assign population ID numbers to individuals
    }
    
    genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, ploidy, format))
    
    geno <- cbind(genoLHS, geno) #combine genotype data with labels to form stampp geno file
    
    geno[,2]=as.character(pop.names)
    geno[,4]=geno2@ploidy
    
    row.names(geno)=NULL
    
  } 
  
  geno <- geno[,-5] #remove format information  
    
  cl <- makeCluster(nclusters)
  registerDoParallel(cl) #establish clusters for multithreaded computing 
  
  percent <- ((100-percent)/100)/2  
  
  d=1
  res <- numeric(length=nboots)
  
  totalind <- nrow(geno) #number of individuals
  nloc <- ncol(geno)-4 #number of loci/markers
  pops <- unique(geno[,2]) #population IDs
  npops <- length(pops) #number of populations
  
  simple.geno <- geno
  
  p.a <- ph.a <- oh <-ninds <- matrix(NA, ncol=nloc, nrow=npops, dimnames=list(pops)) #matrix to store counts of allele A, level of hetrozygousity and number of individuals without missing data per loci
  
  for (i in 1:npops){ #loop for each population
        
    refpop <- pops[i]
    
    pop.geno <- subset(simple.geno, simple.geno[,2]==refpop)
    
    if(nloc>1){ #if there is 2 or more loci
    
      ninds[i,] <- colSums((pop.geno[,5:(4+nloc)]!="NaN")*(pop.geno[,4]/2), na.rm=TRUE) #number of diploid eqiv genomes at each loci without missing data 
      
      p.a[i,] <- colMeans(pop.geno[5:(4+nloc)], na.rm=TRUE) #total number of alleles in population at each loci multipled by number of diploid eqiv genomes
      
      oh[i,] <- colSums((((pop.geno[,5:(4+nloc)])*(pop.geno[,5:(4+nloc)]<=0.5)) + ((pop.geno[5:(4+nloc)]-0.5)*(pop.geno[5:(4+nloc)]>0.5 & pop.geno[5:(4+nloc)]!=1)))*(pop.geno[,4]/2), na.rm=TRUE) #level of hetrozygousity in population multipled by number of diploid eqiv genomes
      oh[i,] <- oh[i,]*2
      
    }else{ #if the genotype dataset only has 1 locus
      
      ninds[i,] <- sum((pop.geno[,5:(4+nloc)]!="NaN")*(pop.geno[,4]/2), na.rm=TRUE) #number of diploid eqiv genomes at each loci without missing data 
      
      p.a[i,] <- colMeans(pop.geno[5:(4+nloc)], na.rm=TRUE) #total number of alleles in population at each loci multipled by number of diploid eqiv genomes
      
      oh[i,] <- sum((((pop.geno[,5:(4+nloc)])*(pop.geno[,5:(4+nloc)]<=0.5)) + ((pop.geno[5:(4+nloc)]-0.5)*(pop.geno[5:(4+nloc)]>0.5 & pop.geno[5:(4+nloc)]!=1)))*(pop.geno[,4]/2), na.rm=TRUE) #level of hetrozygousity in population multipled by number of diploid eqiv genomes
      oh[i,] <- oh[i,]*2
      
    }
    
  }
  

  p <- p.a #frequency of allele A based on the number of diploid eqiv genomes at each loci without missing data
  oh <- oh/ninds #frequency of hetrozygous genotypes at each loci based on the number of diploid eqiv genomes without missing data

  
  index2 <- index1 <- NULL
  step <- 2
  step2 <- 1
  for(i in step:npops){
    index1 <- c(index1, i:npops)
    index2 <- c(index2, rep(step2, length(step:npops)))
    step=step+1
    step2=step2+1
  }
  
  #### Fst calculation ####
  
  r=2 
    
  ninds.dup.1 <- ninds[index1,]
  ninds.dup.2 <- ninds[index2,]
  p1 <- p[index1,]
  p2 <- p[index2,]
  oh1 <- oh[index1,]
  oh2 <- oh[index2,]
     
  n.bar <- (ninds.dup.1+ninds.dup.2)/r    
  nc <- (r*n.bar)-(((ninds.dup.1^2)+(ninds.dup.2^2)) / (r*n.bar))   
  p.bar <- ((ninds.dup.1*p1)/(r*n.bar)) + ((ninds.dup.2*p2)/(r*n.bar))  
  s.square <- ((ninds.dup.1*((p1-p.bar)^2))/n.bar) + ((ninds.dup.2*((p2-p.bar)^2))/n.bar)     
  h.bar <- ((ninds.dup.1*oh1)/(r*n.bar)) + ((ninds.dup.2*oh2)/(r*n.bar))
      
  a <- (n.bar/nc) * (s.square - (1/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - ((1/4)*h.bar) ))     
  b <- (n.bar/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - (((2*n.bar-1)/(4*n.bar))*h.bar))
  c <- (1/2)*h.bar
      
  index.inf <- which(!is.finite(a) | !is.finite(b) | !is.finite(c)) #index of infinite values
  a[index.inf]=NA #remove infinite values
  b[index.inf]=NA #remove infinite values
  c[index.inf]=NA #remove infinite values
  
  
  if(nloc>1){ #if there is more than 1 locus in the genotype dataset
  
    if(npops>2){ #if there are greater than 2 populations, ie. greater than 1 pairwise comparision 
      
      fst <- rowSums(a, na.rm=TRUE) / (rowSums(a, na.rm=TRUE) + rowSums(b, na.rm=TRUE) + rowSums(c, na.rm=TRUE)) #Fst results 
      
    }else{ #if there is only 2 populations, ie. 1 pairwise comparision
      
      fst <- sum(a, na.rm=TRUE) / (sum(a, na.rm=TRUE) + sum(b, na.rm=TRUE) + sum(c, na.rm=TRUE)) #Fst results 
      
    }
  
  }else{ # if there is only 1 locus in the genotype dataset
    
    fst <- a/(a+b+c) #Fst results
    
  }
  
  fstmat <- matrix(NA, nrow=npops, ncol=npops, dimnames=list(pops, pops)) 
  fstmat[lower.tri(fstmat)]=fst #assign Fst values to matrix of pairwise comparisions
  
  #### end of Fst calculation ####
  ################################
  #### Bootstrapped Fst calc  ####
  
  if(nboots>1){
  
    boots <- matrix(NA, ncol=(nboots+4), nrow=((npops*npops)-npops)/2) #matrix to store Fst results from bootstrapping
    pvalues <- matrix(NA, ncol=npops, nrow=npops, dimnames=list(pops, pops)) #matrix to store p-values of pairwise Fsts
    
    name.index <- matrix(NA, ncol=2, nrow=((npops*npops)-npops)/2)
    
    line=0
    
    for(i in 1:(npops-1)){
      for(y in (i+1):npops){
        line=line+1
        name.index[line,c(1,2)]=c(i,y) #all pairwise population names
      }
    }
    
    boot.names1 <- pops[name.index[,1]]
    boot.names2 <- pops[name.index[,2]]
    boot.names <- cbind.data.frame(boot.names1, boot.names2)
    
    popnames <- array(unique(geno[,2]))
    
    
    p.master <- p
    oh.master <- oh
    ninds.master <- ninds
    
    res <- foreach(d = 1:nboots, .combine=cbind) %dopar% { #loop for the specified number of bootstraps
      
      boot.index <- sample((1:nloc), nloc, replace=TRUE) #bootstrap loci names
      
      p <- p.master[,boot.index] #frequency of allele A based on bootstrapped loci
      oh <- oh.master[,boot.index] #frequecny of heterozygousity based on bootstrapped loci
      ninds <- ninds.master[,boot.index] #number of inds without missing data based on bootstrapped loci
      
      r=2 
      
      ninds.dup.1 <- ninds[index1,]
      ninds.dup.2 <- ninds[index2,]
      p1 <- p[index1,]
      p2 <- p[index2,]
      oh1 <- oh[index1,]
      oh2 <- oh[index2,]    
      
      n.bar <- (ninds.dup.1+ninds.dup.2)/r   
      nc <- (r*n.bar)-(((ninds.dup.1^2)+(ninds.dup.2^2)) / (r*n.bar)) 
      p.bar <- ((ninds.dup.1*p1)/(r*n.bar)) + ((ninds.dup.2*p2)/(r*n.bar))
      s.square <- ((ninds.dup.1*((p1-p.bar)^2))/n.bar) + ((ninds.dup.2*((p2-p.bar)^2))/n.bar)
      h.bar <- ((ninds.dup.1*oh1)/(r*n.bar)) + ((ninds.dup.2*oh2)/(r*n.bar))
      a <- (n.bar/nc) * (s.square - (1/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - ((1/4)*h.bar) ))
      b <- (n.bar/(n.bar-1)) * ( (p.bar*(1-p.bar)) - (((r-1)/r)*s.square) - (((2*n.bar-1)/(4*n.bar))*h.bar))
      c <- (1/2)*h.bar
      
      index.inf <- which(!is.finite(a) | !is.finite(b) | !is.finite(c)) #index of infinite values
      
      a[index.inf]=NA #remove infinite values
      b[index.inf]=NA #remove infinite values
      c[index.inf]=NA #remove infinite values
      
      if(nloc>1){ #if there is more than 1 locus in the genotype dataset
      
        if(npops>2){ #if there are greater than 2 populations, ie. greater than 1 pairwise comparision     
          
          rowSums(a, na.rm=TRUE) / (rowSums(a, na.rm=TRUE) + rowSums(b, na.rm=TRUE) + rowSums(c, na.rm=TRUE))  #bootstrapped Fst results    
          
        }else{#if there is only 2 populations, ie. 1 pairwise comparision
          
          sum(a, na.rm=TRUE) / (sum(a, na.rm=TRUE) + sum(b, na.rm=TRUE) + sum(c, na.rm=TRUE))  #bootstrapped Fst results    
          
        }
      
      }else{ #if there is only 1 locus in the genotype dataset
        
        a/(a+b+c) #bootstrapped Fst results 
        
      }
      
    }
    
    
    boots[,(1:nboots)]=res #add bootstrapped results to matrix
    
    boots <- cbind(boot.names, boots) #add pairwise population names to the matrix of bootstrap results
    bootrows <- length(boots[,1]) #number of pairwise comparisions
    
    order.boots <- boots[,3:(2+nboots)] 
    order.boots <- data.matrix(order.boots)
    order.boots <- t(apply(order.boots, 1, sort)) #sort bootstrapped Fst values from smallest to largest
    lowerper <- order.boots[,(ceiling(percent*nboots))] #identify the Fst values on the lower CI boundry
    upperper <- order.boots[,(floor((1-percent)*nboots))] #identify the Fst values on the upper CI boundry
    pval <- (rowSums(order.boots<=0, na.rm=TRUE))/nboots #calculate p-values
    
    boots[,3:(2+nboots)]=order.boots
    boots[,(3+nboots):(5+nboots)]=cbind(lowerper, upperper, pval)
    
    
    line=0
    
    for(i in 1:(npops-1)){ 
      for(y in (i+1):npops){
        
        #store the calculated p-value results in a matrix format 
        
        line=line+1
        pvalues[y,i]=boots[line,(nboots+5)]
        fstcol <- match(boots[line,1],pops)
        fstrow <- match(boots[line,2],pops)
        boots[line, (6+nboots)]=fstmat[fstrow, fstcol]
      }
    }
    
    row.names(pvalues)=pops #label the rows of the p-value matrix with the population names
    colnames(pvalues)=pops #label the columns of the p-value matrix with the population names
    
    bootscol <- cbind("Population1", "Population2", t(c(1:nboots)), "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst") 
    colnames(boots) <- bootscol #label the columns of the matrix of bootstrap results
    
    results <- list(Fsts=fstmat, Pvalues=pvalues, Bootstraps=boots)  
  
  }else{
    results <- fstmat
  }
  
  stopCluster(cl)
    
  return(results)
  
}
