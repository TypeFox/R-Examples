####################################
#
# Nei's Genetic Distance
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
###################################

stamppNeisD <-
function (geno, pop=TRUE){
    
  
    if(class(geno)=="genlight"){  #if input file is a genlight object convert to a data.frame
            
      geno2 <- geno
      
      geno <- as.matrix(geno2) #extract genotype data from genlight object
      sample <- row.names(geno) #individual names
      pop.names <- pop(geno2) #population names
      ploidy <- ploidy(geno2) #ploidy level
      geno=geno*(1/ploidy) #convert genotype data (number of allele 2) to precentage allele frequency
      geno[is.na(geno)]=NaN
      format <- vector(length=length(geno[,1])) 
      format[1:length(geno[,1])]="genlight"
      
      
      pops <- unique(pop.names) #population names
      
      pop.num <- vector(length=length(geno[,1])) #create vector of population ID numbers
      
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
    
    totalind <- nrow(geno) #number of individuals
    nloc <- ncol(geno)-4 #number of loci
    pops <- unique(geno[,2])  #population IDs
    npops <- length(pops) #number of populations
    
    if(pop==TRUE){ #if calculating genetic distance between populations
      
      p <- matrix(NA, ncol=nloc, nrow=npops, dimnames=list(pops)) #matrix to store population allele frequencies 
      
      for (i in 1:npops){
        
        #calculate population allele frequencies
        
        ref.pop <- pops[i]
        
        pop.geno <- subset(geno, geno[,2]==ref.pop)
        
        p[i,] <- colMeans(pop.geno[5:(4+nloc)], na.rm=TRUE)
        
      }
      
      colnames(p)=colnames(geno[5:(4+nloc)])
      row.names(p)=pops
      
      m <- matrix(1, nrow=dim(p)[1], ncol=dim(p)[2])
      m[is.na(p)]=0 #matrix of true false (1/0) missing genotype of each ind & loci
            
      p <- as.matrix(p)
      p[is.na(p)]=0
      
      tmp.x <- p
      tmp.y <- (1-p)*m
      jxy <- (tmp.x%*%t(tmp.x))+(tmp.y%*%t(tmp.y))
      
      ###### stepwise calculation of jxjy due to large matrix and memory restrictions ######
                  
      split.size <- ceiling((50/as.numeric(object.size(p)/1048576))*npops) #size of duplicated p matrix under 50mb
      
      if(split.size==0){
        split.size=npops #attempt split size the same as original p-matrix size --- may require computer with big memory
      }
      splits <- ceiling((npops*npops)/split.size) #number of 50mb splits to perform
      pre.split <- 0
      
      p.dup.ids <- (sort((rep(1:npops, npops)))) #duplication ids for p
      m.dup.ids <- rep(1:npops, npops) #duplication ids for m
      
      
      jxjy=NULL
      
      if(split.size > 1){ #if split size is greater than 1, therefore p.dup will be a 2dim matrix
      
      if(splits > 1){ #if dataset is too large and needs to be split      
        for(i in 1:(splits-1)){ #stepwise calculation of jxjy based on 50mb split chuncks

          p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
          m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
          m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],])      
          
          jxjy <- c(jxjy, rowSums((p.dup*m.dup)^2)+rowSums(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
          
          pre.split <- i
          
        }
      } 
      
      if( length(((pre.split*split.size)+1):(npops*npops))>1 ){ #if final split size is >1 and therefore calculations are on a 2dim matrix     
      
      p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
      m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
      m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],])
      
      jxjy <- c(jxjy, rowSums((p.dup*m.dup)^2)+rowSums(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
      
      }else{ #if final split size is 1 and therefore calculations are on a vector
      
      p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
      m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
      m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],])
        
      jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
        
      }
      
      }else{ #if split size is 1, therefore p.dup and m.dup will be a vector, therefore rowSums do not work
        
        if(splits > 1){ #if dataset is too large and needs to be split      
          for(i in 1:(splits-1)){ #stepwise calculation of jxjy based on 50mb split chuncks
            
            p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
            m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
            m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],])      
            
            jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
            
            pre.split <- i
          }
        } 
        
        p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
        m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(npops*npops)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
        m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(npops*npops)],])
        
        jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
        
      }
      
      ######################### 
            
      sq.jxjy <- sqrt(matrix(jxjy, nrow=npops, ncol=npops)) #square root of the sum of squared allele freq

      i <- jxy/sq.jxjy
      i <- i/t(sq.jxjy) #normalised identity averaged across loci
      
      d <- -log(i) #Nei's genetic distance 
      neis.d <- (matrix(as.numeric(sprintf("%.6f", d)), nrow=npops)) #summarise to six decimal places
      row.names(neis.d)=pops
      
    }
    
    
    if(pop==FALSE){ #if calculating genetic distance between individuals
      
      p <- geno[,5:(nloc+4)] #matrix of allele frequencies
      row.names(p)=geno[,1]
      colnames(p)=colnames(geno[5:(4+nloc)])
      
      m <- matrix(1, nrow=dim(p)[1], ncol=dim(p)[2])
      m[is.na(p)]=0 #matrix of true false (1/0) missing genotype of each ind & loci
      
      p <- as.matrix(p)
      p[is.na(p)]=0
      
      tmp.x <- p
      tmp.y <- (1-p)*m
      jxy <- (tmp.x%*%t(tmp.x))+(tmp.y%*%t(tmp.y))
            
      ###### stepwise calculation of jxjy due to large matrix and memory restrictions ######
          
      split.size <- ceiling((50/as.numeric(object.size(p)/1048576))*totalind) #size of duplicated p matrix under 50mb
      
      if(split.size==0){
        splite.size=totalind #attempt split size the same as original p-matrix size --- may require computer with big memory
      }
      
      splits <- ceiling((totalind*totalind)/split.size) #number of 50mb splits to perform
      
      pre.split <- 0
      
      p.dup.ids <- (sort((rep(1:totalind, totalind)))) #duplication ids for p
      m.dup.ids <- rep(1:totalind, totalind) #duplication ids for m
      
      jxjy=NULL
      
      if(split.size > 1){ #if split size is greater than 1, therefore p.dup will be a 2dim matrix
      
        if(splits > 1){ #if dataset is too large and needs to be split      
          for(i in 1:(splits-1)){ #stepwise calculation of jxjy based on 50mb split chuncks
            
            p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
            m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
            m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],])      
  
            jxjy <- c(jxjy, rowSums((p.dup*m.dup)^2)+rowSums(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
            
            pre.split <- i
            
          }
        } 
        
        if( length(((pre.split*split.size)+1):(totalind*totalind))>1 ){ #if final split size is >1 and therefore calculations are on a 2dim matrix
        
        p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
        m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
        m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],])
              
        jxjy <- c(jxjy, rowSums((p.dup*m.dup)^2)+rowSums(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
        
        }else{ #if final split size is 1 and therefore calculations are on a vector
        
        p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
        m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
        m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],])
          
        jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
                  
        }
        
        
      }else{ #if split size is 1, therefore p.dup and m.dup will be a vector, therefore rowSums do not work
        
        if(splits > 1){ #if dataset is too large and needs to be split      
          for(i in 1:(splits-1)){ #stepwise calculation of jxjy based on 50mb split chuncks
            
            p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
            m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(i*split.size)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
            m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(i*split.size)],])      
            
            jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
            
            pre.split <- i
          }
        } 
        
        p.dup <- (p[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of allele frequencies to calculate squared allele freq, jx & jy
        m.dup <- (m[m.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],]) #matrix of true false (1/0) missing genotype of each ind & loci to adjust p.dup matrix for missing alleles
        m.dup2 <- (m[p.dup.ids[((pre.split*split.size)+1):(totalind*totalind)],])
        
        jxjy <- c(jxjy, sum((p.dup*m.dup)^2)+sum(((1-p.dup)*m.dup*m.dup2)^2)) #sum of square allele freq
        
      }
     #########################     

      sq.jxjy <- sqrt(matrix(jxjy, nrow=totalind, ncol=totalind)) #square root of the sum of squared allele freq
      
      i <- jxy/sq.jxjy 
      i <- i/t(sq.jxjy) #normalised identity averaged across loci
      
      d <- -log(i) #Nei's genetic distance 
      neis.d <- (matrix(as.numeric(sprintf("%.6f", d)), nrow=totalind)) #summarise to six decimal places
      row.names(neis.d)=geno[,1]

    }
    
    return(neis.d)
    
  }
