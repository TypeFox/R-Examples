####################################
#
# Genomic Relationship Matrix
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
###################################

stamppGmatrix <-
function(geno){

  
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
  
  totalind <- length(geno[,1]) #number of individuals
  
  simple.geno <- 1-geno[5:(length(geno[1,]))] #inverse genotype data, ie from 0=BB to 0=AA
  
  p <- colMeans(simple.geno, na.rm=TRUE) #allele frequency at each loci
  
  p <- as.matrix(t(p), nrow=1, ncol=length(p))
  
  simple.geno <- simple.geno*2 #Convert allele frequency format to 0-to-2 format, ie 0=AA, 1=AB, 2=BB, 0.5=AAAB... etc.
  
  ### Calculation of genomic relationship and inbreeding coefficient ###
  
  p.dup <- p[rep(1, totalind),] #vector of allele freq of each snp in the dataset duplicated for each ind
  
  ###### stepwise calculation of genomic relationship values due to large matrix and memory restrictions ######
  
  #split.size <- floor(50/as.numeric(object.size(simple.geno)/1048576))*totalind #size of duplicated p matrix under 50mb
  
  split.size <- ceiling((50/as.numeric(object.size(simple.geno)/1048576))*totalind) #size of duplicated p matrix under 50mb
  
  
  if(split.size==0){
    splite.size=totalind #attempt split size the same as original p-matrix size --- may require computer with big memory
  }
  splits <- ceiling((totalind*totalind)/split.size) #number of 50mb splits to perform
  pre.split <- 0
  
  xj.m2p <- as.matrix(simple.geno-(2*p.dup)) #xj - 2p
  xk.m2p <- xj.m2p #xk - 2p
  
  twop.1mp  <- (2*p*(1-p)) #2p(1 - p)
  
  xj.ids <- rep(1:totalind, totalind) #duplication ids for xj
  xk.ids <- (sort((rep(1:totalind, totalind)))) #duplication ids for xk
    
  a=NULL
  
  if(split.size > 1){  #if split size is greater than 1, therefore xj,m2p.dup will be a 2dim matrix
    if(splits > 1){ #if dataset is too large and needs to be split      
      for(i in 1:(splits-1)){ #stepwise calculation of genomic relationships based on 50mb split chuncks
        
        xj.m2p.dup <- (xj.m2p[xj.ids[((pre.split*split.size)+1):(i*split.size)],])
        xk.m2p.dup <- (xk.m2p[xk.ids[((pre.split*split.size)+1):(i*split.size)],])
        
        twop.1mp.dup <- twop.1mp[rep(1, split.size),]
        
        a <- c(a, rowMeans(((xj.m2p.dup*xk.m2p.dup)/twop.1mp.dup), na.rm=TRUE)) #estimted genomic relationships
        
        pre.split <- i
        
      }
    }
    
    if( length(((pre.split*split.size)+1):(totalind*totalind))>1 ){ #if final split size is >1 and therefore calculations are on a 2dim matrix
      
      xj.m2p.dup <- (xj.m2p[xj.ids[((pre.split*split.size)+1):(totalind*totalind)],])
      xk.m2p.dup <- (xk.m2p[xk.ids[((pre.split*split.size)+1):(totalind*totalind)],])
      
      twop.1mp.dup <- twop.1mp[rep(1, length(((pre.split*split.size)+1):(totalind*totalind))),]
      
      a <- c(a, rowMeans(((xj.m2p.dup*xk.m2p.dup)/twop.1mp.dup), na.rm=TRUE)) #estimted genomic relationships
      
    }else{ #if final split size is 1 and therefore calculations are on a vector
      
      xj.m2p.dup <- (xj.m2p[xj.ids[((pre.split*split.size)+1):(totalind*totalind)],])
      xk.m2p.dup <- (xk.m2p[xk.ids[((pre.split*split.size)+1):(totalind*totalind)],])
      
      twop.1mp.dup <- twop.1mp[rep(1, length(((pre.split*split.size)+1):(totalind*totalind))),]
      
      a <- c(a, mean(((xj.m2p.dup*xk.m2p.dup)/twop.1mp.dup), na.rm=TRUE)) #estimted genomic relationships
      
    }
  }else{ #if split size is 1, therefore xj,m2p.dup will be a vector, therefore rowSums do not work
    if(splits > 1){ #if dataset is too large and needs to be split      
      for(i in 1:(splits-1)){ #stepwise calculation of genomic relationships based on 50mb split chuncks
        
        xj.m2p.dup <- (xj.m2p[xj.ids[((pre.split*split.size)+1):(i*split.size)],])
        xk.m2p.dup <- (xk.m2p[xk.ids[((pre.split*split.size)+1):(i*split.size)],])
        
        twop.1mp.dup <- twop.1mp[rep(1, split.size),]
        
        a <- c(a, mean(((xj.m2p.dup*xk.m2p.dup)/twop.1mp.dup), na.rm=TRUE)) #estimted genomic relationships
        
        pre.split <- i
        
      }
    }
    
    xj.m2p.dup <- (xj.m2p[xj.ids[((pre.split*split.size)+1):(totalind*totalind)],])
    xk.m2p.dup <- (xk.m2p[xk.ids[((pre.split*split.size)+1):(totalind*totalind)],])
    
    twop.1mp.dup <- twop.1mp[rep(1, length(((pre.split*split.size)+1):(totalind*totalind))),]
    
    a <- c(a, mean(((xj.m2p.dup*xk.m2p.dup)/twop.1mp.dup), na.rm=TRUE)) #estimted genomic relationships
        
  }
  
  
  #########################     

  a <- matrix(a, nrow=totalind, ncol=totalind)
          
  f <- ((simple.geno^2)-((1+(2*p.dup))*simple.geno)+(2*(p.dup^2))) / ((2*p.dup)*(1-p.dup)) #inbreeding coefficient (F), j=k
  f <- rowMeans(f, na.rm=TRUE)
        
  diag(a)=(1+f) #inset 1+F into the diagonal of the relationship matrix
  
  row.names(a)=geno[,1]
  
  return(a)
  
}
