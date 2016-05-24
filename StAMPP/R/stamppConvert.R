###########################################################
#
# Import genotype data and convert it into StAMPP format
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
###########################################################

stamppConvert <-
function (genotype.file, type="csv"){
    
    if(type=="csv" | type=="r"){ #import genotype data from csv file or R workspace
    
      if(type=="csv"){
        geno <- read.csv(genotype.file)  #import from csv file
      }else{
        geno <- genotype.file #import from R workspace
      }  
      
      
      totalind <- nrow(geno) #number of individuals
      nloc <- ncol(geno)-4 #number of loci/markers
      
      pops <- unique(geno[,2]) #population names
      npops <- length(pops) #number of populations
      
      pop.num <- vector(length=totalind) #create vector of population ID numbers
      
      for (i in 1:totalind){
        pop.num[i]=which(geno[i,2]==pops) #assign population ID numbers to individuals
      }
      
      format <- geno[,4] #genotype format
      
      ploidy <- geno[,3] #ploidy levels
      
      geno <- cbind(geno[,1:2], pop.num, ploidy, format, geno[,5:(4+nloc)]) #combine genotype data with labels to form stampp geno file
  
      ab.geno <- subset(geno, geno[,5]=="BiA") #subset individuals with AB coded genotypes
      nind.ab.geno <- length(ab.geno[,2])
      
      freq.geno <- subset(geno, geno[,5]=="freq") #subset individuals with genotypes stored as allele frequencies
      nind.freq.geno <- length(freq.geno[,2])
      
      if(nind.ab.geno > 0){
      
        tmp <- ab.geno[,-c(1:5)]
        tmp <- gsub("-9", "", as.matrix(tmp), fixed=TRUE)
        tmp.a <- gsub("B", "", as.matrix(tmp), fixed=TRUE)
        tmp.a <- nchar(as.matrix(tmp.a))
        tmp.b <- gsub("A", "", as.matrix(tmp), fixed=TRUE)
        tmp.b <- nchar(as.matrix(tmp.b))
        
        res <- matrix(NA, nrow=nind.ab.geno, ncol=nloc)
        
        for(i in 1:nloc){
          res[,i]=(tmp.a[,i]/(tmp.a[,i]+tmp.b[,i]))
        }
        
        rm(tmp.a, tmp.b, tmp)
        
        ab.geno.pt1 <- as.data.frame(ab.geno[,c(1:5)])
        ab.geno.pt2 <- as.data.frame(res)
        ab.geno <- cbind(ab.geno.pt1, ab.geno.pt2)    
        
        rm(ab.geno.pt2, ab.geno.pt1, res)
      }
      
      
      colnames(ab.geno)=colnames(geno)
      
      freq.geno.pt1 <- freq.geno[,c(1:5)]
      freq.geno.pt2 <- as.matrix(freq.geno[,-c(1:5)])
      class(freq.geno.pt2)="numeric"
      freq.geno.pt2[freq.geno.pt2==-9]=NA
      freq.geno <- cbind(freq.geno.pt1, freq.geno.pt2)
      colnames(freq.geno)=colnames(geno)
      
      comb.geno <- rbind(ab.geno, freq.geno)
      rm(ab.geno, freq.geno, geno)
      
      comb.geno[,1]=as.character(comb.geno[,1])
      comb.geno[,2]=as.character(comb.geno[,2])
      comb.geno[,3]=as.integer(as.character(comb.geno[,3]))
      comb.geno[,4]=as.integer(as.character(comb.geno[,4]))
      comb.geno[,5]=as.character(comb.geno[,5])
      
      comb.geno <- comb.geno[ order(comb.geno[,3]),]
                
      return(comb.geno)
    
    }
    
    if(type=="genlight"){
                  
      geno2 <- genotype.file
      
      geno <- as.matrix(geno2) #extract genotype data from genlight object
            
      sample <- row.names(geno) #individual names
      pop.names <- pop(geno2) #population names
      ploidy <- ploidy(geno2) #ploidy level
      geno=geno*(1/ploidy) #convert genotype data (number of allele B) to precentage allele frequency
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
      
      
      return(geno)
      
      
    }
    
    
  }
