# for cores

concatenate <- function(obj,n.cores){
    
genome       <- new("GENOME")
region.data  <- new("region.data")
region.stats <- new("region.stats")

# Init 
populations  		<- NULL
populations2 		<- NULL
popmissing   		<- NULL
outgroup     		<- NULL
outgroup2    		<- NULL
CodingSNPS   		<- NULL
Coding.matrix           <- NULL
Coding.matrix2          <- NULL 
UTRSNPS      		<- NULL
UTR.matrix              <- NULL
IntronSNPS   		<- NULL
Intron.matrix           <- NULL
ExonSNPS                <- NULL
Exon.matrix             <- NULL
GeneSNPS                <- NULL
Gene.matrix             <- NULL
reading.frame           <- NULL
rev.strand		<- NULL
transitions  		<- NULL
biallelic.matrix        <- NULL
biallelic.sites         <- NULL
biallelic.sites2        <- NULL
matrix_codonpos         <- NULL
synonymous              <- NULL
matrix_freq             <- NULL
n.singletons            <- NULL
polyallelic.sites       <- NULL
n.nucleotides           <- NULL
biallelic.compositions  <- NULL
biallelic.substitutions <- NULL
minor.alleles           <- NULL
codons                  <- as.list(NULL)
sites.with.gaps         <- NULL
sites.with.unknowns     <- as.list(NULL)

# region.stats
nucleotide.diversity   <- NULL
haplotype.diversity    <- NULL
haplotype.counts       <- NULL       # sfreqh
minor.allele.freqs     <- NULL       # JFD
biallelic.structure    <- NULL       # SXX
linkage.disequilibrium <- NULL  

n.sites             <- NULL
n.sites2            <- NULL
n.valid.sites       <- NULL    
n.gaps              <- NULL    
n.unknowns          <- NULL   
n.polyallelic.sites <- NULL
n.biallelic.sites   <- NULL
trans.transv.ratio  <- NULL
region.names        <- NULL

### Progress
progr <- progressBar()
###

for (xx in 1:n.cores) {

dat  <- obj[[xx]]@region.data 
stat <- obj[[xx]]@region.stats

# region.data

      populations                  <- c(populations,dat@populations)
      populations2     		   <- c(populations2,dat@populations2)
      popmissing                   <- c(popmissing,dat@popmissing)
      outgroup                     <- c(outgroup,dat@outgroup)
      outgroup2                    <- c(outgroup2,dat@outgroup2)

      CodingSNPS       		   <- c(CodingSNPS,dat@CodingSNPS)
      Coding.matrix                <- c(Coding.matrix,dat@Coding.matrix)
      Coding.matrix2               <- c(Coding.matrix2,dat@Coding.matrix2)
      UTRSNPS          		   <- c(UTRSNPS,dat@UTRSNPS)
      UTR.matrix                   <- c(UTR.matrix,dat@UTR.matrix)
      IntronSNPS       		   <- c(IntronSNPS,dat@IntronSNPS)
      Intron.matrix                <- c(Intron.matrix,dat@Intron.matrix)
      ExonSNPS                     <- c(ExonSNPS,dat@ExonSNPS)
      Exon.matrix                  <- c(Exon.matrix,dat@Exon.matrix)
      GeneSNPS                     <- c(GeneSNPS,dat@GeneSNPS)
      Gene.matrix                  <- c(Gene.matrix,dat@Gene.matrix)
      reading.frame                <- c(reading.frame,dat@reading.frame)
      rev.strand                   <- c(rev.strand,dat@rev.strand)

      transitions      		<- c(transitions,dat@transitions) 

      if(length(obj[[xx]]@BIG.BIAL)!=0){
	# in case of a slide or split object BIG.DATA
        bbb                  		         <- obj[[xx]]@BIG.BIAL[[1]][,obj[[xx]]@SLIDE.POS[[1]]]
	dat@biallelic.matrix    		 <- list( ff(bbb,dim=dim(bbb)) )
        colnames(dat@biallelic.matrix[[1]])      <- obj[[xx]]@region.data@biallelic.sites[[1]]
      
      } 

      biallelic.matrix 		<- c(biallelic.matrix,dat@biallelic.matrix)
      biallelic.sites  		<- c(biallelic.sites,dat@biallelic.sites)
      biallelic.sites2          <- c(biallelic.sites2,dat@biallelic.sites2) 
      matrix_codonpos 	        <- c(matrix_codonpos,dat@matrix_codonpos) 
      synonymous       		<- c(synonymous,dat@synonymous)
      matrix_freq      		<- c(matrix_freq,dat@matrix_freq)
      n.singletons    		<- c(n.singletons,dat@n.singletons) 
      polyallelic.sites 	<- c(polyallelic.sites,dat@polyallelic.sites) 
      n.nucleotides   		<- c(n.nucleotides,dat@n.nucleotides) 
      biallelic.compositions    <- c(biallelic.compositions,dat@biallelic.compositions)  
      biallelic.substitutions   <- c(biallelic.substitutions,dat@biallelic.substitutions) 
      minor.alleles    		<- c(minor.alleles,dat@minor.alleles) 
      sites.with.gaps  		<- c(sites.with.gaps,dat@sites.with.gaps) 
      sites.with.unknowns      <- c(sites.with.unknowns,dat@sites.with.unknowns)

      nucleotide.diversity   <- c(nucleotide.diversity,stat@nucleotide.diversity)
      haplotype.diversity    <- c(haplotype.diversity,stat@haplotype.diversity)
      haplotype.counts       <- c(haplotype.counts,stat@haplotype.counts)      
      minor.allele.freqs     <- c(minor.allele.freqs,stat@minor.allele.freqs)       
      biallelic.structure    <- c(biallelic.structure,stat@biallelic.structure)       
      linkage.disequilibrium <- c(linkage.disequilibrium,stat@linkage.disequilibrium)

      # GENOME data
      
      n.sites             <- c(n.sites,obj[[xx]]@n.sites)
      n.sites2            <- c(n.sites2,obj[[xx]]@n.sites)
      n.valid.sites       <- c(n.valid.sites,obj[[xx]]@n.valid.sites)    
      n.gaps              <- c(n.gaps,obj[[xx]]@n.gaps)     
      n.unknowns          <- c(n.unknowns,obj[[xx]]@n.unknowns)   
      n.polyallelic.sites <- c(n.polyallelic.sites,obj[[xx]]@n.polyallelic.sites)
      n.biallelic.sites   <- c(n.biallelic.sites,obj[[xx]]@n.biallelic.sites)
      trans.transv.ratio  <- c(trans.transv.ratio,obj[[xx]]@trans.transv.ratio)
      region.names        <- c(region.names,obj[[xx]]@region.names) 	

## Progress
progr <- progressBar(xx,n.cores, progr)
####

}

# FILL THE NEW OBJECT OF CLASS GENOME

region.data@populations  		<- populations
region.data@populations2 		<- populations2 
region.data@popmissing   		<- popmissing
region.data@outgroup     		<- outgroup 
region.data@outgroup2    		<- outgroup2 
region.data@CodingSNPS   		<- CodingSNPS
region.data@Coding.matrix               <- Coding.matrix
region.data@Coding.matrix2              <- Coding.matrix2
region.data@UTRSNPS      		<- UTRSNPS
region.data@UTR.matrix                  <- UTR.matrix
region.data@IntronSNPS   		<- IntronSNPS
region.data@Intron.matrix               <- Intron.matrix
region.data@GeneSNPS                    <- GeneSNPS
region.data@Gene.matrix                 <- Gene.matrix
region.data@ExonSNPS                    <- ExonSNPS
region.data@Exon.matrix                 <- Exon.matrix
region.data@reading.frame               <- reading.frame
region.data@rev.strand                  <- rev.strand
region.data@transitions  		<- transitions
region.data@biallelic.matrix            <- biallelic.matrix 
region.data@biallelic.sites             <- biallelic.sites  
region.data@biallelic.sites2            <- biallelic.sites2
region.data@matrix_codonpos             <- matrix_codonpos 
region.data@synonymous                  <- synonymous 
region.data@matrix_freq                 <- matrix_freq
region.data@n.singletons                <- n.singletons 
region.data@polyallelic.sites           <- polyallelic.sites 
region.data@n.nucleotides               <- n.nucleotides
region.data@biallelic.compositions      <- biallelic.compositions
region.data@biallelic.substitutions     <- biallelic.substitutions
region.data@minor.alleles               <- minor.alleles
region.data@codons                      <- codons
region.data@sites.with.gaps             <- sites.with.gaps
region.data@sites.with.unknowns         <- sites.with.unknowns 

# region.stats
region.stats@nucleotide.diversity   <- nucleotide.diversity
region.stats@haplotype.diversity    <- haplotype.diversity 
region.stats@haplotype.counts       <- haplotype.counts         # sfreqh
region.stats@minor.allele.freqs     <- minor.allele.freqs       # JFD
region.stats@biallelic.structure    <- biallelic.structure      # SXX
region.stats@linkage.disequilibrium <- linkage.disequilibrium  

genome@region.stats <- region.stats
genome@region.data  <- region.data 


genome@populations         <- vector("list",1)
genome@n.sites             <- n.sites
genome@n.sites2            <- n.sites2
genome@n.valid.sites       <- n.valid.sites  
genome@n.gaps              <- n.gaps  
genome@n.unknowns          <- n.unknowns    
genome@n.polyallelic.sites <- n.polyallelic.sites
genome@n.biallelic.sites   <- n.biallelic.sites
genome@trans.transv.ratio  <- trans.transv.ratio
genome@region.names        <- region.names

genome@big.data                  <- obj[[1]]@big.data
genome@snp.data                  <- obj[[1]]@snp.data
genome@gff.info                  <- obj[[1]]@gff.info

genome@Pop_Neutrality$calculated <- FALSE
genome@Pop_FSTN$calculated       <- FALSE
genome@Pop_FSTH$calculated       <- FALSE
genome@Pop_MK$calculated         <- FALSE
genome@Pop_Linkage$calculated    <- FALSE
genome@Pop_Slide$calculated      <- TRUE


genome@Pop_Detail$calculated     <- FALSE



return(genome)

}

# for whole genome ##############################################

concatenate_to_whole_genome <- function(obj,n.chunks){
    

snp.data     <- obj@snp.data[[1]]
genome       <- new("GENOME")
region.data  <- new("region.data")
region.stats <- new("region.stats")


# Init 
populations  		<- NULL
populations2 		<- NULL
popmissing   		<- NULL
outgroup     		<- NULL
outgroup2    		<- NULL
CodingSNPS   		<- NULL
UTRSNPS      		<- NULL
IntronSNPS   		<- NULL
ExonSNPS                <- NULL
GeneSNPS                <- NULL
transitions  		<- NULL
biallelic.substitutions <- NULL
biallelic.matrix        <- NULL
reading.frame           <- NULL
rev.strand              <- NULL
Coding.matrix           <- NULL
Coding.matrix2          <- NULL
Intron.matrix           <- NULL
UTR.matrix              <- NULL
Exon.matrix             <- NULL
Gene.matrix             <- NULL
biallelic.sites         <- NULL
biallelic.sites2        <- NULL
matrix_codonpos         <- NULL
synonymous              <- NULL
matrix_freq             <- NULL
n.singletons            <- NULL
polyallelic.sites       <- NULL
n.nucleotides           <- NULL
biallelic.compositions  <- NULL
biallelic.substitutions <- NULL
minor.alleles           <- NULL
codons                  <- as.list(NULL)
sites.with.gaps         <- NULL
sites.with.unknowns     <- as.list(NULL)

# region.stats
nucleotide.diversity   <- NULL
haplotype.diversity    <- NULL
haplotype.counts       <- NULL       # sfreqh
minor.allele.freqs     <- NULL       # JFD
biallelic.structure    <- NULL       # SXX
linkage.disequilibrium <- NULL  

n.sites             <- NULL
n.valid.sites       <- NULL    
n.gaps              <- NULL    
n.unknowns          <- NULL   
n.polyallelic.sites <- NULL
n.biallelic.sites   <- NULL
trans.transv.ratio  <- NULL
region.names        <- NULL


start <- 1

if(obj@big.data){

  iddd             <- which(obj@n.biallelic.sites>1)[1] # first region with data
  rows_bial        <- dim(obj@region.data@biallelic.matrix[[iddd]])[1]
  cols_bial        <- sum(obj@n.biallelic.sites)

  if(cols_bial==0){
  stop("No biallelic positions available for concatenation")
  }

  if(rows_bial*cols_bial>.Machine$integer.max){
    print("Warning: Matrix too big for ff package --> using the bigmemory package! ")
    #require(bigmemory)
    options(bigmemory.allow.dimnames=TRUE)
    biallelic.matrix <- bigmemory::filebacked.big.matrix(backingfile="BIGBIAL",ncol=cols_bial,nrow=rows_bial, type="double")
  }else{
    biallelic.matrix <- ff(0,dim=c(rows_bial,cols_bial))
  }

  rownames(biallelic.matrix) <- rownames(obj@region.data@biallelic.matrix[[iddd]])

  # init the gff informations
  if(obj@gff.info){
   if(obj@snp.data){
   # reading.frame + rev.strand
    cols.reading     <- sapply(obj@region.data@reading.frame,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.reading        <- sum(cols.reading)
    if(R.reading>0){
    reading.frame    <- ff(0,dim=c(R.reading,3))
    rev.strand       <- ff(0,R.reading)
    start.reading    <- 1
    }
   }

   # Coding CDS
    cols.coding      <- sapply(obj@region.data@Coding.matrix,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.coding         <- sum(cols.coding)
    if(R.coding>0){
    Coding.matrix    <- ff(0,dim=c(R.coding,2))
    start.coding     <- 1
    }
   
    # Coding2 CDS
    cols.coding2      <- sapply(obj@region.data@Coding.matrix2,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.coding2         <- sum(cols.coding2)
    if(R.coding2>0){
    Coding.matrix2    <- ff(0,dim=c(R.coding2,2))
    start.coding2     <- 1
    }

   # Exons
    cols.exons       <- sapply(obj@region.data@Exon.matrix,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.exons          <- sum(cols.exons)
    if(R.exons>0){
    Exon.matrix      <- ff(0,dim=c(R.exons,2))
    start.exons      <- 1
    }
    # Introns
    cols.introns     <- sapply(obj@region.data@Intron.matrix,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.introns        <- sum(cols.introns)
    if(R.introns){
    Intron.matrix    <- ff(0,dim=c(R.introns,2))
    start.introns    <- 1
    }
    # UTR
    cols.utrs        <- sapply(obj@region.data@UTR.matrix,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.utrs           <- sum(cols.utrs)
    if(R.utrs){
    UTR.matrix       <- ff(0,dim=c(R.utrs,2))
    start.utrs       <- 1
    }
    # UTR
    cols.genes       <- sapply(obj@region.data@Gene.matrix,function(x){if(length(x)==0){return(0)};return(dim(x)[1])})
    R.genes          <- sum(cols.genes)
    if(R.genes){
    Gene.matrix      <- ff(0,dim=c(R.genes,2))
    start.genes      <- 1
    }

  }
  
}

ADD.GFF  <- 0
ADD.GFF2 <- 0
ADD.GFF3 <- 0
ADD.GFF4 <- 0
ADD.GFF5 <- 0
ADD.GFF6 <- 0

# Important for GFF
if(snp.data){ 
  ADD.SITES <- obj@n.sites2
}else{
  ADD.SITES <- obj@n.sites
} 


cat("\n")
cat("Concatenate regions \n")
### Progress
 progr <- progressBar()
###




for (xx in 1:n.chunks) {

dat       <- obj@region.data
stat      <- obj@region.stats

if(obj@n.biallelic.sites[xx]==0){next}

# region.data

      # populations                  <- c(populations,dat@populations[[xx]])
      # populations2                 <- c(populations2,dat@populations2)
      # popmissing                   <- c(popmissing,dat@popmissing)
      # outgroup                     <- c(outgroup,dat@outgroup)
      # outgroup2                    <- c(outgroup2,dat@outgroup2)

        CodingSNPS       		   <- c(CodingSNPS,dat@CodingSNPS[[xx]])
        UTRSNPS          		   <- c(UTRSNPS,dat@UTRSNPS[[xx]])
        IntronSNPS       		   <- c(IntronSNPS,dat@IntronSNPS[[xx]])
        ExonSNPS                           <- c(ExonSNPS,dat@ExonSNPS[[xx]])
        GeneSNPS                           <- c(GeneSNPS,dat@GeneSNPS[[xx]])
      
      transitions      		<- c(transitions,dat@transitions[[xx]]) 
      synonymous                <- c(synonymous,dat@synonymous[[xx]])
      biallelic.substitutions   <- cbind(biallelic.substitutions,dat@biallelic.substitutions[[xx]])
       
      if(obj@big.data){

       end   <- start + (obj@n.biallelic.sites[xx]-1) 
       open(dat@biallelic.matrix[[xx]]) 
       biallelic.matrix[,start:end]     <- dat@biallelic.matrix[[xx]][,]
       close(dat@biallelic.matrix[[xx]])
       start <- end + 1
       
       ### GFF Infos
       if(obj@gff.info){
	  
        if(obj@snp.data){
          # reading.frame + rev.strand
         if(obj@snp.data){
          if(R.reading){
           if(cols.reading[xx]!=0){
            end.reading                               <- start.reading + cols.reading[xx] - 1
            open(dat@reading.frame[[xx]]) 
            open(dat@rev.strand[[xx]])      
            reading.frame[start.reading:end.reading,] <- dat@reading.frame[[xx]][,] # + ADD.GFF6
            rev.strand[start.reading:end.reading]     <- dat@rev.strand[[xx]][,]
            close(dat@reading.frame[[xx]])
            close(dat@rev.strand[[xx]]) 
            start.reading  <- end.reading + 1
            # ADD.GFF      <- end
            # ADD.GFF6       <- ADD.GFF6 + obj@n.sites[xx]
           }
          }
         }

          # Coding CDS
          if(R.coding){
           if(cols.coding[xx]!=0){
            end.coding                              <- start.coding + cols.coding[xx] - 1 
            open(dat@Coding.matrix[[xx]])       
            Coding.matrix[start.coding:end.coding,] <- dat@Coding.matrix[[xx]][,] #+ ADD.GFF
	    close(dat@Coding.matrix[[xx]]) 
            start.coding <- end.coding + 1
            #ADD.GFF      <- ADD.GFF + ADD.SITES[xx]
           }
          }
           # Coding2 CDS
          if(R.coding2){
           if(cols.coding2[xx]!=0){
            end.coding2                                <- start.coding2 + cols.coding2[xx] - 1
            open(dat@Coding.matrix2[[xx]]) 
            if(xx>1){       
            Coding.matrix2[start.coding2:end.coding2,] <- dat@Coding.matrix2[[xx]][,] + sum(ADD.SITES[1:(xx-1)]) #+ ADD.GFF
            }else{
	    Coding.matrix2[start.coding2:end.coding2,] <- dat@Coding.matrix2[[xx]][,]	
	    }
            close(dat@Coding.matrix2[[xx]])
            start.coding2 <- end.coding2 + 1           
           }
           #ADD.GFF       <- ADD.GFF + ADD.SITES[xx]
          }

          if(R.exons){
           if(cols.exons[xx]!=0){
            # Exon
            end.exons                               <- start.exons + cols.exons[xx] - 1 
            open(dat@Exon.matrix[[xx]])       
            Exon.matrix[start.exons:end.exons,]     <- dat@Exon.matrix[[xx]][,] #+ ADD.GFF2
            close(dat@Exon.matrix[[xx]]) 
            start.exons <- end.exons + 1
            # ADD.GFF2 <- end
            # ADD.GFF2      <- ADD.GFF2 + obj@n.sites[xx]
            #ADD.GFF2      <- ADD.GFF2 + ADD.SITES[xx]
           }
          }
          if(R.introns){
           if(cols.introns[xx]!=0){
            end.introns                                   <- start.introns + cols.introns[xx] - 1 
            # Intron
            open(dat@Intron.matrix[[xx]]) 
            Intron.matrix[start.introns:end.introns,]     <- dat@Intron.matrix[[xx]][,] #+ ADD.GFF3
            close(dat@Intron.matrix[[xx]]) 
            start.introns <- end.introns + 1
            #ADD.GFF3      <- end
            #ADD.GFF3        <- ADD.GFF3 + obj@n.sites[xx]
            #ADD.GFF3      <- ADD.GFF3  + ADD.SITES[xx]
           }
          }
          if(R.utrs){
           if(cols.utrs[xx]!=0){
            # UTR
            end.utrs                               <- start.utrs + cols.utrs[xx] - 1 
            open(dat@UTR.matrix[[xx]]) 
            UTR.matrix[start.utrs:end.utrs,]       <- dat@UTR.matrix[[xx]][,] #+ ADD.GFF4
            close(dat@UTR.matrix[[xx]]) 
            start.utrs <- end.utrs + 1
            # ADD.GFF4   <- end
            # ADD.GFF4      <- ADD.GFF4 + obj@n.sites[xx]
            #ADD.GFF4   <- ADD.GFF4 + ADD.SITES[xx]
           }
          }
          if(R.genes){
           if(cols.genes[xx]!=0){
            # Gene
            end.genes                               <- start.genes + cols.genes[xx] - 1
            open(dat@Gene.matrix[[xx]])  
            Gene.matrix[start.genes:end.genes,]     <- dat@Gene.matrix[[xx]][,] # + ADD.GFF5
            close(dat@Gene.matrix[[xx]]) 
            start.genes <- end.genes + 1
            # ADD.GFF5    <- end
            # ADD.GFF5     <- ADD.GFF5 + obj@n.sites[xx]
            #ADD.GFF5     <- ADD.GFF5 + ADD.SITES[xx]
           }
          }
        }# end of snp.data
        else{
	  # reading.frame
         if(obj@snp.data){
          if(R.reading){
           if(cols.reading[xx]!=0){
            end.reading                               <- start.reading + cols.reading[xx] - 1
            open(dat@reading.frame[[xx]])      
            open(dat@rev.strand[[xx]])  
            reading.frame[start.reading:end.reading,] <- dat@reading.frame[[xx]][,] # + ADD.GFF6
            rev.strand[start.reading:end.reading]     <- dat@rev.strand[[xx]][,] # + ADD.GFF6
            close(dat@reading.frame[[xx]])
            close(dat@rev.strand[[xx]]) 
            start.reading  <- end.reading + 1
            # ADD.GFF      <- end
            # ADD.GFF6       <- ADD.GFF6 + obj@n.sites[xx]
           }
          }
         }

          # Coding CDS
          if(R.coding){
           if(cols.coding[xx]!=0){
            end.coding                              <- start.coding + cols.coding[xx] - 1 
            open(dat@Coding.matrix[[xx]])       
            Coding.matrix[start.coding:end.coding,] <- dat@Coding.matrix[[xx]][,] + ADD.GFF
            close(dat@Coding.matrix[[xx]]) 
            start.coding <- end.coding + 1
            ADD.GFF      <- ADD.GFF + ADD.SITES[xx]
           }
          }
           # Coding2 CDS
          #if(R.coding2){
          # if(cols.coding2[xx]!=0){
          #  end.coding2                                <- start.coding2 + cols.coding2[xx] - 1       
          #  Coding.matrix2[start.coding2:end.coding2,] <- dat@Coding.matrix2[[xx]][,] + ADD.GFF
          #  start.coding2 <- end.coding2 + 1
          #  ADD.GFF       <- ADD.GFF + ADD.SITES[xx]
          # }
          #}

          if(R.exons){
           if(cols.exons[xx]!=0){
            # Exon
            end.exons                               <- start.exons + cols.exons[xx] - 1
            open(dat@Exon.matrix[[xx]])        
            Exon.matrix[start.exons:end.exons,]     <- dat@Exon.matrix[[xx]][,] + ADD.GFF2
            close(dat@Exon.matrix[[xx]]) 
            start.exons <- end.exons + 1
            # ADD.GFF2 <- end
            # ADD.GFF2      <- ADD.GFF2 + obj@n.sites[xx]
            ADD.GFF2      <- ADD.GFF2 + ADD.SITES[xx]
           }
          }
          if(R.introns){
           if(cols.introns[xx]!=0){
            end.introns                                   <- start.introns + cols.introns[xx] - 1 
            # Intron
            open(dat@Intron.matrix[[xx]]) 
            Intron.matrix[start.introns:end.introns,]     <- dat@Intron.matrix[[xx]][,] + ADD.GFF3
            close(dat@Intron.matrix[[xx]]) 
            start.introns <- end.introns + 1
            # ADD.GFF3      <- end
            # ADD.GFF3        <- ADD.GFF3 + obj@n.sites[xx]
            ADD.GFF3      <- ADD.GFF3  + ADD.SITES[xx]
           }
          }
          if(R.utrs){
           if(cols.utrs[xx]!=0){
            # UTR
            end.utrs                               <- start.utrs + cols.utrs[xx] - 1
            open(dat@UTR.matrix[[xx]])  
            UTR.matrix[start.utrs:end.utrs,]       <- dat@UTR.matrix[[xx]][,] + ADD.GFF4
            close(dat@UTR.matrix[[xx]]) 
            start.utrs <- end.utrs + 1
            # ADD.GFF4   <- end
            # ADD.GFF4      <- ADD.GFF4 + obj@n.sites[xx]
            ADD.GFF4   <- ADD.GFF4 + ADD.SITES[xx]
           }
          }
          if(R.genes){
           if(cols.genes[xx]!=0){
            # Gene
            end.genes                               <- start.genes + cols.genes[xx] - 1
            open(dat@Gene.matrix[[xx]])  
            Gene.matrix[start.genes:end.genes,]     <- dat@Gene.matrix[[xx]][,]  + ADD.GFF5
            close(dat@Gene.matrix[[xx]]) 
            start.genes <- end.genes + 1
            # ADD.GFF5    <- end
            # ADD.GFF5     <- ADD.GFF5 + obj@n.sites[xx]
            ADD.GFF5     <- ADD.GFF5 + ADD.SITES[xx]
           }
          }	
         }# end of else snp.data

       } # end of gff.info      
 
       
      }else{ # end of big.data

          end                           <- start + (obj@n.biallelic.sites[xx]-1)  
          biallelic.matrix 		<- cbind(biallelic.matrix,dat@biallelic.matrix[[xx]])
          Coding.matrix                 <- rbind(Coding.matrix,dat@Coding.matrix[[xx]]+ADD.GFF)
          Exon.matrix                   <- rbind(Exon.matrix,dat@Exon.matrix[[xx]]+ADD.GFF2)
          UTR.matrix                    <- rbind(UTR.matrix,dat@UTR.matrix[[xx]]+ADD.GFF3)
          Gene.matrix                   <- rbind(Gene.matrix,dat@Gene.matrix[[xx]]+ADD.GFF4)
          Intron.matrix                 <- rbind(Intron.matrix,dat@Intron.matrix[[xx]]+ADD.GFF5)
         
        #  ADD.GFF         <- ADD.GFF  + obj@n.sites[xx]       
        #  ADD.GFF2        <- ADD.GFF2 + obj@n.sites[xx] 
        #  ADD.GFF3        <- ADD.GFF3 + obj@n.sites[xx]       
        #  ADD.GFF4        <- ADD.GFF4 + obj@n.sites[xx]
        #  ADD.GFF5        <- ADD.GFF5 + obj@n.sites[xx]
  
	  ADD.GFF         <- ADD.GFF  + ADD.SITES[xx]       
          ADD.GFF2        <- ADD.GFF2 + ADD.SITES[xx] 
          ADD.GFF3        <- ADD.GFF3 + ADD.SITES[xx]       
          ADD.GFF4        <- ADD.GFF4 + ADD.SITES[xx]
          ADD.GFF5        <- ADD.GFF5 + ADD.SITES[xx] 
       
      }



      if(xx>1){
          biallelic.sites2  		<- c(biallelic.sites2,(dat@biallelic.sites2[[xx]] + sum(ADD.SITES[1:(xx-1)]))) #biallelic.sites2[length(biallelic.sites2)]))
          if(snp.data){
            biallelic.sites  		<- c(biallelic.sites,(dat@biallelic.sites[[xx]]))  
          }else{
            biallelic.sites  		<- c(biallelic.sites,(dat@biallelic.sites[[xx]]   + sum(ADD.SITES[1:(xx-1)])))  #biallelic.sites[length(biallelic.sites)])) 
          }
      }else{
	  biallelic.sites               <- c(biallelic.sites,dat@biallelic.sites[[xx]]) 
          biallelic.sites2              <- c(biallelic.sites2,dat@biallelic.sites2[[xx]]) 
      }


     # matrix_codonpos 	        <- c(matrix_codonpos,dat@matrix_codonpos[[xx]]) 
     # synonymous       		<- c(synonymous,dat@synonymous[[xx]])
     # matrix_freq      		<- c(matrix_freq,dat@matrix_freq[[xx]])
     # n.singletons    		<- c(n.singletons,dat@n.singletons[[xx]]) 
     # polyallelic.sites 	<- c(polyallelic.sites,dat@polyallelic.sites[[xx]]) 
     # n.nucleotides   		<- c(n.nucleotides,dat@n.nucleotides[[xx]]) 
     # biallelic.compositions    <- c(biallelic.compositions,dat@biallelic.compositions[[xx]])  
     # biallelic.substitutions   <- c(biallelic.substitutions,dat@biallelic.substitutions[[xx]]) 
     # minor.alleles             <- c(minor.alleles,dat@minor.alleles[[xx]]) 
     # sites.with.gaps  		<- c(sites.with.gaps,dat@sites.with.gaps[[xx]]) 
     # rsites.with.unknowns      <- c(sites.with.unknowns,dat@sites.with.unknowns[[xx]])

     # nucleotide.diversity   <- c(nucleotide.diversity,stat@nucleotide.diversity)
     # haplotype.diversity    <- c(haplotype.diversity,stat@haplotype.diversity)
     # haplotype.counts       <- c(haplotype.counts,stat@haplotype.counts)      
     # minor.allele.freqs     <- c(minor.allele.freqs,stat@minor.allele.freqs)       
     # biallelic.structure    <- c(biallelic.structure,stat@biallelic.structure)       
     # linkage.disequilibrium <- c(linkage.disequilibrium,stat@linkage.disequilibrium)

      # GENOME data
      
      n.sites              <- c(n.sites,obj@n.sites[xx])
      #n.valid.sites       <- c(n.valid.sites,obj@n.valid.sites[xx])    
      #n.gaps              <- c(n.gaps,obj@n.gaps[xx])     
      #n.unknowns          <- c(n.unknowns,obj@n.unknowns[xx])   
      #n.polyallelic.sites <- c(n.polyallelic.sites,obj@n.polyallelic.sites[xx])
       n.biallelic.sites   <- c(n.biallelic.sites,obj@n.biallelic.sites[xx])
       trans.transv.ratio  <- c(trans.transv.ratio,obj@trans.transv.ratio[xx])
      #region.names        <- c(region.names,obj@region.names[xx]) 	

## Progress
 progr <- progressBar(xx,n.chunks, progr)
####

}
cat("\n")

# FILL THE NEW OBJECT OF CLASS GENOME

#region.data@populations  		<- populations
#region.data@populations2 		<- populations2 
#region.data@popmissing   		<- popmissing
#region.data@outgroup     		<- outgroup 
#region.data@outgroup2    		<- outgroup2 
region.data@CodingSNPS   		<- list(CodingSNPS)
region.data@UTRSNPS      		<- list(UTRSNPS)
region.data@IntronSNPS   		<- list(IntronSNPS)
region.data@ExonSNPS                    <- list(ExonSNPS)
region.data@GeneSNPS                    <- list(GeneSNPS)

colnames(biallelic.matrix)               <- biallelic.sites
region.data@transitions                  <- list(transitions)
region.data@biallelic.substitutions      <- list(biallelic.substitutions)
region.data@biallelic.matrix             <- list(biallelic.matrix) 
region.data@biallelic.sites              <- list(biallelic.sites)
region.data@biallelic.sites2             <- list(biallelic.sites2)
region.data@synonymous                   <- list(synonymous)
region.data@Coding.matrix                <- list(Coding.matrix)
region.data@Coding.matrix2               <- list(Coding.matrix2)
region.data@Exon.matrix                  <- list(Exon.matrix)
region.data@Intron.matrix                <- list(Intron.matrix)
region.data@UTR.matrix                   <- list(UTR.matrix)
region.data@Gene.matrix                  <- list(Gene.matrix)
if(obj@snp.data){
 region.data@reading.frame               <- list(reading.frame)
 region.data@rev.strand                  <- list(rev.strand)  
}
#region.data@matrix_codonpos         <- matrix_codonpos 
#region.data@synonymous              <- synonymous 
#region.data@matrix_freq             <- matrix_freq
#region.data@n.singletons            <- n.singletons 
#region.data@polyallelic.sites       <- polyallelic.sites 
#region.data@n.nucleotides           <- n.nucleotides
#region.data@biallelic.compositions   <- biallelic.compositions
#region.data@biallelic.substitutions <- biallelic.substitutions
#region.data@minor.alleles           <- minor.alleles
#region.data@codons                  <- codons
#region.data@sites.with.gaps         <- sites.with.gaps
#region.data@sites.with.unknowns     <- sites.with.unknowns 

## region.stats
#region.stats@nucleotide.diversity   <- nucleotide.diversity
#region.stats@haplotype.diversity    <- haplotype.diversity 
#region.stats@haplotype.counts       <- haplotype.counts         # sfreqh
#region.stats@minor.allele.freqs     <- minor.allele.freqs       # JFD
#region.stats@biallelic.structure    <- biallelic.structure      # SXX
#region.stats@linkage.disequilibrium <- linkage.disequilibrium  

#genome@region.stats <- region.stats
 
 genome@region.data  <- region.data

 if(obj@snp.data){
 genome@n.sites      <- biallelic.sites[length(biallelic.sites)]
 genome@n.sites2     <- sum(obj@n.sites2,na.rm=TRUE)
 }else{
 genome@n.sites      <- sum(obj@n.sites,na.rm=TRUE)
 }

 genome@n.biallelic.sites       <- sum(obj@n.biallelic.sites,na.rm=TRUE)
 genome@n.gaps                  <- sum(obj@n.gaps,na.rm=TRUE)
 genome@n.unknowns              <- sum(obj@n.unknowns,na.rm=TRUE)
 genome@n.valid.sites           <- sum(obj@n.valid.sites,na.rm=TRUE)
 genome@n.polyallelic.sites     <- sum(obj@n.polyallelic.sites,na.rm=TRUE)
 # genome@trans.transv.ratio      <- sum(obj@trans.transv.ratio,na.rm=TRUE)
 genome@trans.transv.ratio      <- NaN
 
 
 genome@region.data@populations <- list(obj@region.data@populations[[1]]) 
 genome@region.data@outgroup    <- list(obj@region.data@outgroup[[1]])
 genome@region.data@popmissing  <- list(NULL)
 genome@populations             <- obj@populations
 genome@genelength              <- 1 

#genome@n.sites             <- n.sites
#genome@n.valid.sites       <- n.valid.sites  
#genome@n.gaps              <- n.gaps  
#genome@n.unknowns          <- n.unknowns    
#genome@n.polyallelic.sites <- n.polyallelic.sites
#genome@n.biallelic.sites   <- n.biallelic.sites
#genome@trans.transv.ratio  <- trans.transv.ratio
#genome@region.names        <- region.names

genome@region.names              <- "Concatenate"
genome@snp.data                  <- obj@snp.data
genome@big.data                  <- obj@big.data
genome@gff.info                  <- obj@gff.info
genome@Pop_Slide$calculated      <- TRUE
genome@Pop_Neutrality$calculated <- FALSE
genome@Pop_FSTN$calculated       <- FALSE
genome@Pop_Recomb$calculated     <- FALSE
genome@Pop_FSTH$calculated       <- FALSE
genome@Pop_MK$calculated         <- FALSE
genome@Pop_Linkage$calculated    <- FALSE
genome@Pop_Detail$calculated     <- FALSE




return(genome)

}



