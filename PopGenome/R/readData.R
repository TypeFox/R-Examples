
readData <- function(path,populations=FALSE,outgroup=FALSE,include.unknown=FALSE,gffpath=FALSE,format="fasta",parallized=FALSE,progress_bar_switch=TRUE,
FAST=FALSE,big.data=FALSE,SNP.DATA=FALSE){


## CHECK INPUT
if(!file.exists(path)){stop("Cannot find path !")}
if(!file.info(path)[2]){stop("Put your file/files in a folder !")}

if(gffpath[1]!=FALSE){
	if(!file.exists(gffpath)){stop("Cannot find gff path !")}
	if(!file.info(gffpath)[2]){stop("Put your file/files in a folder ! (GFF)")}
}

if(format=="HapMap"){SNP.DATA=TRUE;FAST=TRUE}
if(format=="VCF")   {SNP.DATA=TRUE;FAST=TRUE}
#if(format=="VCFhap"){SNP.DATA=TRUE;FAST=TRUE}

if(SNP.DATA){big.data=TRUE}

# Parallized version of readData

if(parallized){
 
	#library(parallel)
	n.cores <- parallel::detectCores() #multicore:::detectCores()
	#n.cores <- n.cores - 1 
	#options(cores=n.cores)
        #getOption('cores')

	files       <- list.files(path)
        split_files <- split(files,sort(rep(1:n.cores,ceiling(length(files)/n.cores))))
       
        if(gffpath[1]!=FALSE){
        gff_files       <- list.files(gffpath)
        split_files_gff <- split(gff_files,sort(rep(1:n.cores,ceiling(length(gff_files)/n.cores))))
        }
 
	
   xxx <- NULL
   yyy <- NULL

   for (i in 1:n.cores){
   
       command <- paste("mkdir split",i,sep="")      
       system(command)
       filename      <- paste("split",i,sep="")
       filename_path <- file.path(getwd(),filename)
       sapply(split_files[[i]],function(x){
       command <- file.path(path,x)
       #command <- paste("mv",command,filename_path,sep=" ")
       command <- paste("cp",command,filename_path,sep=" ")
       system(command)
       
       })
       
       xxx <- c(xxx,filename_path)
     
      if(gffpath[1]!=FALSE){

       command <- paste("mkdir GFFRObjects_split",i,sep="")      
       system(command)
       filename      <- paste("GFFRObjects_split",i,sep="")
       filename_path <- file.path(getwd(),filename)
       sapply(split_files_gff[[i]],function(x){
       command <- file.path(gffpath,x)
       #command <- paste("mv",command,filename_path,sep=" ")
        command <- paste("cp",command,filename_path,sep=" ")
       system(command)
       })
       
       yyy <- c(yyy,filename_path)

      }
 
    }


if(gffpath[1]==FALSE){

#print(xxx)
#print(format)
#print(progress_bar_switch)
#stop("Jo")

 rueck    	      <- parallel::mclapply(xxx,readData,
			 format=format,parallized=FALSE,progress_bar_switch=FALSE,
			 FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA,
			 include.unknown=include.unknown, mc.cores = n.cores, mc.silent = TRUE)

}else{

cat("Calculation ... \n")
INPUT <- vector("list",n.cores)
for(iii in 1:n.cores){
 INPUT[[iii]] <- c(xxx[iii],yyy[iii])
}

rueck    	       <- parallel::mclapply(INPUT,function(x){
                          daten    <- x[1]
                          gffdaten <- x[2]
                          TT <- readData(path=daten,
                                 gffpath=gffdaten,format=format,parallized=FALSE,progress_bar_switch=FALSE,
                                 FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA,include.unknown=include.unknown)
                          
                          return(TT) 
                          },mc.cores = n.cores, mc.silent = TRUE, mc.preschedule = TRUE)

}


genome   	       <- concatenate(rueck,n.cores)
genome@basepath        <- file.path(path)  
genome@project         <- file.path(path) 
liste    	       <- list.files(path,full.names=TRUE)
genome@genelength      <- length(liste)

if(FAST){
genome@Pop_Slide$calculated      <- TRUE # nur wegen Einschraenkungen, die auch im Slid Modus gelten
}else{
genome@Pop_Slide$calculated      <- FALSE
}

if(!is.list(populations)){genome@populations <- list(NULL)}else{genome@populations <- populations}
genome@outgroup        <- outgroup

## Move splits back to original folder
## and Delete splits afterwards #TODO

for (i in 1:n.cores){  
        #command <- paste("mv split",i,"/*"," ",path,sep="")
	#system(command)
        command <- paste("rm -r split",i,sep="")      
        system(command)
       if(gffpath[1]!=FALSE){
	#command <- paste("mv split",i,"/*"," ",gffpath,sep="")
	#system(command)
        command <- paste("rm -r GFFRObjects_split",i,sep="")      
        system(command)
       }
}# end of delete
cat("\n")
return(genome)

   # SNOWFALL

   #def.cpus <- 3
   #sfInit(parallel=TRUE,cpus=def.cpus,type="SOCK")
   #sfLibrary(ape)
   #sfExportAll() # notwendig !
   #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/GENOME.R")
   #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/GEN.R")
   #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/DATA.R")
   
   # split the Data
   #files   <- list.files(path)
   #split_files <- split(files,1:def.cpus)
   
   #xxx <- NULL
   #for (i in 1:def.cpus){
   
   #    command <- paste("md split",i,sep="")      
   #    shell(command)
   #    filename      <- paste("split",i,sep="")
   #    filename_path <- file.path(getwd(),filename)
   #    sapply(split_files[[i]],function(x){
   #      command <- file.path(path,x)
   #      command <- paste("cp",command,filename_path,sep=" ")
   #      shell(command)
   #    })
   #    
   #    xxx <- c(xxx,filename_path) 
   #}
   
   #genome <- sfLapply(xxx,readData,parallized=FALSE)
   #sfStop()
   #return(genome)
   #genome <- sfLapply()

}



#### End of parallization

methods <- "DATA"


########################
# if (.Platform$OS.type == "unix") {
# path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.so")
# }else{
# ## muss je nach 32 64 Win geaendert werden !
# path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.dll")
# }
# dyn.load(path_C_code,PACKAGE="PopGenome")
########################
  
 
 npops    <- length(populations)
 popnames <- paste("pop",1:npops)
 
 # Get the alignments
 liste    <- list.files(path,full.names=TRUE)
 liste2   <- list.files(path)
 liste3   <- gsub("\\.[a-zA-Z]*","",liste2)
 

 # sort the list -------------
 ordered            <- as.numeric(gsub("\\D", "", liste))
 if(!any(is.na(ordered))){
  names(ordered)    <- 1:length(ordered)
  ordered           <- sort(ordered)
  ids               <- as.numeric(names(ordered))
  liste             <- liste [ids]
  liste2            <- liste2[ids]
  liste3            <- liste3[ids]  
 }
 # --------------------------- 
 

 gff_objects  <- vector("list",length(liste))
 SNP.GFF      <- FALSE

 
 GFF.BOOL        <- FALSE
 ## GET THE GFF FILES -----------------------------
 if(gffpath[1]!=FALSE){ 
 
    GFF.BOOL     <- TRUE
    gff_liste    <- list.files(gffpath,full.names=TRUE)
    gff_liste2   <- list.files(gffpath)
    gff_liste3   <- gsub("\\.[a-zA-Z]*","",gff_liste2)
    # print(gff_liste3)
    # print("-------------")
    # print(liste3)             
    #treffer      <- match(gff_liste3,liste3)
    treffer       <- match(liste3,gff_liste3)
    #print(liste3)
    #print(gff_liste3)
	
    #print("------------")
    #print(treffer)
    gff_liste    <- gff_liste[treffer]
    #print(liste)
    #print(gff_liste)
    #stop("")
    
 }
 

 # ---------------------------------------------


################# Get the poppairs names ##################################
if(npops>1){
 if(outgroup[1]!=FALSE){
  poppairs <- choose(npops+1,2) # Outgroup is included !!
  pairs    <- combn(1:(npops+1),2)
 }else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 } 
###########################################################################

 
#### --- Names of population pairs --- ####################################### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ######################################### 


 
 ### --- ------------------------------ #### 
 sizeliste              <- length(liste)
 genome                 <- new("GENOME")  
 genome@basepath        <- file.path(path)  
 genome@project         <- file.path(path) 
 genome@genelength      <- sizeliste
 #populationsX           <- as.matrix(populations)
 #rownames(populationsX) <- popnames
 #colnames(populationsX) <- "Number of Samples"
 if(!is.list(populations)){genome@populations <- list(NULL)}else{genome@populations <- populations}
 genome@poppairs        <- nn
 genome@outgroup        <- outgroup
 genome@region.names    <- liste2
 DATABOOL               <- is.element("DATA",methods)

 nsites                 <- integer(sizeliste)
 nsites2                <- integer(sizeliste)

 # if(ALL_METHODS){
 # genelist  <- vector("list",length(liste))      # list of objects (class GEN )
 # datalist  <- vector("list",length(liste))    # list of objects (class DATA)
 # } 

 # INIT
 # do this calculation every time
 region.data  <- new("region.data")
 region.stats <- new("region.stats")
 init         <- vector("list",sizeliste)
 init2        <- numeric(sizeliste)
 init3        <- rep(NaN,sizeliste)

      # region.data init
      populationsX      <- init
      populations2      <- init
      popmissing        <- init 
      outgroupX         <- init
      outgroup2         <- init

      CodingSNPS        <- init
      UTRSNPS           <- init
      IntronSNPS        <- init
      ExonSNPS          <- init
      GeneSNPS          <- init
  
      reading.frame     <- init
      rev.strand        <- init 
      Coding.matrix     <- init
      Coding.matrix2    <- init
      UTR.matrix        <- init
      Intron.matrix     <- init
      Exon.matrix       <- init
      Gene.matrix       <- init

      #Coding.info       <- init
      #UTR.info          <- init
      #Intron.info       <- init
      #Exon.info         <- init
      #Gene.info         <- init
      

      Coding.region     <- init2
      UTR.region        <- init2
      Intron.region     <- init2
      Exon.region       <- init2
      Gene.region       <- init2

      transitions         <- init  # matrix_sv  transition war eine 1
      biallelic.matrix    <- init  # matrix_pol
      biallelic.sites     <- init  # matrix_pos
      biallelic.sites2    <- init
      
      reference           <- init
      matrix_codonpos     <- init # codonpos. of biallelics
      synonymous          <- init # synnonsyn
      matrix_freq         <- init
      n.singletons        <- init # unic
      polyallelic.sites   <- init # mhitbp
      n.nucleotides       <- init # sum_sam
      biallelic.compositions  <- init # TCGA
      biallelic.substitutions <- init # subst
      minor.alleles       <- init # mutations
      codons              <- init
      sites.with.gaps     <- init # gaps
      sites.with.unknowns <- init
      
      # GENOME data  init
      n.valid.sites       <- init2
      n.gaps              <- init2
      n.unknowns          <- init2
      n.polyallelic.sites <- init2
      n.biallelic.sites   <- init2
      trans.transv.ratio  <- init3


## PROGRESS #########################
if(progress_bar_switch){ ### because of parallized 
   progr <- progressBar()
}
#####################################

 # -----------------------------#
 #if(!is.list(populations)){
 #allseq <- T}else{allseq <- F}
 # -----------------------------#

 for(xx in 1:sizeliste){ # 
      
    if(!progress_bar_switch){
	print(liste[xx])
    }	

    #print(liste[xx])
    CCC <- try(PopGenread(liste[xx],format),silent=TRUE) 
    if(is.na(CCC$matrix[1])){next}

    gen <- CCC$matrix
    pos <- CCC$positions
    ref <- CCC$reference 

# GFF stuff ----------------------------------------------------------------------
     
 gff.object.exists    <- FALSE
 FIT                  <- FALSE

 if(GFF.BOOL){

    gff.object.exists <- TRUE
    # if(length(grep("GFFRObjects",gffpath))!=0){
 
    if(SNP.DATA){
     if(!is.na(gff_liste[xx])){
       if(length(grep("GFFRObjects",gffpath))!=0){
        
          Robj                       <- load(gff_liste[xx])
          gff_object                 <- get(Robj[1])
       }else{
          gff_object                 <- gffRead(gff_liste[xx])        
       }

          gff_object_fit             <- fitting_gff_fast(pos,gff_object)              
          FIT                        <- TRUE
     }else{gff.object.exists <- FALSE} 
    }else{
      
        if(!is.na(gff_liste[xx])){
         gff_object                  <- gffRead(gff_liste[xx])
        }else{gff.object.exists <- FALSE}
    }   
   # Parsing GFF file
        if(gff.object.exists){

          gff_object                 <- parse_gff(gff_object,SNP.DATA=SNP.DATA)
          if(FIT){
          gff_object_fit             <- parse_gff(gff_object_fit,SNP.DATA=SNP.DATA)   
          GLOBAL.GFF$GFF             <- NULL                  
          }else{
          gff_object_fit             <- gff_object # in case of non SNP.DATA
          }

        }else{
          gff_object     <- FALSE
          gff_object_fit <- FALSE
        }

 }else{gff_object <- FALSE;gff_object_fit <- FALSE}

# GFF stuff END -------------------------------------------------

    if(is.matrix(gen)){
    
    nsites[xx]  <- dim(gen)[2] 
    nsites2[xx] <- dim(gen)[2] # important for SNPDATA
    
    #------------------------------------------------------#
    result    <- popgen(gen,Populations=populations,outgroup=outgroup,methods=methods,include.unknown=include.unknown,gff=gff_object_fit,FAST,SNP.DATA)
    #------------------------------------------------------#
    rm(gen) # delete gen
    
    }else{result <- NA;next}

# PROGRESS #######################################################
    if(progress_bar_switch){ # wegen parallized
    progr <- progressBar(xx,sizeliste, progr)
    }
   
    if(xx==ceiling(sizeliste/2)){gc()}
###################################################################
 
    if(is.list(result)){ # biallelic.sites exist
      
      # fill region.data
      populationsX[[xx]]        <- result$populations
      populations2[[xx]]        <- result$populations2
      popmissing[[xx]]          <- result$popmissing
      outgroupX[[xx]]           <- result$outgroup
      outgroup2[[xx]]           <- result$outgroup2

      datt                      <- result$data.sum
      
      CodingSNPS[[xx]]          <- datt$CodingSNPS
      UTRSNPS[[xx]]             <- datt$UTRSNPS
      IntronSNPS[[xx]]          <- datt$IntronSNPS
      ExonSNPS[[xx]]            <- datt$ExonSNPS
      GeneSNPS[[xx]]            <- datt$GeneSNPS	

      if(GFF.BOOL & !gff.object.exists){
       fillinit                  <- vector(,length(datt$biallelic.sites))
       CodingSNPS[[xx]]          <- fillinit
       UTRSNPS[[xx]]             <- fillinit
       IntronSNPS[[xx]]          <- fillinit
       ExonSNPS[[xx]]            <- fillinit
       GeneSNPS[[xx]]            <- fillinit
      }

      transitions[[xx]]         <- datt$transitions       # matrix_sv  transition war eine 1
      
     
      if(is.na(pos[1])){
       biallelic.sites[[xx]]     <- datt$biallelic.sites   # matrix_pos
       polyallelic.sites[[xx]]   <- datt$polyallelic.sites
       sites.with.gaps[[xx]]     <- datt$sites.with.gaps
       sites.with.unknowns[[xx]] <- datt$sites.with.unknowns
      }else{
       biallelic.sites2[[xx]]    <- datt$biallelic.sites
       biallelic.sites[[xx]]     <- pos[datt$biallelic.sites]
       nsites[xx]                <- biallelic.sites[[xx]][length(biallelic.sites[[xx]])] 
       polyallelic.sites[[xx]]   <- pos[datt$polyallelic.sites] # mhitbp
       sites.with.gaps[[xx]]     <- pos[datt$sites.with.gaps]
       sites.with.unknowns[[xx]] <- pos[datt$sites.with.unknowns]          
      }

     if(!big.data){

      biallelic.matrix[[xx]]           <- datt$biallelic.matrix  # matrix_pol
      colnames(biallelic.matrix[[xx]]) <- biallelic.sites[[xx]]
      
      if(GFF.BOOL & gff.object.exists){
       
       reading.frame[[xx]]              <- gff_object$reading.frame
       rev.strand[[xx]]                 <- gff_object$rev.strand
       Coding.matrix[[xx]]              <- gff_object$Coding
       UTR.matrix[[xx]]                 <- gff_object$UTR
       Intron.matrix[[xx]]              <- gff_object$Intron
       Exon.matrix[[xx]]                <- gff_object$Exon
       Gene.matrix[[xx]]                <- gff_object$Gene

       #
       #Coding.info[[xx]]                <- gff_object$Coding.info
       #UTR.info[[xx]]                   <- gff_object$UTR.info
       #Intron.info[[xx]]                <- gff_object$Intron.info
       #Exon.info[[xx]]                  <- gff_object$Exon.info
       #Gene.info[[xx]]                  <- gff_object$Gene.info
     
      }

     }else{ # BIG DATA FF

      biallelic.matrix[[xx]]           <- ff(datt$biallelic.matrix,dim=dim(datt$biallelic.matrix))  # matrix_pol
      close(biallelic.matrix[[xx]])
      dimnames(biallelic.matrix[[xx]]) <- list(rownames(datt$biallelic.matrix),biallelic.sites[[xx]])
      
     if(GFF.BOOL & gff.object.exists){
 
      if(dim(gff_object$Coding)[1]>0){
        
        fill                    <- as.matrix(gff_object$Coding)
        Coding.matrix[[xx]]     <- ff(fill,dim=dim(fill))
        close(Coding.matrix[[xx]])

        if(dim(gff_object_fit$Coding)[1]>0){
         reading.frame[[xx]]     <- gff_object_fit$reading.frame
         rev.strand[[xx]]        <- gff_object_fit$rev.strand
         fill                    <- as.matrix(gff_object_fit$Coding)
         Coding.matrix2[[xx]]    <- ff(fill,dim=dim(fill)) # wegen set.synnonsyn
         close(Coding.matrix2[[xx]])
        }
      #  Coding.info[[xx]]       <- ff(as.factor(gff_object$Coding.info))
      }
      if(dim(gff_object$UTR)[1]>0){
        fill                    <- as.matrix(gff_object$UTR)
        UTR.matrix[[xx]]        <- ff(fill,dim=dim(fill))
        close(UTR.matrix[[xx]])
      #  UTR.info[[xx]]          <- ff(as.factor(gff_object$UTR.info))
      }
      if(dim(gff_object$Intron)[1]>0){
        fill                    <- as.matrix(gff_object$Intron)   
        Intron.matrix[[xx]]     <- ff(fill,dim=dim(fill))
        close(Intron.matrix[[xx]])
      #  Intron.info[[xx]]       <- ff(as.factor(gff_object$Intron.info))
      }
      if(dim(gff_object$Exon)[1]>0){
        fill                    <- as.matrix(gff_object$Exon) 
        Exon.matrix[[xx]]       <- ff(fill,dim=dim(fill))
        close(Exon.matrix[[xx]])
      #  Exon.info[[xx]]         <- ff(as.factor(gff_object$Exon.info))
      } 
      if(dim(gff_object$Gene)[1]>0){
        fill                    <- as.matrix(gff_object$Gene)  
        Gene.matrix[[xx]]       <- ff(fill,dim=dim(fill))
        close(Gene.matrix[[xx]])
      #  Gene.info[[xx]]         <- ff(as.factor(gff_object$Gene.info))
      }    
     }
    }

      

      if(!is.na(ref[1])){
      reference[[xx]]           <- ref[datt$biallelic.sites]
      }

      matrix_codonpos[[xx]]     <- datt$matrix_codonpos   # codonpos. of biallelics
      synonymous[[xx]]          <- datt$synonymous        # synnonsyn
      matrix_freq[[xx]]         <- datt$matrix_freq
      n.singletons[[xx]]        <- datt$n.singletons      # unic
     
      n.nucleotides[[xx]]       <- datt$n.nucleotides     # sum_sam
      biallelic.compositions[[xx]]  <- datt$biallelic.compositions # TCGA
      biallelic.substitutions[[xx]] <- datt$biallelic.substitutions# subst
      minor.alleles[[xx]]       <- datt$minor.alleles # mutations
      codons[[xx]]              <- datt$codons
     
      
      # fill Genome data
      n.valid.sites[xx]       <- datt$n.valid.sites
      n.gaps[xx]              <- length(datt$sites.with.gaps)
      n.unknowns[xx]          <- length(datt$sites.with.unknowns)
      n.polyallelic.sites[xx] <- length(datt$polyallelic.sites)
      n.biallelic.sites[xx]   <- length(datt$biallelic.sites)
      trans.transv.ratio[xx]  <- datt$trans.transv.ratio

      Coding.region[xx]       <- datt$Coding_region_length
      UTR.region[xx]          <- datt$UTR_region_length
      Intron.region[xx]       <- datt$Intron_region_length
      Exon.region[xx]         <- datt$Exon_region_length
      Gene.region[xx]         <- datt$Gene_region_length    

    }# else{warnings("No biallelic position !")}

        
 }# End of For  
 
      # region.data
      region.data@populations      <- populationsX
      region.data@populations2     <- populations2
      region.data@popmissing       <- popmissing
      region.data@outgroup         <- outgroupX
      region.data@outgroup2        <- outgroup2

     
      region.data@CodingSNPS       <- CodingSNPS
      region.data@UTRSNPS          <- UTRSNPS
      region.data@IntronSNPS       <- IntronSNPS
      region.data@ExonSNPS         <- ExonSNPS
      region.data@GeneSNPS         <- GeneSNPS
      
      region.data@reading.frame    <- reading.frame
      region.data@rev.strand       <- rev.strand
      region.data@Coding.matrix    <- Coding.matrix
      region.data@Coding.matrix2   <- Coding.matrix2
      region.data@UTR.matrix       <- UTR.matrix
      region.data@Intron.matrix    <- Intron.matrix
      region.data@Exon.matrix      <- Exon.matrix
      region.data@Gene.matrix      <- Gene.matrix
      
      #region.data@Coding.info      <- Coding.info
      #region.data@UTR.info         <- UTR.info
      #region.data@Intron.info      <- Intron.info
      #region.data@Exon.info        <- Exon.info
      #region.data@Gene.info        <- Gene.info
      
      region.data@transitions      <- transitions  # matrix_sv  transition war eine 1
      region.data@biallelic.matrix <- biallelic.matrix# matrix_pol
      region.data@biallelic.sites  <- biallelic.sites # matrix_pos
      region.data@biallelic.sites2 <- biallelic.sites2
      region.data@reference        <- reference
      region.data@matrix_codonpos  <- matrix_codonpos # codonpos. of biallelics
      region.data@synonymous       <- synonymous# synnonsyn
      region.data@matrix_freq      <- matrix_freq
      region.data@n.singletons     <- n.singletons # unic
      region.data@polyallelic.sites <- polyallelic.sites # mhitbp
      region.data@n.nucleotides    <- n.nucleotides # sum_sam
      region.data@biallelic.compositions  <- biallelic.compositions  # TCGA
      region.data@biallelic.substitutions <- biallelic.substitutions # subst
      region.data@minor.alleles    <- minor.alleles # mutations
      region.data@codons           <- codons
      region.data@sites.with.gaps  <- sites.with.gaps # gaps
      region.data@sites.with.unknowns <- sites.with.unknowns

      region.stats@nucleotide.diversity   <- init
      region.stats@haplotype.diversity    <- init
      region.stats@haplotype.counts       <- init       # sfreqh
      region.stats@minor.allele.freqs     <- init       # JFD
      region.stats@biallelic.structure    <- init       # SXX
      region.stats@linkage.disequilibrium <- init

      # GENOME data
      genome@big.data            <- big.data
      names(nsites)              <- liste2
      genome@n.sites             <- nsites
      genome@n.sites2            <- nsites2
      genome@n.valid.sites       <- n.valid.sites    
      genome@n.gaps              <- n.gaps     
      genome@n.unknowns          <- n.unknowns   
      genome@n.polyallelic.sites <- n.polyallelic.sites
      genome@n.biallelic.sites   <- n.biallelic.sites
      genome@trans.transv.ratio  <- trans.transv.ratio

      genome@Coding.region     <- Coding.region
      genome@UTR.region        <- UTR.region
      genome@Intron.region     <- Intron.region
      genome@Exon.region       <- Exon.region
      genome@Gene.region       <- Gene.region


genome@region.data  <- region.data
genome@region.stats <- region.stats


genome@Pop_Neutrality$calculated <- FALSE
genome@Pop_FSTN$calculated       <- FALSE
genome@Pop_FSTH$calculated       <- FALSE
genome@Pop_MK$calculated         <- FALSE
genome@Pop_Linkage$calculated    <- FALSE
genome@Pop_Recomb$calculated     <- FALSE
genome@Pop_Slide$calculated      <- FALSE
genome@Pop_Detail$calculated     <- FALSE

if(FAST){genome@Pop_Slide$calculated <- TRUE}
if(GFF.BOOL){genome@gff.info<-TRUE}else{genome@gff.info <- FALSE}

genome@snp.data <- SNP.DATA

cat("\n")
return(genome)

}# End of Function

