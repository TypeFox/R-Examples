setGeneric("splitting.data", function(object,subsites=FALSE, positions=FALSE, type=1, whole.data=TRUE) standardGeneric("splitting.data"))
 setMethod("splitting.data", "GENOME",

 function(object,subsites,positions,type,whole.data){


# if(subsites[1]!=FALSE & !object@gff.info){
# stop("No GFF/GTF informations !")
# }

## SNP DATA SCAN DATA SEPERATELY !
if(whole.data==FALSE){
   object <- splitting.data.sep(object,subsites,positions,type)
return(object)
}


if( object@genelength > 1){ 

 object <- concatenate_to_whole_genome(object,object@genelength)

}


genomeobj               <-  new("GENOME") 
ddatt                   <-  new("region.data")
XXX                     <-  object@region.data 


genomeobj@BIG.BIAL[[1]] <- object@region.data@biallelic.matrix[[1]]
 
if(!is.list(positions)){

 if(subsites=="coding"){
    Gene.matrix <- object@region.data@Coding.matrix[[1]][,]
    Gene.matrix <- unique(Gene.matrix)
    LENGTH      <- GLOBAL.GFF$Coding
 }
 if(subsites=="exon"){
    Gene.matrix <- object@region.data@Exon.matrix[[1]][,]
    Gene.matrix <- unique(Gene.matrix)	
    LENGTH      <- GLOBAL.GFF$Exon
 }
 if(subsites=="utr"){
    Gene.matrix <- object@region.data@UTR.matrix[[1]][,]
    Gene.matrix <- unique(Gene.matrix)
    LENGTH      <- GLOBAL.GFF$UTR
 }
 if(subsites=="intron"){
    Gene.matrix <- object@region.data@Intron.matrix[[1]][,]
    Gene.matrix <- unique(Gene.matrix)  
    LENGTH      <- GLOBAL.GFF$Intron
 }
 if(subsites=="gene"){
    Gene.matrix <- object@region.data@Gene.matrix[[1]][,]
    Gene.matrix <- unique(Gene.matrix)
    LENGTH      <- GLOBAL.GFF$Gene
 }
 
 if(length(dim(Gene.matrix)[1])==0 || dim(Gene.matrix)[1]==0){
  cat("\n")
  stop(paste("No ",subsites, "defined !"))
 }
 
}
 
 # init -------------------------------------------------------
  if(!is.list(positions)){
     init             <- vector("list",dim(Gene.matrix)[1])
     region.names     <- character(dim(Gene.matrix)[1])
     n.sites          <- numeric(dim(Gene.matrix)[1])
  }else{
     init             <- vector("list",length(positions))    
     region.names     <- character(length(positions))  
     n.sites          <- numeric(length(positions))
  }

  SLIDE.POS        <- init            
  biallelic.matrix <- init
  biallelic.sites  <- init
  outgroup         <- init
  populations      <- init
  popmissing       <- init
  transitions      <- init
  synonymous       <- init
  
  biallelic.substitutions <- init
  
   SPLITGENE        <- FALSE

  if(subsites=="gene"){

   SPLITGENE        <- TRUE
   CodingSNPS       <- init
   ExonSNPS         <- init
   IntronSNPS       <- init
 
  }
  #-------------------------------------------------------------

    if(!is.list(positions) & !object@snp.data){
    
     ### Progress
     progr <- progressBar()
     ###
  
     for(zz in 1:dim(Gene.matrix)[1]){
   
        bialpos                       <- Gene.matrix[zz,1]:Gene.matrix[zz,2]
        bialpos                       <- is.element(XXX@biallelic.sites[[1]],bialpos) # wird hier der Vek neu erzeugt ?
        #bialpos                      <- is.element(object@region.data@biallelic.sites[[1]],bialpos)
        bialpos                       <- which(bialpos)
   

        region.names[zz]              <- paste(Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
        n.sites[zz]                   <- (Gene.matrix[zz,2]-Gene.matrix[zz,1]) + 1
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
       # popmissing[[zz]]             <- XXX@popmissing[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos]  


     ## Progress
     progr <- progressBar(zz,dim(Gene.matrix)[1], progr)
     ####
       
     }
    }
 
    if(!is.list(positions) & object@snp.data){

     ### Progress
     progr <- progressBar()
     ###

      

     for(zz in 1:dim(Gene.matrix)[1]){
   
        
         #bialpos                       <- Gene.matrix[zz,1]:Gene.matrix[zz,2]
         #bialpos                       <- match(bialpos,XXX@biallelic.sites2[[1]])
         #bialpos                       <- bialpos[!is.na(bialpos)]
           
         bialpos                       <- .Call("find_windowC",XXX@biallelic.sites[[1]],Gene.matrix[zz,1],Gene.matrix[zz,2],1)
         # bialpos                       <- XXX@biallelic.sites2[[1]][bialpos[1]:bialpos[2]] 
         if(length(bialpos)!=0){ 
         bialpos                       <- bialpos[1]:bialpos[2]     
         }else{
          region.names[zz]             <- paste(Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
          n.sites[zz]                  <- Gene.matrix[zz,2] - Gene.matrix[zz,1] + 1 
          next
         }

         #start.region                  <- XXX@biallelic.sites[[1]][bialpos[1]]
         #end.region                    <- XXX@biallelic.sites[[1]][bialpos[length(bialpos)]]
         
        region.names[zz]              <- paste(Gene.matrix[zz,1],"-",Gene.matrix[zz,2])
        #region.names[zz]              <- paste(start.region,"-",end.region)
        #n.sites[zz]                   <- LENGTH[zz]
        n.sites[zz]                   <- Gene.matrix[zz,2] - Gene.matrix[zz,1] + 1 
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
       # popmissing[[zz]]              <- XXX@popmissing[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos]  

   
        if(SPLITGENE){

        CodingSNPS[[zz]]                    <- XXX@CodingSNPS[[1]][bialpos]        
        ExonSNPS[[zz]]                      <- XXX@ExonSNPS[[1]][bialpos]
        IntronSNPS[[zz]]                    <- XXX@IntronSNPS[[1]][bialpos]

        }
 
        
     ## Progress
     progr <- progressBar(zz,dim(Gene.matrix)[1], progr)
     ####
       
     }
    }
   


   if(is.list(positions) & type==1) {

     ### Progress
     progr <- progressBar()
     ###
    
    if(object@gff.info){
     SPLITGENE        <- TRUE
     CodingSNPS       <- init
     ExonSNPS         <- init
     IntronSNPS       <- init  
    }
    

      for(zz in 1:length(positions)){

        bialpos                       <- positions[[zz]]
        region.names[zz]              <- paste(positions[[zz]][1],"-",positions[[zz]][length(positions[[zz]])])
        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
       # popmissing[[zz]]              <- XXX@popmissing[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos] 
        n.sites[zz]                   <- biallelic.sites[[zz]][length(biallelic.sites[[zz]])] - biallelic.sites[[zz]][1] + 1
        biallelic.substitutions[[zz]] <- XXX@biallelic.substitutions[[1]][, bialpos, drop=FALSE]
 
     # GFF stuff
      if(SPLITGENE){
        CodingSNPS[[zz]]                    <- XXX@CodingSNPS[[1]][bialpos]        
        ExonSNPS[[zz]]                      <- XXX@ExonSNPS[[1]][bialpos]
        IntronSNPS[[zz]]                    <- XXX@IntronSNPS[[1]][bialpos]
      }       

     ## Progress
     progr <- progressBar(zz,length(positions), progr)
     ####

      }
   }

   if(is.list(positions) & type==2) {

     ### Progress
     progr <- progressBar()
     ###
   
      if(object@gff.info){
       SPLITGENE        <- TRUE
       CodingSNPS       <- init
       ExonSNPS         <- init
       IntronSNPS       <- init  
      }
    

     gc()
     for(zz in 1:length(positions)){
            
        bialpos                       <- is.element(XXX@biallelic.sites[[1]],positions[[zz]])
        bialpos                       <- which(bialpos)
               
        region.names[zz]              <- paste(positions[[zz]][1],"-",positions[[zz]][length(positions[[zz]])])
	n.sites[zz]                   <- length(positions[[zz]])

        ## Progress
        progr <- progressBar(zz,length(positions), progr)
        ####
          
        if(length(bialpos)==0){next}

        SLIDE.POS[[zz]]               <- bialpos
        outgroup[[zz]]                <- XXX@outgroup[[1]]
        populations[[zz]]             <- XXX@populations[[1]]
       # popmissing[[zz]]              <- XXX@popmissing[[1]]
        synonymous[[zz]]              <- XXX@synonymous[[1]][bialpos]
        transitions[[zz]]             <- XXX@transitions[[1]][bialpos]
        biallelic.sites[[zz]]         <- XXX@biallelic.sites[[1]][bialpos] 
        

        #positions[[zz]][length(positions[[zz]])] - positions[[zz]][1] + 1 

        biallelic.substitutions[[zz]] <- XXX@biallelic.substitutions[[1]][, bialpos, drop=FALSE]

        # GFF stuff
	if(SPLITGENE){
        CodingSNPS[[zz]]                    <- XXX@CodingSNPS[[1]][bialpos]        
        ExonSNPS[[zz]]                      <- XXX@ExonSNPS[[1]][bialpos]
        IntronSNPS[[zz]]                    <- XXX@IntronSNPS[[1]][bialpos]
        }                 
     }
   } 


# BUG vector("list",1)[[1]] <- NULL  => list()
# if(length(popmissing)==0){
##  popmissing <- vector("list",1)
#}
# also change in splitting.data.sep !
#


ddatt@UTR.matrix       <- object@region.data@UTR.matrix
ddatt@Exon.matrix      <- object@region.data@Exon.matrix
ddatt@Gene.matrix      <- object@region.data@Gene.matrix
ddatt@Intron.matrix    <- object@region.data@Intron.matrix
ddatt@Coding.matrix    <- object@region.data@Coding.matrix
ddatt@biallelic.matrix <- biallelic.matrix
ddatt@biallelic.sites  <- biallelic.sites
ddatt@populations      <- populations
ddatt@outgroup         <- outgroup
ddatt@popmissing       <- popmissing
ddatt@transitions      <- transitions
ddatt@synonymous       <- synonymous
ddatt@biallelic.substitutions <- biallelic.substitutions


 if(SPLITGENE){

  ddatt@CodingSNPS     <- CodingSNPS
  ddatt@ExonSNPS       <- ExonSNPS
  ddatt@IntronSNPS     <- IntronSNPS
 
 }
 

genomeobj@SLIDE.POS                 <- SLIDE.POS
genomeobj@populations               <- object@populations
genomeobj@region.names              <- region.names
genomeobj@n.sites                   <- n.sites
genomeobj@genelength                <- length(SLIDE.POS)
genomeobj@region.data               <- ddatt
genomeobj@Pop_Neutrality$calculated <- FALSE
genomeobj@Pop_FSTN$calculated       <- FALSE
genomeobj@Pop_FSTH$calculated       <- FALSE
genomeobj@Pop_MK$calculated         <- FALSE
genomeobj@Pop_Linkage$calculated    <- FALSE
genomeobj@Pop_Recomb$calculated     <- FALSE
genomeobj@Pop_Slide$calculated      <- TRUE
genomeobj@Pop_Detail$calculated     <- FALSE
genomeobj@big.data                  <- object@big.data
genomeobj@snp.data                  <- object@snp.data

rm(object)
gc()

return(genomeobj)

})
  
## SUBFUNCTIONS
my_chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
#------------------------------------------------------------

