##################################################################################
## Sliding Window Transform FAST VERSION for scanning SNP chunks seperately
##################################################################################

setGeneric("slide.snp.sep", function(object,width=7,jump=5,type=1) standardGeneric("slide.snp.sep"))
 setMethod("slide.snp.sep", "GENOME",

 function(object,width,jump,type){
  
# Informations we need from the old object
N.Bial.Sites <- object@n.biallelic.sites
N.SITES      <- object@n.sites
CHUNK.NAMES  <- object@region.names


cat("Prepare ..." ,"\n")
object       <- concatenate_to_whole_genome(object,length(object@region.names))



  genomeobj               <-  new("GENOME") 
  ddatt                   <-  new("region.data")
  XXX                     <-  object@region.data 


 if(object@big.data){
  #bial <- object@region.data@biallelic.matrix[[1]]
   open(object@region.data@biallelic.matrix[[1]])
   genomeobj@BIG.BIAL[[1]] <- object@region.data@biallelic.matrix[[1]] 
   bial                    <- object@region.data@biallelic.matrix[[1]]
   #close(object@region.data@biallelic.matrix[[1]])
 }else{
   bial <- popGetBial(object,1)
 }

# Init ---------------------
  init              <- NULL
  CSLIDE.POS        <- init            
  Cbiallelic.matrix <- init
  Cbiallelic.sites  <- init
  Coutgroup         <- init
  Cpopulations      <- init
  Cpopmissing       <- init
  Ctransitions      <- init
  Csynonymous       <- init
  CNAMES            <- init
# ---------------------------

repeatlength <- numeric(length(N.Bial.Sites))

# Check the blocks of the baillelic.matrix
bial.pos.ids <- matrix(,length(N.Bial.Sites),2)
start        <- 1
ADDX         <- 0

for(xx in 1:length(N.Bial.Sites)){
   ADDX              <- ADDX + N.Bial.Sites[xx]
   bial.pos.ids[xx,] <- c(start,ADDX)
   start             <- ADDX + 1
}
###########################################

# we have to add values to the SLIDE Pos, because we use the BIGMATRIX, which includes all pos
ADD <- 0


for (xx in 1: length(N.Bial.Sites)){

  bial.sites               <- object@region.data@biallelic.sites[[1]][bial.pos.ids[xx,1]:bial.pos.ids[xx,2]]      
  # Calculate type 3 = reference Biallelic Matrix. type 5 = reference GEN   
  if(type==1){repeatlength[xx] <- ceiling( (N.Bial.Sites[xx]-width+1)/jump)}
  if(type==2){repeatlength[xx] <- ceiling( (N.SITES[xx]-width+1)/jump)} 
 
     
  # init -------------------------------------------------------
  init             <- vector("list",repeatlength[xx])
  SLIDE.POS        <- init            
  biallelic.matrix <- init
  biallelic.sites  <- init
  outgroup         <- init
  populations      <- init
  popmissing       <- init
  transitions      <- init
  synonymous       <- init
  NAMES            <- vector(,repeatlength[xx])
  #-------------------------------------------------------------

cat(CHUNK.NAMES[xx],"\n")
## PROGRESS #########################
progr <- progressBar()
#####################################
   
   MERKEN <- 1 
   for(zz in 1:repeatlength[xx]){
 
        
        start      <- ((zz-1) * jump + 1)
        end 	   <- ((zz-1) * jump + width) 
        NAMES[zz]  <- paste(CHUNK.NAMES[xx],start,"-",end,":")
        #  NAMES[zz]  <- CHUNK.NAMES[xx]

       if(type==1){
     
         window	   <- start:end   
                 
    
        #genomeobj@DATA[[count]] <- new("DATA")
        #ddatt@biallelic.matrix[[count]] <- object@region.data@biallelic.matrix[[xx]][,window,drop=FALSE]

         if(object@big.data){

                SLIDE.POS[[zz]]         <- window + ADD
               
#               biallelic.matrix[[zz]]  <- ff(bial[,window,drop=FALSE],dim=dim(bial#[,window,drop=FALSE]))
               
         }else{
               biallelic.matrix[[zz]]  <- bial[,window,drop=FALSE]
         }

        outgroup[[zz]]         <- XXX@outgroup[[1]]
        populations[[zz]]      <- XXX@populations[[1]]
        # popmissing[[zz]]       <- XXX@popmissing[[1]]
        synonymous[[zz]]       <- XXX@synonymous[[1]][window+ADD]
        transitions[[zz]]      <- XXX@transitions[[1]][window+ADD]
        biallelic.sites[[zz]]  <- XXX@biallelic.sites[[1]][window+ADD]
        
       }
        
             
        if(type==2){
      
            
             # ids            <- (bial.sites >= start) & (bial.sites<=end) 
             # bialpos        <- which(ids)
             # if(length(bialpos)==0){next}

 
             # FASTER ! but check if runs correctly again  !
                bialpos            <- .Call("find_windowC",bial.sites,start,end,MERKEN)
                if(length(bialpos)>0){
                   bialpos        <- bialpos[1]:bialpos[2] 
                   MERKEN         <- bialpos[1] 
                }else{next}
                  
              
               
              

              if(object@big.data){
               SLIDE.POS[[zz]]         <- bialpos  + ADD
               
      
       #biallelic.matrix[[zz]]  <- ff(bial[,bialpos,drop=FALSE],dim=dim(bial[,bialpos,drop=FALSE]))               
              }else{
               biallelic.matrix[[zz]]  <- bial[,bialpos,drop=FALSE]
              }

              outgroup[[zz]]          <- XXX@outgroup[[1]]
              populations[[zz]]       <- XXX@populations[[1]]
            #  popmissing[[zz]]        <- XXX@popmissing[[1]]
              synonymous[[zz]]        <- XXX@synonymous[[1]][bialpos + ADD]
              transitions[[zz]]       <- XXX@transitions[[1]][bialpos + ADD]
              biallelic.sites[[zz]]   <- XXX@biallelic.sites[[1]][bialpos + ADD]
           
         }
  
 # PROGRESS #######################################################
    progr <- progressBar(zz,repeatlength[xx], progr)
 ################################################################### 

} # end for Sliding

 
ADD               <- ADD    + N.Bial.Sites[xx] 

CSLIDE.POS        <- c(CSLIDE.POS,SLIDE.POS)          
Cbiallelic.matrix <- c(Cbiallelic.matrix,biallelic.matrix)
Cbiallelic.sites  <- c(Cbiallelic.sites,biallelic.sites)
Coutgroup         <- c(Coutgroup,outgroup)
Cpopulations      <- c(Cpopulations,populations)
Cpopmissing       <- c(Cpopmissing,popmissing)
Ctransitions      <- c(Ctransitions,transitions)
Csynonymous       <- c(Csynonymous,synonymous)
CNAMES            <- c(CNAMES,NAMES)
 

} # end over all chunks/chr

ddatt@biallelic.matrix <- Cbiallelic.matrix
ddatt@biallelic.sites  <- Cbiallelic.sites
ddatt@populations      <- Cpopulations
ddatt@outgroup         <- Coutgroup
ddatt@popmissing       <- Cpopmissing
ddatt@transitions      <- Ctransitions
ddatt@synonymous       <- Csynonymous


genomeobj@SLIDE.POS                 <- CSLIDE.POS
genomeobj@populations               <- object@populations
genomeobj@region.names              <- CNAMES
genomeobj@genelength                <- length(CNAMES) # repeatlength
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

return(genomeobj)

 })
 
 
