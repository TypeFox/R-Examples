##################################################################################
## Sliding Window Transform FAST VERSION (sliding.window.transform.fast)
##################################################################################
## slightly modified for Jack.knife mechanisms
# the only change is genomeobj@jack.knife <- TRUE
# see last line 

setGeneric("jack.knife.transform", function(object,width=7,jump=5,type=1,start.pos=FALSE, end.pos=FALSE) standardGeneric("jack.knife.transform"))
 setMethod("jack.knife.transform", "GENOME",

 function(object,width,jump,type,start.pos,end.pos){
  

## PROGRESS #########################
 progr <- progressBar()
#####################################

# When you want to scan a specific region

if(start.pos[1]!=FALSE){

  # Suche beg. SNP of original data
  if(type==2){
  snp.begin  <- .Call("whichbigger_C",start.pos,object@region.data@biallelic.sites[[1]])
  snp.begin  <- snp.begin - 1
  }else{#type=1
  snp.begin  <- start.pos - 1 
  }
 
  start.pos <- start.pos - 1
 
  # split the data 
  object     <- splitting.data(object, positions=list(start.pos:end.pos), type=type)
  cat("\n")
   

}else{

  if(length(object@keep.start.pos)!=0){ # wenn bei readVCF z.b nur ne bestimmte region !
  start.pos  <- object@keep.start.pos
  }else{
  start.pos  <- 0
  }

snp.begin <- 0

}


genomeobj               <-  new("GENOME") 
ddatt                   <-  new("region.data")
XXX                     <-  object@region.data 


 if(object@big.data){
  # bial <- object@region.data@biallelic.matrix[[1]]

   if(is(object@region.data@biallelic.matrix[[1]])[1]=="ff_matrix"){
      open(object@region.data@biallelic.matrix[[1]]) # FIXME
   }
 
     if(length(object@BIG.BIAL)==0){
     genomeobj@BIG.BIAL[[1]]   <- object@region.data@biallelic.matrix[[1]]
     }else{
     genomeobj@BIG.BIAL[[1]]   <- object@BIG.BIAL[[1]] # when splitting the data before
     }     

     colsbial                  <- length(object@region.data@biallelic.sites[[1]])
    
   
   #close(object@region.data@biallelic.matrix[[1]])

 }else{

   colsbial                  <-  length(object@region.data@biallelic.sites[[1]])
  
    if(length(object@BIG.BIAL)==0){
     bial   <- object@region.data@biallelic.matrix[[1]]
     }else{
     bial   <- object@BIG.BIAL[[1]] # when splitting the data before
    }    # when big.data is FALSE --> copy biallelic matrix to each region !    
 }

 bial.sites        <- object@region.data@biallelic.sites[[1]] 

                    
  # Calculate type 3 = reference Biallelic Matrix. type 5 = reference GEN   
  if(type==1){repeatlength <- ceiling( (colsbial-width+1)/jump)}
  if(type==2){repeatlength <- ceiling( (object@n.sites[1]-width+1)/jump)} 
      
  # init -------------------------------------------------------
  init             <- vector("list",repeatlength)
  SLIDE.POS        <- init            
  biallelic.matrix <- init
  biallelic.sites  <- init
  outgroup         <- init
  populations      <- init
  popmissing       <- init
  transitions      <- init
  synonymous       <- init
  NAMES            <- vector(,repeatlength)
  n.sites          <- numeric(repeatlength)
  #-------------------------------------------------------------

   MERKEN <- 1 
   for(zz in 1:repeatlength){
 
        
        start      <- ((zz-1) * jump + 1)
        end 	   <- ((zz-1) * jump + width) 
        
        

       if(type==1){
     
         NAMES[zz]   <- paste(start + snp.begin ,"-",end + snp.begin ,":")
         window	     <- start:end
         n.sites[zz] <- XXX@biallelic.sites[[1]][end] -  XXX@biallelic.sites[[1]][start] + 1     
    
        #genomeobj@DATA[[count]] <- new("DATA")
        #ddatt@biallelic.matrix[[count]] <- object@region.data@biallelic.matrix[[xx]][,window,drop=FALSE]

         if(object@big.data){

                SLIDE.POS[[zz]]         <- window + snp.begin
#               biallelic.matrix[[zz]]  <- ff(bial[,window,drop=FALSE],dim=dim(bial#[,window,drop=FALSE]))
               
         }else{

               biallelic.matrix[[zz]]  <- bial[,window + snp.begin,drop=FALSE]

         }

        outgroup[[zz]]         <- XXX@outgroup[[1]]
        populations[[zz]]      <- XXX@populations[[1]]
        #popmissing[[zz]]       <- XXX@popmissing[[1]]
        synonymous[[zz]]       <- XXX@synonymous[[1]][window]
        transitions[[zz]]      <- XXX@transitions[[1]][window]
        biallelic.sites[[zz]]  <- XXX@biallelic.sites[[1]][window]
        
       }
        
             
        if(type==2){
      	    	    
             n.sites[zz]          <- (end + start.pos) - (start + start.pos) + 1       

             NAMES[zz]            <- paste(( start + start.pos ),"-", (end + start.pos) ,":")

             # ids            <- (bial.sites >= start) & (bial.sites<=end) 
             # bialpos        <- which(ids)
             # if(length(bialpos)==0){next}

            
 
             # FASTER ! but check if runs correctly again  !
                bialpos            <- .Call("find_windowC",bial.sites,(start+start.pos),(end+start.pos),MERKEN)
                if(length(bialpos)>0){
                   bialpos        <- bialpos[1]:bialpos[2] 
                   MERKEN         <- bialpos[1] 
                }else{next}
                  
              
               
              

              if(object@big.data){
              SLIDE.POS[[zz]]         <- bialpos + snp.begin       
       #biallelic.matrix[[zz]]  <- ff(bial[,bialpos,drop=FALSE],dim=dim(bial[,bialpos,drop=FALSE]))               
              }else{
               biallelic.matrix[[zz]]  <- bial[,bialpos + snp.begin,drop=FALSE]
              }

              outgroup[[zz]]          <- XXX@outgroup[[1]]
              populations[[zz]]       <- XXX@populations[[1]]
              #popmissing[[zz]]        <- XXX@popmissing[[1]]
              synonymous[[zz]]        <- XXX@synonymous[[1]][bialpos]
              transitions[[zz]]       <- XXX@transitions[[1]][bialpos]
              biallelic.sites[[zz]]   <- XXX@biallelic.sites[[1]][bialpos]
           
         }
  
  # PROGRESS #######################################################
    progr <- progressBar(zz,repeatlength, progr)
  ###################################################################

} # end for Sliding


 


ddatt@biallelic.matrix <- biallelic.matrix
ddatt@biallelic.sites  <- biallelic.sites
ddatt@populations      <- populations
ddatt@outgroup         <- outgroup
ddatt@popmissing       <- popmissing
ddatt@transitions      <- transitions
ddatt@synonymous       <- synonymous


genomeobj@SLIDE.POS                 <- SLIDE.POS
genomeobj@populations               <- object@populations
genomeobj@region.names              <- NAMES
genomeobj@n.sites                   <- n.sites
genomeobj@genelength                <- length(NAMES)
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
genomeobj@jack.knife                <- TRUE

return(genomeobj)

 })
 
 
