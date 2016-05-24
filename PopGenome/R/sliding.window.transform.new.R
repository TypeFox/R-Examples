setGeneric("sliding.window.transform.new", function(object,width=7,jump=5,type=1,start.pos=FALSE,end.pos=FALSE) standardGeneric("sliding.window.transform.new"))
 setMethod("sliding.window.transform.new", "GENOME",

 function(object,width,jump,type,start.pos,end.pos){

n.region.names  <- length(object@region.names)

## PROGRESS #########################
 progr <- progressBar()
#####################################


if(start.pos[1]!=FALSE){
  object     <- splitting.data(object, positions=list(start.pos:end.pos), type=type)
  cat("\n")
}else{
  start.pos  <- 0
}



genomeobj               <-  new("GENOME") 
ddatt                   <-  new("region.data")
NAM                     <-  NULL
count                   <-  1

init.stuff       <- check_init_length(object,width,jump,type)
init.length      <- sum(init.stuff,na.rm=TRUE) + sum(is.na(init.stuff))
init             <- vector("list",init.length)

biallelic.matrix <- init
outgroup         <- init
popmissing       <- init 
populations      <- init
transitions      <- init 
synonymous       <- init
biallelic.sites  <- init

n.sites          <- numeric(init.length)



for(xx in 1:n.region.names){

  if(is.na(init.stuff[xx])){
    # biallelic.matrix[[count]]       <- NULL
    NAM                             <- c(NAM,object@region.names[xx])
    count                           <- count + 1
    next
  }


 if(object@big.data){
  open(object@region.data@biallelic.matrix[[xx]])
  bial <- object@region.data@biallelic.matrix[[xx]] 
 }else{
  bial <- popGetBial(object,xx)
 }

    bial.sites            <- as.numeric(colnames(bial))                 
    repeatlength          <- init.stuff[xx]  
   
   MERKEN <- 1
   for(zz in 1:repeatlength){
 
        
        start      <- ((zz-1) * jump + 1)
        end 	   <- ((zz-1) * jump + width) 
       

       if(type==1){
        
       
	 window	   <- start:end 
      
         n.sites[count] <- object@region.data@biallelic.sites[[xx]][end] -  object@region.data@biallelic.sites[[xx]][start]   +  1     

         if(object@big.data){
               biallelic.matrix[[count]]  <- ff(bial[,window,drop=FALSE],dim=dim(bial[,window,drop=FALSE]))  
               close(biallelic.matrix[[count]])             
         }else{
               biallelic.matrix[[count]]  <- bial[,window,drop=FALSE]
         }

        outgroup[[count]]         <- object@region.data@outgroup[[xx]]
        populations[[count]]      <- object@region.data@populations[[xx]]
        popmissing[[count]]       <- object@region.data@popmissing[[xx]]
        synonymous[[count]]       <- object@region.data@synonymous[[xx]][window]
        transitions[[count]]      <- object@region.data@transitions[[xx]][window]
        biallelic.sites[[count]]  <- object@region.data@biallelic.sites[[xx]][window]
        
        count                     <- count + 1     
        NAM                       <- c(NAM,object@region.names[xx])

       }
        
       
        if(type==2){
                   
                 bialpos            <- .Call("find_windowC",bial.sites,(start + start.pos),(end + start.pos),MERKEN)
                 
                 n.sites[count]     <- (end + start.pos) - (start + start.pos) + 1
                 
                #if(length(bialpos)>0){
                #   bialpos        <- bialpos[1]:bialpos[2] 
                #   MERKEN         <- bialpos[1] 
                #}else{next}
              
            # bialpos        <- as.numeric(colnames(bial))
            #ids            <- (bial.sites >= start) & (bial.sites<=end) 
            # ids            <- is.element(bialpos,window)
            #bialpos        <- which(ids)

           if(length(bialpos)!=0){

               bialpos        <- bialpos[1]:bialpos[2] 
               MERKEN         <- bialpos[1]
                                    


              if(object@big.data){
               biallelic.matrix[[count]]  <- ff(bial[,bialpos,drop=FALSE],dim=dim(bial[,bialpos,drop=FALSE]))  
               #close(biallelic.matrix[[count]]) # FIXME             
              }else{
               biallelic.matrix[[count]]  <- bial[,bialpos,drop=FALSE]
              }

              outgroup[[count]]          <- object@region.data@outgroup[[xx]]
              populations[[count]]       <- object@region.data@populations[[xx]]
              popmissing[[count]]        <- object@region.data@popmissing[[xx]]
              synonymous[[count]]        <- object@region.data@synonymous[[xx]] [bialpos]
              transitions[[count]]       <- object@region.data@transitions[[xx]][bialpos]
              biallelic.sites[[count]]   <- object@region.data@biallelic.sites[[xx]][bialpos]
              count                      <- count + 1     
              NAM                        <- c(NAM,object@region.names[xx])

           }else{

             # biallelic.matrix[[count]]       <- NULL
             NAM                             <- c(NAM,object@region.names[xx])
             count                           <- count + 1

           }
       } 
 
# slidenames[zz]  <- paste(start,"-",end,sep="")
 

} # end for Sliding

 
  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################
 if(object@big.data){
    close(bial)
 }

}# End over all regions

genomeobj@populations               <- object@populations
genomeobj@region.names              <- paste(1:length(NAM),":",NAM,sep="")
genomeobj@genelength                <- length(NAM)
genomeobj@n.sites                   <- n.sites

 if(length(ddatt@popmissing)==0){ddatt@popmissing <- vector("list",length(genomeobj@region.names))}
        # weil liste fuellen mit NULL nicht funktioniert, eh besser die Sachen schon vorher zu definieren


ddatt@biallelic.matrix <- biallelic.matrix
ddatt@outgroup         <- outgroup
ddatt@popmissing       <- popmissing
ddatt@populations      <- populations
ddatt@transitions      <- transitions
ddatt@synonymous       <- synonymous
ddatt@biallelic.sites  <- biallelic.sites



genomeobj@region.data               <- ddatt
genomeobj@Pop_Neutrality$calculated <- FALSE
genomeobj@Pop_FSTN$calculated       <- FALSE
genomeobj@Pop_FSTH$calculated       <- FALSE
genomeobj@Pop_MK$calculated         <- FALSE
genomeobj@Pop_Recomb$calculated     <- FALSE
genomeobj@Pop_Linkage$calculated    <- FALSE
genomeobj@Pop_Slide$calculated      <- TRUE
genomeobj@Pop_Detail$calculated     <- FALSE
genomeobj@big.data                  <- object@big.data
genomeobj@snp.data                  <- object@snp.data

return(genomeobj)

 })
 
 
