save.session <- function(object,folder){

curr <- getwd()
dir.create(folder)
setwd(folder)

save(object,file="GENOME")

if(object@big.data){

  if(length(object@BIG.BIAL)>0){

    bial <- object@BIG.BIAL[[1]]
    ffsave(bial,file="bial_1")		
    #ffsave(object@BIG.BIAL,file="bial_1")

    #if(object@gff.info){
       
    #   if(is(object@region.data@Coding.matrix[[1]])[1]=="ff_matrix"){   
    #   ffsave(object@region.data@Coding.matrix,file="Coding.matrix")
    #   }
    #   if(is(object@region.data@Gene.matrix[[1]])[1]=="ff_matrix"){ 
    #   ffsave(object@region.data@Gene.matrix,file="Gene.matrix")
    #   }
    #   if(is(object@region.data@Exon.matrix[[1]])[1]=="ff_matrix"){ 
    #   ffsave(object@region.data@Exon.matrix,file="Exon.matrix")
    #   }
    #   if(is(object@region.data@Intron.matrix[[1]])[1]=="ff_matrix"){ 
    #   ffsave(object@region.data@Intron.matrix,file="Intron.matrix")
    #   }
    #   if(is(object@region.data@UTR.matrix[[1]])[1]=="ff_matrix"){ 
    #   ffsave(object@region.data@UTR.matrix,file="UTR.matrix")
    #   }

    #}

  }else{

    for(xx in 1:length(object@region.data@biallelic.matrix)){

      file.name <- paste("bial_",xx,sep="")
      bial      <- object@region.data@biallelic.matrix[[xx]]
      ffsave(bial,file=file.name)

      if(object@gff.info){
 
        if(is(object@region.data@Coding.matrix[[xx]])[1]=="ff_matrix"){ 
        file.name <- paste("Coding.matrix_",xx,sep="")
        Coding.matrix   <- object@region.data@Coding.matrix[[xx]]
        ffsave(Coding.matrix,file=file.name)
        }
        if(is(object@region.data@Gene.matrix[[xx]])[1]=="ff_matrix"){ 
        file.name <- paste("Gene.matrix_",xx,sep="")
        Gene.matrix       <- object@region.data@Gene.matrix[[xx]]
        ffsave(Gene.matrix,file=file.name)
        }
        if(is(object@region.data@Exon.matrix[[xx]])[1]=="ff_matrix"){ 
        file.name <- paste("Exon.matrix_",xx,sep="")
        Exon.matrix       <- object@region.data@Exon.matrix[[xx]]
        ffsave(Exon.matrix,file=file.name)
        }
        if(is(object@region.data@Intron.matrix[[xx]])[1]=="ff_matrix"){ 
        file.name <- paste("Intron.matrix_",xx,sep="")
        Intron.matrix     <- object@region.data@Intron.matrix[[xx]]
        ffsave(Intron.matrix,file=file.name)
        }
	if(is(object@region.data@UTR.matrix[[xx]])[1]=="ff_matrix"){ 
        file.name <- paste("UTR.matrix_",xx,sep="")
        UTR.matrix        <- object@region.data@UTR.matrix[[xx]]
        ffsave(UTR.matrix,file=file.name)
        }

      }
    }
  }  

} 

setwd(curr)
}

##########  ---- LOAD ------ #############

load.session <- function(folder){

curr       <- getwd()
setwd(folder)

name       <- load("GENOME")
object     <- get(name[1])
change     <- object@region.data

 if(object@big.data){

   if(length(object@BIG.BIAL)>0){

#   bial                 <- object@BIG.BIAL[[1]]
#   ff.object.name       <- ffload(bial,file="bial_1")	
#   ff.object.name       <- ffload(object@BIG.BIAL,file="bial_1")
#   ff.object            <- get(ff.object.name[1])
#   object@BIG.BIAL[[1]] <- ff.object
 
    ffload("bial_1",overwrite=TRUE)
    vcv                  <- get("bial")
    object@BIG.BIAL[[1]] <- vcv 	

    #if(object@gff.info){

    #  if(file.exists("Coding.matrix.RData")){
    # ff.object.name     <- ffload(object@region.data@Coding.matrix,file="Coding.matrix")
    #  ff.object          <- get(ff.object.name[1])
    #  change@Coding.matrix[[1]] <- ff.object
    #  }
    #  if(file.exists("Gene.matrix.RData")){
    #  ff.object.name     <- ffload(object@region.data@Gene.matrix,file="Gene.matrix")
    #  ff.object          <- get(ff.object.name[1])
    #  change@Gene.matrix[[1]] <- ff.object
    #  }
    #  if(file.exists("Exon.matrix.RData")){
    #  ff.object.name     <- ffload(object@region.data@Exon.matrix,file="Exon.matrix")
    #  ff.object          <- get(ff.object.name[1])
    #  change@Exon.matrix[[1]] <- ff.object
    #  }
    #  if(file.exists("Intron.matrix.RData")){
    #  ff.object.name     <- ffload(object@region.data@Intron.matrix,file="Intron.matrix")
    #  ff.object          <- get(ff.object.name[1])
    #  change@Intron.matrix[[1]] <- ff.object
    #  }
    #  if(file.exists("UTR.matrix.RData")){
    #  ff.object.name     <- ffload(object@region.data@UTR.matrix,file="UTR.matrix")
    #  ff.object          <- get(ff.object.name[1])
    #  change@UTR.matrix[[1]] <- ff.object
    #  }
    #}

  }else{

    for(xx in 1:length(object@region.data@biallelic.matrix)){
      
      file.name      <-  paste("bial_",xx,sep="")
      ffload(file.name)     
      ff.object      <-  get("bial")  
      change@biallelic.matrix[[xx]] <- ff.object
      open(change@biallelic.matrix[[xx]])
      
      if(object@gff.info){

       file.name    <- paste("Coding.matrix_",xx,sep="")
       if(file.exists(paste(file.name,".RData",sep=""))){
       ffload(file.name)
       ff.object    <- get("Coding.matrix")
       change@Coding.matrix[[xx]] <- ff.object
       open(change@Coding.matrix[[xx]]) 
       }
       file.name    <- paste("Gene.matrix_",xx,sep="")
       if(file.exists(paste(file.name,".RData",sep=""))){
       ffload(file.name)
       ff.object    <- get("Gene.matrix")
       change@Gene.matrix[[xx]] <- ff.object
       open(change@Gene.matrix[[xx]])
       }
       file.name    <- paste("Exon.matrix_",xx,sep="")
       if(file.exists(paste(file.name,".RData",sep=""))){
       ffload(file.name)
       ff.object    <- get("Exon.matrix")
       change@Exon.matrix[[xx]] <- ff.object
       open(change@Exon.matrix[[xx]])
       }
       file.name    <- paste("Intron.matrix_",xx,sep="")
       if(file.exists(paste(file.name,".RData",sep=""))){
       ffload(file.name)
       ff.object    <- get("Intron.matrix")
       change@Intron.matrix[[xx]] <- ff.object
       open(change@Intron.matrix[[xx]])
       }
       file.name    <- paste("UTR.matrix_",xx,sep="")
       if(file.exists(paste(file.name,".RData",sep=""))){
       ffload(file.name)
       ff.object    <- get("UTR.matrix")
       change@UTR.matrix[[xx]] <- ff.object
       open(change@UTR.matrix[[xx]]) 
       }
      }
    }
  }  

 }

object@region.data <- change
setwd(curr)
assign(folder,object,envir=parent.frame())
#return(object)
}
