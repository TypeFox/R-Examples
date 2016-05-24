create.PopGenome.method <- function(function.name, population.specific=TRUE){

zz <- file(paste(function.name,".R",sep=""),"w")

if(population.specific){

cat("setGeneric(\"",function.name,"\",function(object,new.populations=FALSE,subsites=FALSE) standardGeneric(\"",
function.name,"\"))\n",sep="",file=zz)

cat("setMethod(\"",function.name,"\" , \"GENOME\" ,\n " ,sep="",file=zz)
cat("function(object,new.populations,subsites){ \n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)

cat("# BEGIN: DO NOT EDIT\n",file=zz)

cat("region.names                       <- object@region.names\n",file=zz)
cat("n.region.names                     <- length(region.names)\n",file=zz)

cat("\n",file=zz)
cat("\n",file=zz)

cat("if(missing(new.populations)){\n",file=zz)
cat(" npops          <- length(object@populations)\n",file=zz)
cat(" NEWPOP         <- FALSE\n",file=zz)
cat("}else{\n",file=zz)
cat(" NEWPOP         <- TRUE\n",file=zz)
cat(" populations    <- vector(\"list\",npops)\n",file=zz)
cat(" npops          <- length(new.populations)\n",file=zz)
cat("}\n",file=zz)


cat("\n",file=zz)
cat("\n",file=zz)

 cat("new.value          <- matrix(,n.region.names,npops)\n",file=zz)
 
cat("\n",file=zz)
cat("\n",file=zz)

 cat("colnames(new.value) <- paste(\"pop\",1:npops) \n",file=zz)
 cat("rownames(new.value) <- region.names\n",file=zz)

cat("\n",file=zz)

#cat("progr <- progressBar()\n",file=zz)

cat("\n",file=zz)
cat("for(xx in 1:n.region.names){\n",file=zz)
cat("\n",file=zz)

cat("bial <- popGetBial(object,xx)\n",file=zz)
cat("\n",file=zz)

cat("if(subsites[1]!=FALSE){\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"transitions\" & length(bial!=0)){\n",file=zz)
cat("            tran       <- which(object@region.data@transitions[[xx]]==TRUE)\n",file=zz)
cat("            bial       <- bial[,tran,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"transversions\" & length(bial!=0)){\n",file=zz)
cat("            transv     <- which(object@region.data@transitions[[xx]]==FALSE)\n",file=zz)
cat("            bial       <- bial[,transv,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"syn\" & length(bial!=0)){\n",file=zz)
cat("            syn        <- which(object@region.data@synonymous[[xx]]==TRUE)\n",file=zz)
cat("            bial       <- bial[,syn,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"nonsyn\" & length(bial!=0)){\n",file=zz)
cat("            nonsyn     <- which(object@region.data@synonymous[[xx]]==FALSE)\n",file=zz)
cat("            bial       <- bial[,nonsyn,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)


cat("     if(subsites==\"intron\" & length(bial!=0)){\n",file=zz)
cat("           intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)\n",file=zz)
cat("            if(length(intron)==0){\n",file=zz)
cat("              intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]\n",file=zz)	  
cat("            }\n",file=zz)
cat("           bial          <- bial[,intron,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)


cat("     if(subsites==\"utr\" & length(bial!=0)){\n",file=zz)
cat("          utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)\n",file=zz)
cat("          bial          <- bial[,utr,drop=FALSE]\n",file=zz)  
cat("     }\n",file=zz)

cat("     if(subsites==\"exon\" & length(bial!=0)){\n",file=zz)
cat("          exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)\n",file=zz)
cat("          bial           <- bial[,exon,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("     if(subsites==\"coding\" & length(bial!=0)){\n",file=zz)
cat("          coding           <- which(!is.na(object@region.data@synonymous[[xx]]))\n",file=zz)
cat("          bial             <- bial[,coding,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("     if(subsites==\"gene\" & length(bial!=0)){\n",file=zz)
cat("         gene             <- which(object@region.data@GeneSNPS[[xx]])\n",file=zz)
cat("         bial             <- bial[,gene,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("}\n",file=zz)

cat("\n",file=zz)

cat(" if(length(bial)!=0){ # if a biallelic position exists \n",file=zz)  
cat("\n",file=zz)       
cat("    if(NEWPOP){\n",file=zz) 
cat("       for(yy in 1:npops){\n",file=zz)
cat("           if(is.character(new.populations[[yy]])){\n",file=zz)
cat("              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))\n",file=zz)
cat("              naids             <- which(!is.na(populations[[yy]]))\n",file=zz)
cat("              populations[[yy]] <- populations[[yy]][naids]\n",file=zz)
cat("           }else{\n",file=zz)
cat("              populations[[yy]] <- new.populations[[yy]]\n",file=zz)
cat("              ids               <- which(populations[[yy]]>dim(bial)[1])\n",file=zz)
cat("              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}\n",file=zz)   
cat("           }\n",file=zz)   
cat("       }\n",file=zz)
cat("\n",file=zz)       
cat("       temp         <- delNULLpop(populations)\n",file=zz)
cat("       populations  <- temp$Populations\n",file=zz)
cat("       popmissing   <- temp$popmissing\n",file=zz)
cat("\n",file=zz)       
cat("       if(length(populations)==0){next}\n",file=zz)
cat("\n",file=zz)       
cat("    }else{\n",file=zz)
cat("     populations <- object@region.data@populations[[xx]]\n",file=zz)
cat("    }\n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)    
cat("    if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}\n",file=zz)
cat("    if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}\n",file=zz)
    
cat("\n",file=zz)
cat("# END: DO NOT EDIT\n",file=zz)
cat("\n",file=zz)
cat("# define here your own function in the PopGenome framework.\n",file=zz)
cat("# default: bial (biallelic matrix), populations (the defined populations)\n",file=zz)
cat("# ... chose here everthing you want from the GENOME class\n",file=zz)
cat("\n",file=zz)
     
cat("new.value[xx,respop]         <- your_intern_function(bial,populations,...)\n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)
     
    
#cat("    progr <- progressBar(xx,n.region.names, progr)\n",file=zz)
 
cat(" }\n",file=zz)
cat("}\n",file=zz)
cat("\n",file=zz)
  
cat("  return(new.value)\n",file=zz)
cat(" })\n",file=zz)

}# END of population specific 

if(!population.specific){

cat("setGeneric( \"",function.name,"\",function(object,new.populations=FALSE,subsites=FALSE) standardGeneric(\"",function.name,"\" )) \n",sep="",file=zz)

cat("setMethod(\"",function.name,"\" , \"GENOME\" ,\n ",sep="",file=zz)
cat("function(object,new.populations,subsites){ \n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)

cat("# BEGIN: DO NOT EDIT\n",file=zz)

cat("region.names                       <- object@region.names\n",file=zz)
cat("n.region.names                     <- length(region.names)\n",file=zz)

cat("\n",file=zz)
cat("\n",file=zz)

cat("if(missing(new.populations)){\n",file=zz)
cat(" npops          <- length(object@populations)\n",file=zz)
cat(" NEWPOP         <- FALSE\n",file=zz)
cat("}else{\n",file=zz)
cat(" NEWPOP         <- TRUE\n",file=zz)
cat(" populations    <- vector(\"list\",npops)\n",file=zz)
cat(" npops          <- length(new.populations)\n",file=zz)
cat("}\n",file=zz)


cat("\n",file=zz)
cat("\n",file=zz)

 cat("new.value          <- matrix(,n.region.names,1)\n",file=zz)
 
cat("\n",file=zz)
cat("\n",file=zz)

 #cat("colnames(new.value) <- paste(\"pop\",1:npops) \n",file=zz)
 cat("rownames(new.value) <- region.names\n",file=zz)

cat("\n",file=zz)

#cat("progr <- progressBar()\n",file=zz)

cat("\n",file=zz)
cat("for(xx in 1:n.region.names){\n",file=zz)
cat("\n",file=zz)

cat("bial <- popGetBial(object,xx)\n",file=zz)
cat("\n",file=zz)

cat("if(subsites[1]!=FALSE){\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"transitions\" & length(bial!=0)){\n",file=zz)
cat("            tran       <- which(object@region.data@transitions[[xx]]==TRUE)\n",file=zz)
cat("            bial       <- bial[,tran,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"transversions\" & length(bial!=0)){\n",file=zz)
cat("            transv     <- which(object@region.data@transitions[[xx]]==FALSE)\n",file=zz)
cat("            bial       <- bial[,transv,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"syn\" & length(bial!=0)){\n",file=zz)
cat("            syn        <- which(object@region.data@synonymous[[xx]]==TRUE)\n",file=zz)
cat("            bial       <- bial[,syn,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)

cat("     if(subsites==\"nonsyn\" & length(bial!=0)){\n",file=zz)
cat("            nonsyn     <- which(object@region.data@synonymous[[xx]]==FALSE)\n",file=zz)
cat("            bial       <- bial[,nonsyn,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)
cat("\n",file=zz)


cat("     if(subsites==\"intron\" & length(bial!=0)){\n",file=zz)
cat("           intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)\n",file=zz)
cat("            if(length(intron)==0){\n",file=zz)
cat("              intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]\n",file=zz)	  
cat("            }\n",file=zz)
cat("           bial          <- bial[,intron,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)


cat("     if(subsites==\"utr\" & length(bial!=0)){\n",file=zz)
cat("          utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)\n",file=zz)
cat("          bial          <- bial[,utr,drop=FALSE]\n",file=zz)  
cat("     }\n",file=zz)

cat("     if(subsites==\"exon\" & length(bial!=0)){\n",file=zz)
cat("          exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)\n",file=zz)
cat("          bial           <- bial[,exon,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("     if(subsites==\"coding\" & length(bial!=0)){\n",file=zz)
cat("          coding           <- which(!is.na(object@region.data@synonymous[[xx]]))\n",file=zz)
cat("          bial             <- bial[,coding,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("     if(subsites==\"gene\" & length(bial!=0)){\n",file=zz)
cat("         gene             <- which(object@region.data@GeneSNPS[[xx]])\n",file=zz)
cat("         bial             <- bial[,gene,drop=FALSE]\n",file=zz)
cat("     }\n",file=zz)

cat("}\n",file=zz)

cat("\n",file=zz)

cat(" if(length(bial)!=0){ # if a biallelic position exists \n",file=zz)  
cat("\n",file=zz)       
cat("    if(NEWPOP){\n",file=zz) 
cat("       for(yy in 1:npops){\n",file=zz)
cat("           if(is.character(new.populations[[yy]])){\n",file=zz)
cat("              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))\n",file=zz)
cat("              naids             <- which(!is.na(populations[[yy]]))\n",file=zz)
cat("              populations[[yy]] <- populations[[yy]][naids]\n",file=zz)
cat("           }else{\n",file=zz)
cat("              populations[[yy]] <- new.populations[[yy]]\n",file=zz)
cat("              ids               <- which(populations[[yy]]>dim(bial)[1])\n",file=zz)
cat("              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}\n",file=zz)   
cat("           }\n",file=zz)   
cat("       }\n",file=zz)
cat("\n",file=zz)       
cat("       temp         <- delNULLpop(populations)\n",file=zz)
cat("       populations  <- temp$Populations\n",file=zz)
cat("       popmissing   <- temp$popmissing\n",file=zz)
cat("\n",file=zz)       
cat("       if(length(populations)==0){next}\n",file=zz)
cat("\n",file=zz)       
cat("    }else{\n",file=zz)
cat("     populations <- object@region.data@populations[[xx]]\n",file=zz)
cat("    }\n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)    
#cat("    if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}\n",file=zz)
#cat("    if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)#[-popmissing]}else{respop <- 1:npops}}\n",file=zz)
#    
cat("\n",file=zz)
cat("# END: DO NOT EDIT\n",file=zz)
cat("\n",file=zz)
cat("# define here your own function in the PopGenome framework.\n",file=zz)
cat("# default: bial (biallelic matrix), populations (the defined populations)\n",file=zz)
cat("# ... chose here everthing you want from the GENOME class\n",file=zz)
cat("\n",file=zz)
     
cat("new.value[xx]         <- your_intern_function(bial,populations,...)\n",file=zz)
cat("\n",file=zz)
cat("\n",file=zz)
     
    
#cat("    progr <- progressBar(xx,n.region.names, progr)\n",file=zz)
 
cat(" }\n",file=zz)
cat("}\n",file=zz)
cat("\n",file=zz)
  
cat("  return(new.value)\n",file=zz)
cat(" })\n",file=zz)


}

close(zz)

}
