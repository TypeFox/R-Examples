getsyn <- function(matrix){
 
# matrix            <- code(matrix) #PGE TOOLBOX

if(length(colnames(matrix))==0){colnames(matrix) <- 1:dim(matrix)[2]}


   # matrixone          <- deletecodongaps(matrix) #### changed codonise64
   
    matrixone <- matrix
  
if(length(matrixone)==0){return(list(mono=as.matrix(NaN),syn=as.matrix(NaN),nonsyn=as.matrix(NaN),monoid=NaN,synid=NaN,nonsynid=NaN,Codons=as.list(NaN)))}

matrix                 <- codonise64(matrixone)

#print(matrix)

monoid   <- NULL
synid    <- NULL
nonsynid <- NULL

# -------------------------------------
acc      <- 1
Codons   <- vector("list",dim(matrix)[2])
# -------------------------------------

erg <- apply(matrix,2,function(vek){
 
 uniquevek     <- unique(vek) 
 size          <- length(uniquevek)

 if(size==1){
    # All equal
   Codons[[acc]] <<- uniquevek
   acc           <<- acc + 1 
   return(1) # monomorph
   
 }

 if(size>1){

    Codons[[acc]] <<- uniquevek
    acc           <<- acc + 1

    cc   <- codontable()
    id65 <- which(uniquevek==65) # gap or unknown position
    if(length(id65)>0){vek  <- uniquevek[-id65]}
    tt   <- unique(cc$Protein[1,vek]) # Standard Codierung
    size <- length(tt)
   
    if(size==1){
       return(2) # synonym
    }

    if(size>1){
      return(3) # nonsynonym
    }
 }

}) 
rm(matrix)

mono   <- which(erg==1) # monomorph
syn    <- which(erg==2) # synonym
nonsyn <- which(erg==3) # nonsynonym

# reproduziere die syn/nonsyn/mono Stellen

# monoid      <- vector(,length(mono)*3)
# synid       <- vector(,length(syn)*3)
# nonsynid    <- vector(,length(nonsyn)*3) 
 
#if(length(monoid)>0){
#  yy <- 1
#  for(xx in seq(1,length(monoid),by=3)){
#   monoid[xx:(xx+2)] <- c(3*mono[yy]-2,3*mono[yy]-1,3*mono[yy])   
#   yy <- yy + 1 
#  }
#}

if(length(mono)>0){
 monoid     <- c(3*mono-2,3*mono-1,3*mono)
 monoid     <- sort(monoid)
}

#if(length(synid)>0){
# yy <- 1
# for(xx in seq(1,length(synid),by=3)){
#  synid[xx:(xx+2)] <- c(3*syn[yy]-2,3*syn[yy]-1,3*syn[yy])   
#  yy <- yy + 1
# }
#}

if(length(syn)>0){
  synid     <- c(3*syn-2,3*syn-1,3*syn)
  synid     <- sort(synid)
}


#if(length(nonsynid)>0){
# yy <- 1
# for(xx in seq(1,length(nonsynid),by=3)){
#  nonsynid[xx:(xx+2)] <- c(3*nonsyn[yy]-2,3*nonsyn[yy]-1,3*nonsyn[yy])   
#  yy <- yy + 1 
# }
#}

if(length(nonsyn)>0){
 nonsynid     <- c(3*nonsyn-2,3*nonsyn-1,3*nonsyn)
 nonsynid     <- sort(nonsynid)
}



# monoid   <- monoid  [!is.na(monoid)]
# synid    <- synid   [!is.na(synid)]
# nonsynid <- nonsynid[!is.na(nonsynid)] 

# if(length(monoid)>1)   {  m <- decodonise64(matrix[,mono,drop=F]) }else{m   <- NA}
# if(length(synid) >1)   {  s <- decodonise64(matrix[,syn,drop=F])  }else{s   <- NA}
# if(length(nonsynid)>1) { ns <- decodonise64(matrix[,nonsyn,drop=F])}else{ns <- NA}

# m  <- as.matrix(NaN)
# s  <- as.matrix(NaN)
# ns <- as.matrix(NaN)

monoid   <- as.numeric(colnames(matrixone[,monoid]))
synid    <- as.numeric(colnames(matrixone[,synid]))
nonsynid <- as.numeric(colnames(matrixone[,nonsynid]))


return(list(monoid=monoid,synid=synid,nonsynid=nonsynid,Codons=Codons))
}


