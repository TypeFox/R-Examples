
# check each population. calculate the segragating sites (positions)

get_segsites <- function(matrix_pol,populations){

if(missing(populations)){

populations      <- vector("list",1)
populations[[1]] <- 1:dim(matrix_pol)[1]

}

npops    <- length(populations)
Segpos   <- vector("list",npops)

for(xx in 1:npops){
   popmatrix <- matrix_pol[populations[[xx]],,drop=FALSE]
   erg <- apply(popmatrix,2,function(x){
          check   <- unique(x)
          nongaps <- !is.na(check)
          check   <- check[nongaps]
          if(length(check)==2){return(TRUE)}else{return(FALSE)}
          })
   Segpos[[xx]] <- which(erg)
   
}
return(Segpos)
}

#### Get SegSites FAST

get_segsites_FAST <- function(matrix_pol,populations){

if(missing(populations)){
populations      <- vector("list",1)
populations[[1]] <- 1:dim(matrix_pol)[1]
}

npops    <- length(populations)
Segpos   <- vector("list",npops)

for(xx in 1:npops){

   popmatrix    <- matrix_pol[populations[[xx]],,drop=FALSE]
   # site_length  <- dim(popmatrix)[1]

   erg <- apply(popmatrix,2,function(x){

         check       <- sum(x, na.rm=TRUE)
         site_length <- sum(!is.na(x))

          if(check==0 | check==site_length){return(FALSE)}else{return(TRUE)}
          })

   Segpos[[xx]] <- which(erg)
   
}
return(Segpos)
}
