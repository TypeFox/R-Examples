


## table to dataframe (with names!)
tab2df <- function(x,...){

  is.tabrix <- class(x)[1] %in% c('table','matrix')
  row.nms <- rownames(x)
  
  if(length(dim(x))>1){ #NCOL(x)>1){
    clm.nms <- colnames(x)
    if(is.tabrix)
      x <- lapply(clm.nms, function(i) x[,i]) 
  }else{ # is a vector                                     
    if(is.tabrix){ # single column but matrix output
      clm.nms <- row.nms  # assume we want a 2 clmn dataframe
      row.nms <- NULL
    }else{   # vectors
      clm.nms <- NULL
      if(!is.null(names(x))) #named vector
        row.nms <- names(x) 
      else    # unnamed vector
        row.nms <- 1:length(x)
   }
  }
  
  if(!is.null(clm.nms)){
    x <- as.list(x)
    names(x) <- clm.nms
  }
  df <- data.frame(x,...)
  rownames(df) <- row.nms
  return(df)
  
}




nv <- function(x, name){

  if(class(x)=='data.frame'){
    v <- x[,name[1]]
    if(length(name)==2)
      names(v) <- x[,name[2]]
    else
      names(v) <- rownames(x)
  }else{
    if(NCOL(x)!=1)
      stop('x must be unidimentional')
    v <- x
    names(v) <- name
  }
  v
}


pct <- function(df, clmns){
   for(clmn in clmns)
    df[,paste(clmn,'pct',sep='.')] <- df[,clmn]/sum(df[,clmn])
  return(df) 
}

rerowname <- function(df, old='NA', new ='unknown'){
  tmp.rn <- rownames(df)
  tmp.rn[tmp.rn==old] <- new
  rownames(df) <- tmp.rn
  df
}
