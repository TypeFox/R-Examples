
createInput <- function(y, x, k, z, propscore, control, measurement, data, 
                        fixed.cell, fixed.z, missing, se, bootstrap, mimic,
                        interactions, ids, weights, homoscedasticity,
                        add){
  
  d <- data
  latentz <- z[which(!z %in% names(data))]
  vnames <- list(y=y,x=x,k=k,z=z,propscore=propscore,latentz=latentz)  
  
  ## treatment variable
  if(!is.factor(d[,x])){    
    d[,x] <- as.factor(d[,x])  
  }
  stopifnot(length(levels(d[,x])) <= 10) # test if it works for > 10 (problems with subscripts?)
  
  d[,x] <- relevel(d[,x], control)
  levels.x.original <- levels(d[,x])
  levels(d[,x]) <- paste(0:(length(levels(d[,x]))-1))  
  
  ## categorical covariates
  levels.k.original <- vector("list",length(k))
  names(levels.k.original) <- k
  
  if(!is.null(k)){    
    for(i in 1:length(k)){
      d[,k[i]] <- as.factor(d[,k[i]])
      levels.k.original[[i]] <- levels(d[,k[i]])
      levels(d[,k[i]]) <- paste(0:(length(levels(d[,k[i]]))-1))
    }    
  }
  
  ## unfolded k variable
  levels.kstar.original <- vector("character")
  if(!is.null(k)){
    if(length(k)>1){
      d$kstar <- apply(d[,k],1,paste,collapse="")
      d$kstar <- as.factor(d$kstar)    
    }else{
      d$kstar <- d[,k]
    }
    levels.kstar.original <- levels(d$kstar)
    levels(d$kstar) <- paste(0:(length(levels(d$kstar))-1))
    
    ## check for empty cells
    if(any(table(d$kstar, d[,x]) == 0)){
      stop("EffectLiteR error: Empty cells are currently not allowed.")
    }    
    
  }else{
    d$kstar <- NULL
  }
  
  ## nk
  nk <- 1L
  if(!is.null(k)){
    nk <- length(levels(d$kstar))
  }
  
  ## ng
  ng <- length(levels(d[,x]))  
  
  ## nz
  nz <- length(vnames$z)
  
  ## check for too many cells
  if((nk>10 & ng>10) || (nk>10 & nz>10) || (ng>10 & nz>10)){
    stop("EffectLiteR error: Too many cells")
  }
  
  
  ## cell variable (xk-cells)
  if(!is.null(k)){
    cell <- expand.grid(k=levels(d$kstar), x=levels(d[,x]))
    cell <- with(cell, paste0(x,k))
    dsub <- cbind(d[,x],d$kstar) - 1 # use x=0,1,2... and k=0,1,2,... as labels
    d$cell <- apply(dsub, 1, function(x){
      missing_ind <- sum(is.na(x)) > 0
      if(missing_ind){
        return(NA)
      }else{
        return(paste(x, collapse=""))
      }
    }) 
    d$cell <- as.factor(d$cell)    
  }else{
    cell <- levels(d[,x])
    d$cell <- d[,x]
  }
  
  
  ## observed cell frequencies (fixed.cell only)
  if(!fixed.cell){
    observed.freq <- numeric(0)
    
  }else if(fixed.cell){
    N <- nrow(d)
    observed.freq <- c(table(d$cell)/N)
    
    if(!is.null(weights)){
      message("EffectLiteR message: The observed frequencies have been re-computed taking into account the survey weights.")
      weights_vector <- model.matrix(weights, d)
      if(ncol(weights_vector) > 2){stop("EffectLiteR error: Currently only support for one weights variable")}
      weights_vector <- weights_vector[,-1]
      observed.freq <- c(tapply(weights_vector, d$cell, sum))
      observed.freq <- observed.freq/sum(observed.freq) ## rescale to sum to one
    }
  }
    
  

  ## observed sample means for manifest covariates (fixed.z only)
  if(nz==0){
    sampmeanz <- matrix(nrow=0, ncol=0)
    
  }else if(!fixed.z){
    sampmeanz <- matrix(nrow=0, ncol=0)
    
  }else if(fixed.z){
      
    if(!fixed.cell){
      stop("EffectLiteR error: fixed.z=TRUE requires fixed.cell=TRUE")
      
    }else if(!identical(latentz, character(0))){
      stop("EffectLiteR error: fixed.z=TRUE does not work with latent covariates.")
    }
    
    sampmeanz <- NULL
    for (i in 1:nz) {
      namez <- z[i]
      tmp <- tapply(d[[namez]], d$cell, function(x){mean(x, na.rm=TRUE)})
      sampmeanz <- rbind(sampmeanz, tmp)
    }
    row.names(sampmeanz) <- z
    
  }
  
  
  ## add vlevels for created variables
  vlevels <- list(levels.x.original=levels.x.original,
                  levels.k.original=levels.k.original,
                  levels.kstar.original=levels.kstar.original,
                  x=levels(d[,x]),
                  kstar=levels(d$kstar),
                  cell=levels(d$cell))
  
  
  
  complexsurvey <- list(ids=ids, weights=weights)
  
  ## non-standard se only work with fixed group sizes
  if(se != "standard" & fixed.cell==FALSE){
    
    stop("EffectLiteR error: Non-standard SEs currently only work with fixed cell sizes. Please use fixed.cell=TRUE.")
    
  }
  
  
  res <- new("input",
             vnames=vnames, 
             vlevels=vlevels,
             ng=ng,
             nz=nz,
             nk=nk,
             control=control,
             data=d, 
             measurement=measurement,
             add=add,
             fixed.cell=fixed.cell,
             fixed.z=fixed.z,
             missing=missing,
             observed.freq=observed.freq,
             sampmeanz=sampmeanz,
             se=se,
             bootstrap=bootstrap,
             mimic=mimic,
             interactions=interactions,
             complexsurvey=complexsurvey,
             homoscedasticity=homoscedasticity
  )
  
  return(res)
}

