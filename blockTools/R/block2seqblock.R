block2seqblock <- function(block.obj, assg.obj, data, exact.restr = NULL, covar.restr = NULL, covar.order = NULL, trn = NULL, apstat = "mean", mtrim = 0.1, apmeth = "ktimes", kfac = 2, assgpr = c(0.5, 0.5), distance = NULL, datetime = NULL, orig, seed = NULL, file.name = "sbout.RData", verbose = FALSE){
  
  ## Store call used to run the block function originally as character string:
  block.call <- as.character(deparse(block.obj$call)) 
  ## Remove white space:
  block.call <- gsub("    ", "", block.call, fixed=T)   
  block.call <- paste(block.call, collapse="")
  block.call <- substr(block.call, 7, nchar(block.call)-1)
  block.call <- unlist(strsplit(block.call, " "))
  
  ## Read in id.vars from the provided block object
    first.id <- grep("id.vars", block.call) + 2
    if(length(grep("=", block.call[(first.id + 1):length(block.call)]))!=0){
      last.id <- grep("=", block.call[(first.id + 1):length(block.call)])[1] - 2 + first.id
    } else {last.id <- length(block.call)}
    id.vars <- block.call[first.id:last.id]
    ## Remove c() if it exists:
    if(substr(id.vars[1], 1, 2)=="c("){ 
      id.vars[1] <- substr(id.vars[1], 3, nchar(id.vars[1]))
      ## Drop 2, one for right parenthesis and one for comma --- ),
      id.vars[length(id.vars)] <- substr(id.vars[length(id.vars)], 1, nchar(id.vars[length(id.vars)])-2) 
    }
    for(i in 1:length(id.vars)){
      ## Remove trailing commas 
      if(substr(id.vars[i], nchar(id.vars[i]), nchar(id.vars[i]))==","){ 
        id.vars[i] <- substr(id.vars[i], 1, nchar(id.vars[i])-1)
      }
      ## Remove extra quotes around the word
      if(substr(id.vars[i], 1, 1)=="\""){ 
        id.vars[i] <- substr(id.vars[i], 2, nchar(id.vars[i])-1)
      }
    }
  
  ## Read in exact.vars from the provided block object or set them to NULL
    if(length(grep("groups", block.call)) != 0){
      first.exact <- grep("groups", block.call) + 2
      if(length(grep("=", block.call[(first.exact + 1):length(block.call)]))!=0){
        last.exact <- grep("=", block.call[(first.exact + 1):length(block.call)])[1] - 2 + first.exact
      } else {last.exact <- length(block.call)}
      exact.vars <- block.call[first.exact:last.exact]
      ## Remove c() if it exists
      if(substr(exact.vars[1], 1, 2) == "c("){ 
        exact.vars[1] <- substr(exact.vars[1], 3, nchar(exact.vars[1]))
        exact.vars[length(exact.vars)] <- substr(exact.vars[length(exact.vars)], 1, nchar(exact.vars[length(exact.vars)])-2)
      }
      for(i in 1:length(exact.vars)){
        if(substr(exact.vars[i], nchar(exact.vars[i]), nchar(exact.vars[i])) == ","){
          exact.vars[i] <- substr(exact.vars[i], 1, nchar(exact.vars[i])-1)
        }
        if(substr(exact.vars[i], 1, 1) == "\""){
          exact.vars[i] <- substr(exact.vars[i], 2, nchar(exact.vars[i])-1)
        }
      }
    } else{exact.vars <- NULL}
  
  ## Read in covar.vars from the provided block object or define all non-id-variables as covar.vars   
    if(length(grep("block.vars", block.call)) != 0){
      first.covar <- grep("block.vars", block.call) + 2
      if(length(grep("=", block.call[(first.covar + 1):length(block.call)]))!=0){
        last.covar <- grep("=", block.call[(first.covar + 1):length(block.call)])[1] - 2 + first.covar
      } else {
        last.covar <- length(block.call)
      }
      covar.vars <- block.call[first.covar:last.covar]
      ## Remove c() if it exists
      if(substr(covar.vars[1], 1, 2) == "c("){ 
        covar.vars[1] <- substr(covar.vars[1], 3, nchar(covar.vars[1]))
        covar.vars[length(covar.vars)] <- substr(covar.vars[length(covar.vars)], 1, nchar(covar.vars[length(covar.vars)])-2)
      }
      for(i in 1:length(covar.vars)){
        if(substr(covar.vars[i], nchar(covar.vars[i]), nchar(covar.vars[i])) == ","){
          covar.vars[i] <- substr(covar.vars[i], 1, nchar(covar.vars[i])-1)
        }
        if(last.covar < length(block.call)){
          if(substr(covar.vars[i], 1, 1) == "\""){
            covar.vars[i] <- substr(covar.vars[i], 2, nchar(covar.vars[i])-1)
          }
        } else if(last.covar == length(block.call)){
          if(i < length(covar.vars)){
            if(substr(covar.vars[i], 1, 1) == "\""){
              covar.vars[i] <- substr(covar.vars[i], 2, nchar(covar.vars[i])-1)
            }
          }
          if(i == length(covar.vars)){
            if(substr(covar.vars[i], 1, 1) == "\""){
              covar.vars[i] <- substr(covar.vars[i], 2, nchar(covar.vars[i]))
            }
          }
        }
      }
    } else{
      covar.vars <- colnames(data)[which(colnames(data) != id.vars)]
      covar.vars <- covar.vars[which(covar.vars != exact.vars)]    
    }

  ## Create trn object
  if(is.null(trn)){
    if(length(grep("n.tr", block.call)) != 0){
      temp <- grep("n.tr", block.call) + 2
      numb.tr <- block.call[temp]
      if(substr(numb.tr[1], 1, 2)=="c("){ #get rid of c() if it exists
        numb.tr[1] <- substr(numb.tr[1], 3, nchar(numb.tr[1]))
        ## Drop 2, one for right parenthesis and one for comma --- ),
        numb.tr[length(numb.tr)] <- substr(numb.tr[length(numb.tr)], 1, nchar(numb.tr[length(numb.tr)])-2) 
      }
      if(substr(numb.tr, nchar(numb.tr), nchar(numb.tr))==","){
        numb.tr <- substr(numb.tr, 1, nchar(numb.tr)-1)
      }
      if(substr(numb.tr, 1, 1)=="\""){
        numb.tr <- substr(numb.tr, 2, nchar(numb.tr)-1)
      }
      numb.tr <- as.numeric(numb.tr)
      trn <- paste("Treatment", 1:numb.tr)
    } else{
      numb.tr <- 2
      trn <- paste("Treatment", 1:numb.tr)
    }
  }
  
  ## Create distance object
  if(is.null(distance)){
    if(length(grep("distance", block.call))!=0){
      temporary <- grep("distance", block.call) + 2
      distance <- block.call[temporary]
      ## Remove c() if it exists
      if(substr(distance[1], 1, 2)=="c("){ 
        distance[1] <- substr(distance[1], 3, nchar(distance[1]))
        ## Drop 2, one for right parenthesis and one for comma --- ),
        distance[length(distance)] <- substr(distance[length(distance)], 1, nchar(distance[length(distance)])-2) 
      }
      if(substr(distance, nchar(distance), nchar(distance)) == ","){
        distance <- substr(distance, 1, nchar(distance)-1)
      }
      if(substr(distance, 1, 1) == "\""){
        distance <- substr(distance, 2, nchar(distance)-1)
      }
    }
    else{distance <- "mahalanobis"}
  }
  
  ## Set seed and create x and orig objects 
  if(!is.null(seed)){
    set.seed(seed)
  }
  x1 <- data[, which(colnames(data) %in% id.vars)]
  x2 <- data[, which(colnames(data) %in% exact.vars)]
  x3 <- data[, which(colnames(data) %in% covar.vars)]
  x <- as.data.frame(cbind(x1, x2, x3), stringsAsFactors = FALSE) 
  colnames(x) <- c(id.vars, exact.vars, covar.vars)
  ## Draw treatment assignments from original assignment object:
  assg <- assg.obj[[1]]
  treat.assign <- vector(mode = "list", length = length(trn))
  for(i in 1:length(assg)){
    for(j in 1:length(trn)){
      treat.assign[[j]] <- append(treat.assign[[j]], as.vector(assg[[i]][,j]))
    } 
  }
  ## treat.assign is a list:
  ## first element contains id values for all units that were assigned to TR 1,
  ## second element contains id values for all units assigned to TR 2, etc. 
  treatment <- vector(length = nrow(data))
  for(i in 1:length(trn)){
    treatment[which(x[,colnames(x) == id.vars[1]] %in% treat.assign[[i]])] <- trn[i]
  }
  x <- cbind(x, treatment)
  colnames(x)[ncol(x)] <- "Tr"
  
  for(i in 1:ncol(x)){
    if(class(x[,i]) == "factor"){
      x[,i] <- as.character(x[,i])
    }
  }
  
  bdata <- list()
  bdata$x <- x
  bdata$nid <- id.vars
  bdata$nex <- exact.vars
  bdata$ncv <- covar.vars
  bdata$rex <- exact.restr ## EXACT restricted values list
  if(!is.null(bdata$rex)){
    names(bdata$rex) <- exact.vars
  }
  bdata$rcv <- covar.restr  ## BLOCK restricted values list
  if(!is.null(bdata$rcv)){
    names(bdata$rcv) <- covar.vars
  }
  bdata$ocv <- covar.order   ## block covariates order
  bdata$trn <- trn
  bdata$apstat <- apstat 
  bdata$mtrim <- mtrim
  bdata$apmeth <- apmeth
  bdata$kfac <- kfac ## multiple for method 'ktimes'
  bdata$assgpr <- assgpr
  bdata$distance <- distance
#   bdata$trd <- tr.dist
#   bdata$tr.sort <- tr.sort
#   bdata$p <- p
#   bdata$trcount <- tr.counts
  if(!is.null(datetime)){
    bdata$datetime <- rep(datetime, length=nrow(x))  
  } else(bdata$datetime <- rep(NA, length=nrow(x)))
  bdata$orig <- x

  ## Write bdata to file
  save(bdata, file = file.name)
  if(verbose == TRUE){
    cat("Units' data stored as file ", file.name, ".\nThe current working directory is ", getwd(), "\n", sep="")
    cat("The new data as entered:\n")
    print(x)
  }  
  return(bdata)	
}