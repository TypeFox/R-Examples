mmer2 <- function(fixed=NULL, random=NULL, G=NULL, R=NULL, method="EM", REML=TRUE, iters=50, draw=FALSE, init=NULL, data=NULL, family=gaussian, silent=FALSE, constraint=TRUE, sherman=FALSE, MTG2=FALSE, gss=FALSE){
  if(!is.null(G)){
    cat("With var-cov structures (G) present you may want to try the AI algorithm.\n\n")
  }

  yvar <- gsub(" ", "", as.character(fixed[2]))
  xvar <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+]")[[1]])
  if(!is.null(random)){
    zvar <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    varsss <- c(xvar,zvar)
  }else{varsss <- xvar}
  varsss <- varsss[which(varsss != "1")]
  doto <- grep(":",varsss)
  if(length(doto) >0){varsss <- varsss[-doto]}
  good <- list()
  for(i in 1:length(varsss)){
    xxx <- data[,varsss[i]]
    good[[i]] <- which(!is.na(xxx))
  }
  
  if(length(varsss) == 1){
    keep <- as.vector(unlist(good))
  }else{
    keep <- good[[1]]
    for(j in 1:length(good)){
      keep <- intersect(keep,good[[j]])
    }
  }
  data <- data[keep,] # only good data in model
  if(length(xvar %in% "1") == 1){
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c(xvar), collapse="+"))), data=data))
  }else{
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=data))
  }
  

  if(!is.null(random)){ # ===========

    
    yvars <- data[,yvar] #mox$family$linkfun(mox$y)
    #
    if(is.null(random)){
      Z=NULL
    }else{
      Z <- list()
      for(i in 1:length(zvar)){ # do it factor
        ## 
        vara <- zvar[i]
        data2 <- data.frame(apply(data,2,as.factor))
        zi <- model.matrix((as.formula(paste("~ -1 + ", vara))), data=data2)
        if(dim(zi)[2] == length(yvars)){ # lmer error
          stop("Error: number of levels of each grouping factor must be < number of observations")
        }else{ #
          
          ## 
          ww <- which(names(G) %in% vara)
          if(length(ww) > 0){#
            ki <- G[[ww]]
          }else{ # 
            ki <- diag(dim(zi)[2])
          }
          elem <- list(Z=zi, K=ki)
          Z[[i]] <- elem
        }
      }
      names(Z) <- zvar
      ### 
      res <- mmer(y=yvars, X=X, Z=Z, R=R, method=method, REML=REML, iters=iters, draw=draw, init=init, silent=silent, constraint=constraint, sherman=sherman, MTG2=MTG2, gss=gss)
      #rownames(res$var.comp) <- c(zvar,"Error")
    }
  }else{ # =
    res <- glm(yvars~X, family=family)
  }
  class(res)<-c("mmer")
  return(res)
}
