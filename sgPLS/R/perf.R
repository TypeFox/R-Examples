# ---------------------------------------------------
# perf for sPLS object
# ---------------------------------------------------

perf.sPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        
        #-- spls --#
        spls.res = sPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY)     ## change
        
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), selectVar(spls.res, comp = k)$name.X)
          featuresY[[k]] = c(unlist(featuresY[[k]]), selectVar(spls.res, comp = k)$name.Y)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = selectVar(object, comp = k)$value.X
      features.finalY[[k]] = selectVar(object, comp = k)$value.Y
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

# ---------------------------------------------------
# perf for gPLS object
# ---------------------------------------------------

perf.gPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    ind.block.x = object$ind.block.x
    ind.block.y = object$ind.block.y
    
    
    
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        #-- gpls --#
        spls.res = gPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY,ind.block.x=ind.block.x,ind.block.y=ind.block.y)     ## change
        
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), selectVar(spls.res, comp = k)$name.X)
          featuresY[[k]] = c(unlist(featuresY[[k]]), selectVar(spls.res, comp = k)$name.Y)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict.gPLS(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = selectVar(object, comp = k)$value.X
      features.finalY[[k]] = selectVar(object, comp = k)$value.Y
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }


# ---------------------------------------------------
# perf for sgPLS object
# ---------------------------------------------------




perf.sgPLS <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,setseed=1,
           ...)
  {
    set.seed(setseed)
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from spls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   
    keepY = object$keepY   
    
    mode = object$mode
    ncomp = object$ncomp
    ind.block.x = object$ind.block.x
    ind.block.y = object$ind.block.y
    alpha.x = object$alpha.x
    alpha.y = object$alpha.y
    upper.lambda = object$upper.lambda
    
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        #-- sgpls --#
        spls.res = sgPLS(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY,ind.block.x=ind.block.x,ind.block.y=ind.block.y,alpha.x=alpha.x,alpha.y=alpha.y,upper.lambda=upper.lambda)     ## change
        #  Sparse.Group.spls.BP(X,Y,ncomp=1,mode="regression",max.iter=500,tol=1e-06,keepX=c(4),keepY=c(4),ind.block.y=ind.block.y,ind.block.x=ind.block.x,alpha.x=0.05,alpha.y=0.95,upper.lambda=1000000000)
        #res.sparse <- s
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), selectVar(spls.res, comp = k)$name.X)
          featuresY[[k]] = c(unlist(featuresY[[k]]), selectVar(spls.res, comp = k)$name.Y)
        }
        
        
        #if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict.sgPLS(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = selectVar(object, comp = k)$value.X
      features.finalY[[k]] = selectVar(object, comp = k)$value.Y
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

