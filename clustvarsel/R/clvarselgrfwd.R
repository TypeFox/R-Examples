#############################################################################
## Sequential & parallel backward greedy search
#############################################################################

clvarselgrfwd <- function(X, G = 1:9, 
                          emModels1 = c("E","V"), 
                          emModels2 = mclust.options("emModelNames"),
                          samp = FALSE, sampsize = 2000, 
                          hcModel = "VVV", allow.EEE = TRUE, forcetwo = TRUE, 
                          BIC.diff = 0, itermax = 100, 
                          parallel = FALSE, verbose = FALSE)
{
  X <- as.matrix(X)
  n <- nrow(X) # number of rows = number of observations
  d <- ncol(X) # number of columns = number of variables
  if(is.null(colnames(X))) 
    colnames(X) <- paste("X", 1:d, sep = "")
  G <- setdiff(G, 1)
  
  # If needed, sample the subset of observations for hierarchical clustering
  if(samp) { sub <- sample(1:n, min(sampsize,n), replace = FALSE) }
  else     { sub <- seq.int(1,n) }

  # Start "parallel backend" if needed
  parallel <- if(is.logical(parallel)) 
                { if(parallel) startParallel(parallel) else FALSE }
              else { startParallel(parallel) }
  on.exit(if(parallel)
          parallel::stopCluster(attr(parallel, "cluster")))
  # define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel) `%dopar%` else `%do%`
  i <- NULL # dummy to trick R CMD check 

  # First Step - selecting single variable
  hcModel1 <- if(any(grep("V", hcModel))) "V" else "E"
  if(verbose) print(paste("iter", 1))
  out <- foreach(i = 1:d) %DO% 
  {
    xBIC <- NULL
    # Fit the single variable cluster models
    try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1,
                       initialization = list(subset = sub)),
        silent = TRUE)
    # If we get all NA's from selected starting hierarchical values use "E"
    if((allow.EEE) & sum(is.finite(xBIC$BIC))==0)
      try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1, 
                         initialization = list(hcPairs = hcE(X[sub,i]),
                                               subset = sub)),
          silent = TRUE)
    # maxBIC is the maximum BIC over all clustering models fit
    maxBIC <- if(sum(is.finite(xBIC$BIC))==0) NA
              else max(xBIC$BIC[is.finite(xBIC$BIC)])
    # Fit and get BIC for a single component no-cluster normal model
    try(oneBIC <- Mclust(X[,i], G = 1, modelNames = emModels1,
                         initialization = list(subset = sub))$BIC[1],
        silent = TRUE)
    #
    return(c(maxBIC, oneBIC))
  }
  maxBIC <- sapply(out, "[[", 1)
  oneBIC <- sapply(out, "[[", 2)
  # Difference between maximum BIC for clustering and BIC for no clustering
  maxdiff <- maxBIC - oneBIC
  # Find the single variable with the biggest difference between 
  # clustering and no clustering
  m <- max(maxdiff[is.finite(maxdiff)])
  arg <- which(maxdiff==m,arr.ind=TRUE)[1]
  # This is our first selected variable/S is the matrix of currently selected
  # clustering variables
  S <- X[,arg,drop=FALSE]
  # BICS is the BIC value for the clustering model with the variable(s) in S
  BICS <- maxBIC[arg]
  # NS is the matrix of currently not selected variables
  NS <- X[,-arg,drop=FALSE]
  # info records the proposed variable, BIC for the S matrix and difference 
  # in BIC for clustering versus no clustering on S, whether it was an 
  # addition step and if it was accepted
  info <- data.frame(Var = colnames(S), BIC = BICS, BICdif = maxdiff[arg], 
                     Step = "Add", Decision = "Accepted",
                     stringsAsFactors = FALSE)
  if(verbose) print(info)
  if(verbose) print(paste("iter", 2))
  # Second Step - selecting second variable
  out <- foreach(i = 1:ncol(NS)) %DO%
  {
    # Calculate the BIC for the regression of the proposed variable 
    # on the variable in S
    regBIC <- BICreg(y = NS[,i], x = S)    
    # Fit the cluster model on the two variables
    sBIC <- NULL
    try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, 
                       modelNames = emModels2,
                       initialization = list(hcPairs = hc(hcModel, data = cbind(S,NS[,i])[sub,]), subset = sub)),
        silent = TRUE)
    # If we get all NA's from "VVV" starting hierarchical values use "EEE"
    if((allow.EEE) & sum(is.finite(sBIC$BIC))==0)
       try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, 
                          modelNames = emModels2,
                          initialization = list(hcPairs = hc("EEE", data = cbind(S,NS[,i])[sub,]), subset = sub)),
           silent = TRUE)
    # depBIC is the BIC for the clustering model with both variables
    depBIC <- if(sum(is.finite(sBIC$BIC))==0) NA 
              else max(sBIC$BIC[is.finite(sBIC$BIC)])
    #
    return(c(regBIC, depBIC))
  }
  regBIC <- sapply(out, "[[", 1)
  depBIC <- sapply(out, "[[", 2)
  # cindepBIC is the BIC for the clustering model on S and the 
  # regression model of the new variable on S
  cindepBIC <- regBIC + BICS
  # cdiff is the difference between BIC for the models with variables' being
  # clustering variables versus them being conditionally independent of 
  # the clustering
  cdiff <- depBIC - cindepBIC
  # Choose the variable with the largest difference
  m <- max(cdiff[is.finite(cdiff)])
  arg <- which(cdiff==m,arr.ind=TRUE)[1]

  # if forcetwo is true automatically add the best second variable, 
  # otherwise only add it if its difference is large than BIC.diff
  if(forcetwo || cdiff[arg] > BIC.diff)
    { k <- c(colnames(S),colnames(NS)[arg])
      nks <- c(colnames(NS)[-arg])
      BICS <- depBIC[arg]
      info <- rbind(info, c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Accepted"))
      S <- cbind(S,NS[,arg])
      NS <- as.matrix(NS[,-arg])
      colnames(S) <- k
      colnames(NS) <- nks
    } 
  else
    { info <- rbind(info, c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Rejected")) }

  if(verbose) print(info)

  criterion <- 1
  iter <- 0
  while((criterion == 1) & (iter < itermax))
  {
    iter <- iter+1
    check1 <- colnames(S)
    
    if(verbose) print(paste("iter", iter+2))
    
    # Adding step
    if(verbose) print("adding step")
    # For the special case where we have removed all the clustering 
    # variables/S is empty
    if(ncol(S)==0 || is.null(ncol(S)))
      { 
        # We simply choose the same variable as in the first step and check 
        # whether the difference between the BIC for clustering versus not 
        # clustering is positive or not
        m <- max(maxdiff[is.finite(maxdiff)])
        arg <- which(maxdiff==m,arr.ind=TRUE)[1]
        if(maxdiff[arg] > BIC.diff)
          {
            # if the difference is positive this variable is selected 
            # as a clustering variable
            S <- matrix(c(X[,arg]),n,1)
            BICS <- maxBIC[arg]
            colnames(S) <- colnames(X)[arg]
            NS <- as.matrix(X[,-arg])
            colnames(NS) <- colnames(X)[-arg]
            info <- rbind(info, c(colnames(S),BICS,maxdiff[arg],"Add","Accepted"))
          } 
        else
          { # if the difference is not > BIC.diff no clustering variables exist
            BICS <- NA
            info <- rbind(info,
                          c(colnames(X)[arg],BICS,maxdiff[arg],"Add","Rejected"))
          }
      }
    else
      { # Addition Step in general (for all cases except when S is empty)
        if(ncol(NS) != 0 & !is.null(ncol(NS)))
          { 
            out <- foreach(i = 1:ncol(NS)) %DO%
            {
              # Calculate the BIC for the regression of the proposed
              # variable on the variable(s) in S
              regBIC <- BICreg(y = NS[,i], x = S)
              #
              # Fit the cluster model on the S variables with the proposed 
              # variable for 2 to G groups 
              sBIC <- NULL                 
              try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, 
                                 modelNames = emModels2,
                                 initialization = list(hcPairs = hc(hcModel, data = cbind(S,NS[,i])[sub,]), subset = sub)),
                  silent = TRUE)
              # If we get all NA's from "VVV" starting hierarchical values use "EEE"
              if((allow.EEE) & (sum(is.finite(sBIC$BIC))==0))
                 try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, modelNames = emModels2,
                                    initialization = list(hcPairs = hc("EEE", data = cbind(S,NS[,i])[sub,]),
                                                          subset = sub)),
                     silent = TRUE)
              # depBIC is the BIC for the clustering model with both S
              # and proposed variable
              depBIC <- if(sum(is.finite(sBIC$BIC))==0) NA 
                        else max(sBIC$BIC[is.finite(sBIC$BIC)])
              #
              return(c(regBIC, depBIC))
            }
            regBIC <- sapply(out, "[[", 1)
            depBIC <- sapply(out, "[[", 2)
            # cindepBIC is the BIC for the clustering model on S and the
            # regression model of the new variable on S
            cindepBIC <- regBIC + BICS
            # cdiff is the difference between BIC for the models with  
            # variables' being clustering variables versus them being
            # conditionally independent of the clustering
            cdiff <- depBIC - cindepBIC
            # Choose the variable with the largest difference
            m <- max(cdiff[is.finite(cdiff)])
            arg <- which(cdiff==m,arr.ind=TRUE)[1]
            if(cdiff[arg] > BIC.diff)
              { # if this difference is positive add this variable to S 
                # and update the clustering model's BICS
                BICS <- depBIC[arg] 
                k <- c(colnames(S),colnames(NS)[arg])
                nks <- c(colnames(NS)[-arg])
                info <- rbind(info,
                              c(colnames(NS)[arg], BICS, cdiff[arg], 
                                "Add","Accepted")) 
                S <- cbind(S,NS[,arg])
                NS <- as.matrix(NS[,-arg])
                colnames(S) <- k
                colnames(NS) <- nks
              } 
            else
              {  
                info <- rbind(info,
                              c(colnames(NS)[arg], BICS, cdiff[arg], 
                                "Add","Rejected"))
              }
          }
      }
    # Removal Step for the special case where S contains only 
    # a single variable
    if(verbose) print("removing step")
    if(ncol(S) == 1)
      { 
        cdiff <- 0
        oneBIC <- NA
        try(oneBIC <- Mclust(as.matrix(S), G = 1, modelNames = emModels1,
                             initialization = list(hcPairs = hc(hcModel1, data = S[sub,]), 
                                                   subset = sub))$BIC[1],
            silent = TRUE)
        # Difference between maximum BIC for clustering and BIC 
        # for no clustering
        cdiff <- c(BICS - oneBIC)
        if(is.na(cdiff)) cdiff <- 0
        # Check if difference is less than BIC.diff
        if(cdiff <= BIC.diff)
          { # if negative remove the variable from S and set the BIC 
            # for the model to NA
            BICS <- NA
            info <- rbind(info, c(colnames(S), BICS, cdiff, 
                          "Remove","Accepted"))
            k <- c(colnames(NS),colnames(S))
            NS <- cbind(NS,S)
            S <- NULL
            colnames(NS) <- k
          } 
        else
          { # Otherwise leave S and BICS alone
            info <- rbind(info, c(colnames(S), BICS, cdiff, 
                          "Remove","Rejected"))
          }
      } 
    else
      { # Removal step in general (for all cases except when S is a single 
        # variable or empty)
        if(ncol(S) >= 2)
          {
            # Check if the data is at least 3 dimensional
            name <- if(ncol(S) > 2) emModels2 else emModels1
            hcname <- if(ncol(S) > 2) hcModel else hcModel1
            out <- foreach(i = 1:ncol(S)) %DO%
            {
              # Calculate the BIC for the regression of the proposed
              # variable from S on the other variable(s) in S
              regBIC <- BICreg(y = S[,i], x = S[,-i,drop=FALSE])             
              # Fit the cluster model on the S variables without the 
              # proposed variable for 2 to G groups 
              sBIC <- NULL
              try(sBIC <- Mclust(as.matrix(S[,-i]), G = G, 
                                 modelNames = name,
                                 initialization = list(hcPairs = hc(hcname, data = cbind(S,NS[,i])[sub,]), 
                                                        subset = sub)),
                   silent = TRUE)
              # If we get all NA's from "VVV" starting hierarchical values 
              # use "EEE"
              if((allow.EEE) & ncol(S)>=3 & sum(is.finite(sBIC$BIC))==0)
                { try(sBIC <- Mclust(S[,-i], G = G, modelNames = name,
                                     initialization = list(hcPairs = hc("EEE", data = S[sub,-i]),
                                                           subset = sub)),
                      silent = TRUE) }
               else
                 { if((allow.EEE) & ncol(S)==2 & sum(is.finite(sBIC$BIC))==0)
                     { try(sBIC <- Mclust(as.matrix(S[,-i]), G = G, modelNames = name,
                                          initialization = list(hcPairs = hcE(S[sub,-i,drop=FALSE]), 
                                                                subset = sub)),
                           silent = TRUE) }
                 }
              rdep <- if(sum(is.finite(sBIC$BIC))==0) NA 
                      else max(sBIC$BIC[is.finite(sBIC$BIC)])
              #
              return(c(regBIC, rdep))
            }
            regBIC <- sapply(out, "[[", 1)
            rdep <- sapply(out, "[[", 2)
            # cindepBIC is the BIC for the clustering model on the other 
            # variables in S and the regression model of the proposed
            #  variable on the other variables in S
            cindepBIC <- regBIC + rdep
            # depBIC is the BIC for the clustering model with all 
            # variables in S
            depBIC <- BICS
            # cdiff is the difference between BIC for the models with
            # variables' being clustering variables versus them being
            # conditionally independent of the clustering
            cdiff <- depBIC - cindepBIC
            # Choose the variable with the smallest difference
            m <- min(cdiff[is.finite(cdiff)])
            arg <- which(cdiff==m,arr.ind=TRUE)[1]
            if(cdiff[arg] <= BIC.diff)
              { # if this difference is less than BIC.diff remove this
                # variable from S and update the clustering model's BICS
                BICS <- rdep[arg] 
                k <- c(colnames(NS),colnames(S)[arg])
                nks <- c(colnames(S)[-arg])
                info <- rbind(info,
                              c(colnames(S)[arg], BICS, cdiff[arg], 
                                "Remove","Accepted")) 
                NS <- cbind(NS,S[,arg])
                S <- as.matrix(S[,-arg])
                colnames(S) <- nks
                colnames(NS) <- k
              }
            else
              { info <- rbind(info,
                           c(colnames(S)[arg], BICS, cdiff[arg], 
                             "Remove","Rejected"))
              }
          }
      }

    if(verbose) print(info)
    # Check if the variables in S have changed or not
    check2 <- colnames(S)
    if(is.null(check2)) # all variables have been removed
      { criterion <- 0 }
    else
      # if they have changed (either added one or removed one or changed one)
      # then continue the algorithm (criterion is 1) otherwise stop 
      # (criterion is 0)
      { if(length(check2) != length(check1))
          { criterion <- 1 }
        else
          { criterion <- if(sum(check1==check2) != length(check1)) 1 else 0 }
      }
  }
      
  if(iter >= itermax) 
    warning("Algorithm stopped because maximum number of iterations was reached")
  
  # List the selected variables and the matrix of steps' information
  info[,2] <- as.numeric(info[,2])
  info[,3] <- as.numeric(info[,3])
  colnames(info) <- c("Variable proposed", "BIC", "BIC difference", 
                      "Type of step", "Decision")
  varnames <- colnames(X)
  subset <- if(is.null(S)) NULL 
            else sapply(colnames(S), function(x) which(x == varnames))

  out <- list(variables = varnames,
              subset = subset,
              # sel.var = S, 
              steps.info = info,
              search = "greedy",
              direction = "forward")
  class(out) <- "clustvarsel"  
  return(out)
}



