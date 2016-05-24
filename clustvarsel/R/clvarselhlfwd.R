#############################################################################
## Headlong forward search
#############################################################################

clvarselhlfwd <- function(X, G = 1:9, 
                          emModels1 = c("E","V"), 
                          emModels2 = mclust.options("emModelNames"),
                          samp = FALSE, sampsize = 2000, 
                          hcModel = "VVV",
                          allow.EEE = TRUE, forcetwo = TRUE, 
                          BIC.upper = 0, BIC.lower = -10,
                          itermax = 100)
{
  X <- as.matrix(X)
  n <- nrow(X) # number of rows=number of observations
  d <- ncol(X) # number of columns=number of variables
  G <- setdiff(G, 1)

  # If needed, sample the subset of observations for hierarchical clustering
  if(samp) { sub <- sample(1:n, min(sampsize,n), replace = FALSE) }
  else     { sub <- seq.int(1,n) }

  # First Step - selecting single variable and ordering
  maxBIC <- BICdif <- oneBIC <- rep(NA,d)
  for(i in 1:d)
  { 
    xBIC <- NULL
    # Fit the single variable cluster models
    try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1,
                       initialization = list(subset = sub)),
        silent = TRUE)
    if(is.null(xBIC)) 
       try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1),
           silent = TRUE)
    # If we get all NA's from "V" starting hierarchical values use "E"
    if((allow.EEE) & sum(is.finite(xBIC$BIC))==0)
      try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1, 
                         initialization = list(hcPairs = hcE(X[sub,i]),
                                               subset = sub)),
          silent = TRUE)
    # maxBIC is the maximum BIC over all clustering models fit
    if(sum(is.finite(xBIC$BIC))==0) 
       maxBIC[i] <- NA 
    else 
       maxBIC[i] <- max(xBIC$BIC[is.finite(xBIC$BIC)])
    # Fit and get BIC for a single component no-cluster normal model
    try(oneBIC[i] <- Mclust(X[,i], G = 1, modelNames = emModels1,
                            initialization = list(subset = sub))$BIC[1],
        silent = TRUE)
    # Difference between maximum BIC for clustering and BIC for no clustering
    BICdif[i] <- c(maxBIC[i] - oneBIC[i])
  }
  
  # Find the single variable with the biggest difference between 
  # clustering and no clustering
  m <- max(BICdif[is.finite(BICdif)])
  arg <- which(BICdif==m,arr.ind=TRUE)[1]
  # This is our first selected variable/S is the matrix of currently selected
  # clustering variables
  S <- X[,arg,drop=FALSE]
  # BICS is the BIC value for the clustering model with the variable(s) in S
  BICS <- maxBIC[arg]
  temp <- order(BICdif[-arg], decreasing = TRUE)
  # NS is the matrix of currently not selected variables
  NS <- as.matrix(X[,-arg])
  # This orders NS in terms of strongest evidence of univariate clustering 
  # versus no clustering
  NS <- NS[,temp,drop=FALSE]
  # info records the proposed variable, BIC for the S matrix and difference 
  # in BIC for clustering versus no clustering on S, whether it was an 
  # addition step and if it was accepted
  info <- data.frame(Var = colnames(S), BIC = BICS, BICdif = BICdif[arg], 
                     Step = "Add", Decision = "Accepted",
                     stringsAsFactors = FALSE)

  # Second Step - selecting second variable
  depBIC <- cindepBIC <- cdiff <- rep(NA, ncol(NS))
  crit <- -Inf
  i <- 0
  # We only run until we find a variable whose difference in BIC between 
  # being included in the clustering variables versus conditionally 
  # independent of the clustering is greater than BIC.upper
  while( (crit <= BIC.upper) & (i < ncol(NS)) )
  {
    i <- i+1
    
    # Calculate the BIC for the regression of the proposed variable 
    # on the variable in S
    regBIC <- BICreg(y = NS[,i], x = S)

    sBIC <- NULL
    # Fit the cluster model on the two variables
    try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, modelNames = emModels2,
                       initialization = list(hcPairs = hc(hcModel, data = cbind(S,NS[,i])[sub,]), subset = sub)),                       
        silent = TRUE)
    # If we get all NA's from "VVV" starting hierarchical values use "EEE"
    if((allow.EEE) & sum(is.finite(sBIC$BIC))==0)
       try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, modelNames = emModels2,
                          initialization = list(hcPairs = hc("EEE", data = cbind(S,NS[,i])[sub,]), subset = sub)),
           silent = TRUE)
    # depBIC is the BIC for the clustering model with both variables
    if(sum(is.finite(sBIC$BIC))>0) 
      depBIC[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
    # cindepBIC is the BIC for the clustering model on S and the regression 
    # model of the new variable on S
    cindepBIC[i] <- regBIC + BICS
    cdiff[i] <- depBIC[i] - cindepBIC[i]
    if(!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
    crit <- cdiff[i]
  }
  depBIC <- depBIC[1:i]
  cindepBIC <- cindepBIC[1:i]
  cdiff <- cdiff[1:i]

  # i is the index of those variables not selected but whose evidence of 
  # clustering BIC did not fall below BIC.lower or those not looked at yet
  if(cdiff[i] > BIC.upper)
  {
    # i.e. evidence is stronger for including variable in S
    k <- c(colnames(S),colnames(NS)[i])
    S <- cbind(S,NS[,i])
    colnames(S) <- k
    BICS <- depBIC[i]
    info <- rbind(info, c(colnames(NS)[i],BICS,cdiff[i],"Add","Accepted"))
    ns <- s <- NULL
    if(i < ncol(NS))
      ns <- seq(i+1, ncol(NS))
    if(i > 1) 
      s <- seq(i-1)[which(cdiff[-i] > BIC.lower)]
    ind <- c(s,ns)
    if(!is.null(ind))
      { nks <- c(colnames(NS)[ind])
        # NS is the not selected clustering variables whose recently calculated
        # evidence of clustering BIC was higher than BIC.lower or variables not 
        # yet looked at
        NS <- as.matrix(NS[,ind])
        colnames(NS) <- nks
      } 
    else
      { NS <- NULL }
  } 
  else
  {
    if((cdiff[i] < BIC.upper) & (forcetwo))
      {
        # if the evidence is weaker but we're forcing choice of second variable
        m <- max(cdiff[is.finite(cdiff)])
        i <- which(cdiff==m,arr.ind=TRUE)[1]
        k <- c(colnames(S),colnames(NS)[i])
        S <- cbind(S,NS[,i])
        colnames(S) <- k
        BICS <- depBIC[i]
        info <- rbind(info, c(colnames(NS)[i],BICS,cdiff[i],"Add","Accepted"))
        nks <- c(colnames(NS)[-i])
        # NS is the not selected clustering variables whose recently calculated
        # evidence of clustering BIC was higher than BIC.lower
        NS <- as.matrix(NS[,-i])
        temp <- cdiff[-i]
        if(sum(temp > BIC.lower) != 0)
          { NS <- as.matrix(NS[,c(which(temp > BIC.lower))])
            colnames(NS) <- nks[c(which(temp > BIC.lower))]
          } 
        else
          { NS <- NULL } 
      } 
    else
      { m <- max(cdiff[is.finite(cdiff)])
        i <- which(cdiff==m,arr.ind=TRUE)[1]
        info <- rbind(info, c(colnames(NS)[i],BICS,cdiff[i],"Add","Rejected"))
      }
  }

  criterion <- 1
  iter <- 0
  while((criterion == 1) & (iter < itermax))
  {
    iter <- iter+1
    check1 <- colnames(S)

    # Addition step
    # For the special case where we have removed all the clustering variables/S 
    # is empty [LS: is really needed??]
    if((NCOL(NS) != 0 & !is.null(ncol(NS))) & 
       (ncol(S) == 0) || (is.null(ncol(S))) )
      { # cat("\n really needed ??\n")
        depBIC <- 0
        DepBIC <- NULL
        crit <- -10
        cdiff <- 0
        Cdiff <- NULL
        oneBIC <- rep(NA,d)
        i <- 0
        crit <- -10
        while((crit <= BIC.upper) & (i < ncol(NS)))
        {
          xBIC <- NULL
          i <- i+1
          # Fit the cluster models 
          try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1,
                             initialization = list(subset = sub)),
              silent = TRUE)
          # If we get all NA's from "V" starting hierarchical values use "E"
          if((allow.EEE) & sum(is.finite(xBIC$BIC))==0)
             try(xBIC <- Mclust(X[,i], G = G, modelNames = emModels1,
                                initialization = list(hcPairs = hcE(X[sub,i]),
                                                      subset = sub)),
                 silent = TRUE)
          # depBIC is the maximum BIC over all clustering models fit
          if(sum(is.finite(xBIC$BIC)) == 0) 
            depBIC <- NA 
          else 
            depBIC <- max(xBIC$BIC[is.finite(xBIC$BIC)])
          DepBIC <- c(DepBIC,depBIC)
          # Fit and get BIC for a single component no-cluster normal model
          try(oneBIC <- Mclust(X[,i], G = 1, modelNames = "V")$BIC[1],
              silent = TRUE)
          # Difference between maximum BIC for clustering and BIC for no clustering
          cdiff <- c(depBIC - oneBIC)
          if(!is.finite(cdiff)) cdiff <- BIC.upper
          Cdiff <- c(Cdiff,cdiff)
          crit <- cdiff
        }
        if(cdiff > BIC.upper)
        { # ie. evidence is stronger for including variable in S
          k <- c(colnames(NS)[i])
          S <- as.matrix(NS[,i])
          colnames(S) <- k
          BICS <- depBIC
          info <- rbind(info, c(colnames(NS)[i],BICS,cdiff,"Add","Accepted"))
          ns <- s <- NULL
          # i is the index of those variables not selected but whose evidence of
          # clustering BIC did not fall below BIC.lower or those not looked at yet
          if(i < ncol(NS)) ns <- seq(i+1, ncol(NS))
          if(i > 1) s <- seq(i-1)[which(Cdiff[-i] > BIC.lower)]
          ind <- c(s,ns)
          if(!is.null(ind))
            { nks <- c(colnames(NS)[ind])
              # NS is the not selected clustering variables whose recently
              # calculated evidence of clustering BIC was higher than BIC.lower
              # or variables not yet looked at
              NS <- as.matrix(NS[,ind])
              colnames(NS)<-nks
            } 
          else
            { NS <- NULL }
        } 
        else
        {
          m <- max(Cdiff[is.finite(Cdiff)])
          i <- which(Cdiff==m,arr.ind=TRUE)[1]
          info <- rbind(info, c(colnames(NS)[i],BICS,Cdiff[i],"Add","Rejected"))
          ind <- seq(ncol(NS))[which(Cdiff > BIC.lower)]
          if(!is.null(ind))
            { k <- colnames(NS)[ind]
              # Exclude variables in NS whose evidence of clustering in this step 
              # was lower than BIC.lower
              NS <- as.matrix(NS[,ind])
              colnames(NS) <- k
            } 
          else
            { NS <- NULL }
        }
      } 
      # [LS: is really needed?? -- END]
    else
      {
        # Addition Step in general (for all cases except when S is empty)
        if((NCOL(NS) != 0) & !is.null(ncol(NS)))
          {
            depBIC <- cindepBIC <- cdiff <- rep(NA, ncol(NS))
            crit <- -Inf
            i <- 0
            # We only run until we find a variable whose difference in BIC 
            # between being included in the clustering variables versus 
            # conditionally independent of the clustering is greater than BIC.upper
            while(crit <= BIC.upper & i < ncol(NS))
            { 
              sBIC <- NULL
              i <- i+1
              # Calculate the BIC for the regression of the proposed variable 
              # on the variable(s) in S
              regBIC <- BICreg(y = NS[,i], x = S)
          
              # Fit the cluster model on the S variables with the proposed variable 
              try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, modelNames = emModels2,
                                 initialization = list(hcPairs = hc(hcModel, data = cbind(S,NS[,i])[sub,]), subset = sub)),
                  silent = TRUE)
              # If we get all NA's from "VVV" starting hierarchical values use "EEE"
              if((allow.EEE) & (sum(is.finite(sBIC$BIC))==0))
                try(sBIC <- Mclust(cbind(S,NS[,i]), G = G, modelNames = emModels2,
                                   initialization = list(hcPairs = hc("EEE", data = cbind(S,NS[,i])[sub,]),
                                                        subset = sub)),
                    silent = TRUE)
              # depBIC is the BIC for the clustering model with both S and proposed
              # variable
              if(sum(is.finite(sBIC$BIC))>0) 
                depBIC[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
              # cindepBIC is the BIC for the clustering model on S and the
              # regression model of the new variable on S
              cindepBIC[i] <- regBIC + BICS
              cdiff[i] <- depBIC[i] - cindepBIC[i]
              if(!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
              crit <- cdiff[i]
            }
            depBIC <- depBIC[1:i]
            cindepBIC <- cindepBIC[1:i]
            cdiff <- cdiff[1:i]
            if(cdiff[i] > BIC.upper)
            {
              # i.e. evidence is stronger for including variable in S
              k <- c(colnames(S),colnames(NS)[i])
              nks <- c(colnames(NS)[-i])
              S <- cbind(S,NS[,i])
              colnames(S) <- k
              BICS <- depBIC[i]
              info <- rbind(info,
                            c(colnames(NS)[i],BICS,cdiff[i],"Add","Accepted"))
              ns <- s <- NULL
              # Exclude variables in NS whose evidence of clustering in this step 
              # was lower than BIC.lower
              if(i < ncol(NS)) ns <- seq(i+1,ncol(NS))
              if(i > 1) s <- seq(i-1)[which(cdiff[-i] > BIC.lower)]
              ind <- c(s,ns)
              if(!is.null(ind))
                { nks <- colnames(NS)[ind]
                  NS <- as.matrix(NS[,ind])
                  colnames(NS) <- nks
                } 
              else
                { NS <- NULL }
            } 
            else
            {
              m <- max(cdiff[is.finite(cdiff)])
              i <- which(cdiff==m,arr.ind=TRUE)[1]
              info <- rbind(info,
                            c(colnames(NS)[i],BICS,cdiff[i],"Add","Rejected"))
              ind <- seq(1,ncol(NS))[which(cdiff > BIC.lower)]
              if(!is.null(ind))
                { k <- colnames(NS)[ind]
                  # Exclude variables in NS whose evidence of clustering in this 
                  # step was lower than BIC.lower
                  NS <- as.matrix(NS[,ind])
                  colnames(NS) <- k
                } 
              else
                { NS <- NULL }
            }
          }
      }

    # Removal Step for the special case where S contains only a single variable
    if(ncol(S) == 1)
      { cdiff <- 0
        oneBIC <- NA
        try(oneBIC <- Mclust(S, G = 1, modelNames = "V", 
                             initialization = list(subset = sub))$BIC[1], 
            silent = TRUE)
        # Difference between maximum BIC for clustering and BIC for no clustering
        cdiff <- BICS - oneBIC
        if(is.na(cdiff)) cdiff <- BIC.upper
        # check if difference is negative
        if(cdiff <= BIC.upper)
          {
            # if negative remove the variable from S and set the BIC for 
            # the model to NA
            BICS <- NA
            info <- rbind(info, c(colnames(S),BICS,cdiff,"Remove","Accepted"))
            # Only return variable to NS if difference is greater than BIC.lower
            if(cdiff > BIC.lower)
              { k <- c(colnames(NS),colnames(S))
                NS <- cbind(NS,S)
                colnames(NS) <- k
                S <- NULL
              } 
            else
              { S <- NULL }
          } 
        else
          { info <- rbind(info, c(colnames(S),BICS,cdiff,"Remove","Rejected")) }
      } 
    else
      {
        # Removal step in general (for all cases except when S is a single 
        # variable or empty)
        if(ncol(S) >= 2)
          { depBIC <- BICS
            cindepBIC <- rdep <- cdiff <- rep(NA, ncol(S))
            crit <- Inf
            i <- 0
            # Check if the data is at least 3 dimensional
            name <- if(ncol(S) > 2) emModels2 else emModels1
            # We only run until we find a variable whose difference in BIC between 
            # being included in the clustering variables versus conditionally 
            # independent of the clustering is lower than BIC.upper
            while(crit > BIC.upper & (i<ncol(S)))
            {
              i <- i+1
              # Calculate the BIC for the regression of the proposed variable 
              # from S on the other variable(s) in S
              regBIC <- BICreg(y = S[,i], x = S[,-i])
              # Fit the cluster model on the S variables without the 
              # proposed variable 
              sBIC <- NULL
              try(sBIC <- Mclust(S[,-i], G = G, modelNames = name,
                                 initialization = 
                                 list(hcPairs = hc(hcModel, data = S[sub,-i,drop=FALSE]), subset = sub)),
                  silent = TRUE)
              # If we get all NA's from "VVV" starting hierarchical values use "EEE"
              if(allow.EEE & (ncol(S) >= 3) & sum(is.finite(sBIC$BIC))==0)
                { try(sBIC <- Mclust(S[,-i], G = G, modelNames = name,
                                     initialization = list(hcPairs = hc("EEE", data = S[sub,-i]),
                                                           subset = sub)),
                      silent = TRUE) }
              else
                { if((allow.EEE) & (ncol(S)==2) & sum(is.finite(sBIC$BIC))==0)
                    { try(sBIC <- Mclust(as.matrix(S[,-i]), G = G, modelNames = name,
                                         initialization = list(hcPairs = hcE(S[sub,-i,drop=FALSE]), 
                                                               subset = sub)),
                          silent = TRUE) }
                } 
              if(sum(is.finite(sBIC$BIC))>0) 
                rdep[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
              # cindepBIC is the BIC for the clustering model on the other variables 
              # in S and the regression model of the proposed variable on the other
              # variables in S
              cindepBIC[i] <- regBIC + rdep[i]
              cdiff[i] <- depBIC - cindepBIC[i]
              if(!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
              crit <- cdiff[i]
            }
            if((cdiff[i] < BIC.upper) & (cdiff[i] > BIC.lower))
              { # i.e. evidence is stronger for excluding variable from S but 
                # still including it in NS
                BICS <- rdep[i]
                info <- rbind(info, c(colnames(S)[i], BICS, cdiff[i],
                              "Remove","Accepted"))
                k <- c(colnames(NS),colnames(S)[i])
                nk <- colnames(S)[-i]
                NS <- cbind(NS,S[,i])
                S <- as.matrix(S[,-i])
                colnames(NS) <- k
                colnames(S) <- nk
              }
            else
              { if(cdiff[i] < BIC.lower)
                  { # exclude variable entirely
                    BICS <- rdep[i]
                    info <- rbind(info, c(colnames(S)[i], BICS, cdiff[i],
                                  "Remove","Accepted"))
                    nk <- colnames(S)[-i]
                    S <- as.matrix(S[,-i])
                    colnames(S) <- nk
                  } 
                else
                  { m <- min(cdiff[is.finite(cdiff)])
                    i <- which(cdiff==m,arr.ind=TRUE)[1]
                    info <- rbind(info, c(colnames(S)[i], BICS, cdiff[i],
                                  "Remove","Rejected"))
                  }
              }
          }
        }
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
  info[,2] = as.numeric(info[,2])
  info[,3] = as.numeric(info[,3])
  colnames(info) <- c("Variable proposed", "BIC", "BIC difference", 
                      "Type of step", "Decision")
  # colnames(info) <- c("Variable proposed","BIC of new clustering variables set","BIC difference","Type of step","Decision")
  varnames <- colnames(X)
  subset <- sapply(colnames(S), function(x) which(x == varnames))

  out <- list(variables = varnames,
              subset = subset,
              steps.info = info,
              search = "headlong",
              direction = "forward")
  return(out)

}
