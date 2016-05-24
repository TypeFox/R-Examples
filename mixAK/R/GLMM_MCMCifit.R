##
##  PURPOSE:   Initial (RE)ML fits for GLMM_MCMC function,
##             removal of missing values from design/response matrices,
##             variables derived from design matrices
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    04/08/2009
##
##  FUNCTIONS:  GLMM_MCMCifit
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCifit
## *************************************************************
##
GLMM_MCMCifit <- function(do.init, na.complete,
                          y, dist, id, time, x, z, random.intercept,
                          xempty, zempty, Rc, Rd, p, p_fi, q, q_ri, lalpha, dimb)
{
##### Variables in the resulting object:
#####             Y                 a list of length R with observations really used in fitting process
#####                               (after removal of missing values)
#####             ID                a list of length R with id's corresponding to Y
#####             time              original time after removal of missing values (if na.complete)  
#####             x                 a list resulting from the original argument x after removal of
#####                               observations with some missing information
#####                               additionaly, intercept column is added if fixed intercept included in the model
#####             z                 a list resulting from the original argument z after removal of
#####                               observations with some missing information
#####                               additionaly, intercept column is added if random intercept included in the model
#####
##### ------------------------------------------------------------------------------------------------------------------------------------
#####
#####             I                 number of clusters in the original data (before removing NA's)
#####             n                 list of length R, each component is a vector or length I (may contain zeros
#####                               if some cluster disappears for particular response due to NA's)
#####
##### ------------------------------------------------------------------------------------------------------------------------------------
#####  
#####             Cn                vectorized n
#####             sumCn             sum(Cn) = total number of observations  
#####             Cy_c              vector with continuous response to be passed to C++, equal to 0 if there is no continuous response
#####             Cy_d              vector with discrete response to be passed to C++, equal to 0 if there is no discrete response
#####             CX                vector containing X matrices (without ones for possible intercept)
#####                               to be passed to C++, equal to 0 if there are no X matrices
#####             CZ                vector containing Z matrices (without ones for possible intercept)
#####                               to be passed to C++, equal to 0 if there are no Z matrices
#####             CXtX              vector containing lower triangles of X[s]'X[s] for each response to be passed to C++,
#####                               where X[s] contains also column of ones if there is a fixed intercept in the model,
#####                               equal to 0 if there are no alpha coefficients in the model
#####                        AS OF 21/10/2009, CXtX IS NO MORE COMPUTED
#####                        AND NOT INCLUDED IN THE RESULTING OBJECT  
#####                       
#####             CZitZi            vector containig lower triangles of matrices Z[i,s]'Z[i,s] for each cluster and each response
#####                               to be passed to C++,
#####                               where Z[i,s] contains also column of ones if there is a random intercept in the model,
#####                               Remark:   for cluster without any observations of response s, Z[i,s] is considered
#####                                         to be a matrix of zeros
#####                               => for CZitZi, matrices Z[i,s]'Z[i,s] are for all i of the same dimension 
#####                               ordering: cluster 1: ZitZi matrices for s=1,...,R, cluster 2: ZitZi matrices for s=1,...,R, etc.,
#####                               equal to 0 if there are no random effects in the model
#####                        AS OF 20/10/2009, CZitZi IS NO MORE COMPUTED
#####                        AND NOT INCLUDED IN THE RESULTING OBJECT  
#####
##### ------------------------------------------------------------------------------------------------------------------------------------
#####
#####             iintcpt           data.frame(Est, SE) with estimated intercepts and their SE, R rows,
#####                               row equal to (0, 0) if there is no fixed intercept for particular response  
#####             ifixef            a list of length R, each component is equal to 0 if there are no fixed effects for particular response,
#####                               and is equal to data.frame(Est, SE) if there are fixed effects
#####             isigma            vector of length R, equal to 0 for discrete response, equal to estimated residual
#####                               standard deviation for continuous response  
#####             iEranef           a list of length R, each component is equal to 0 if there are no random effects for particular response,
#####                               and is equal to data.frame(Est, SE) with estimated means of the random effects
#####                               and their std. errors if there are random effects
#####             iSDranef          a list of length R, each component is equal to 0 if there are no random effects for particular response,
#####                               and is equal to a vector with estimated standard deviations of the random effects if there are random effects
#####             ib                a list of length R, each component is equal to 0 if there are no random effects for particular response,
#####                               and a matrix with EB estimates of random effects shifted by their estimated mean if there are random effects
#####              
#####             is.intcpt         logical vector of length R
#####             is.fixef          logical vector of length R
#####             is.ranef          logical vector of length R  
#####             is.sigma          logical vector of length R
#####
##### ------------------------------------------------------------------------------------------------------------------------------------
#####  
#####             ibMat          matrix with initial values of random effects (EB estimates from (RE)ML fits)
#####             ibMat2         matrix with alternative initial values of random effects  
#####             iEranefVec     vector with estimated means of random effects
#####             iSEranefVec    vector with standard errors of estimated means of random effects
#####             iSDranefVec    vector with estimated standard deviations of random effects
#####             ialpha         vector with initial values of alpha's (including fixed intercepts)
#####             ialpha2        vector with alternative initial values of alpha's (including fixed intercepts)  
#####             iSEalpha       vector with standard errors of estimated values of fixed effects  
#####  
##### ------------------------------------------------------------------------------------------------------------------------------------
#####  
##### Other (potentially useful) variables created here:
#####             TAB               table(id)
#####
##### ------------------------------------------------------------------------------------------------------------------------------------    
  R <- Rc + Rd

  ##### Handle NA's if na.complete
  ##### --------------------------------------------------
  if (na.complete){
    dta <- data.frame(id=id)
    dta <- cbind(dta, y)
    for (s in 1:R){
      if (!xempty[s]) dta <- cbind(dta, x[[s]])    
      if (!zempty[s]) dta <- cbind(dta, z[[s]])
    }
    is.NA <- as.logical(apply(is.na(dta), 1, sum))

    if (ncol(y) == 1) y <- data.frame(y1=y[!is.NA,])    ## this avoids that y is changed to a numeric vector
    else              y <- y[!is.NA,]
    
    id <- id[!is.NA]
    time <- time[!is.NA]
    for (s in 1:R){
      if (!xempty[s]){
        if (p[s] == 1) x[[s]] <- data.frame(x1=x[[s]][!is.NA,])
        else           x[[s]] <- x[[s]][!is.NA,]
      }  
      if (!zempty[s]){
        if (q[s] == 1) z[[s]] <- data.frame(z1=z[[s]][!is.NA,])
        else           z[[s]] <- z[[s]][!is.NA,]
      }  
    }
  }  
  
  ##### Remove missing values, create vectors and matrices which will be passed to C++
  ##### Do initial ((RE)ML) fits (this will also check X and Z matrices for linear independence of their columns)
  ##### ---------------------------------------------------------------------------------------------------------
  I <- length(unique(id))
  TAB <- table(id)
  cumTAB <- c(1, cumsum(TAB) + 1)
  n <- list()
  Y <- list()
  ID <- list()
  
  if (Rc) Cy_c <- numeric(0) else Cy_c <- 0
  if (Rd) Cy_d <- numeric(0) else Cy_d <- 0
  Cn     <- numeric(0)
  CX     <- numeric(0)
  CZ     <- numeric(0)
  ##CXtX   <- numeric(0)        ### REMOVED ON 21/10/2009
  ##tmpZitZi <- list()          ### will contain ZitZi matrices separately for each response
                                ### will be reshuffeled afterwards to get main order by cluster
                                ### REMOVED ON 20/10/2009
  LT_q_ri <- (q_ri * (1 + q_ri)) / 2

  
  ##### Objects to keep values from initial fits
  ##### --------------------------------------------------
  iintcpt  <- data.frame(Est=numeric(R), SE=numeric(R))
  ifixef   <- list()
  isigma   <- numeric(R)  
  iEranef  <- list()
  iSDranef <- list()
  ib       <- list()
  rownames(iintcpt) <- names(isigma) <- colnames(y)
  
  is.intcpt <- rep(FALSE, R)
  is.fixef  <- rep(FALSE, R)
  is.sigma  <- rep(FALSE, R)
  is.ranef  <- rep(FALSE, R)
  names(is.intcpt) <- names(is.fixef) <- names(is.sigma) <- names(is.ranef) <- colnames(y)

  for (s in 1:R){
    dta <- data.frame(id=id, y=y[,s])
    if (!xempty[s]) dta <- cbind(dta, x[[s]])    
    if (!zempty[s]) dta <- cbind(dta, z[[s]])
    is.NA <- as.logical(apply(is.na(dta), 1, sum))
    dta <- dta[!is.NA,]

    Y[[s]]  <- dta[,"y"]
    ID[[s]] <- dta[,"id"]
    
    dTAB <- table(dta[,"id"])
    n[[s]] <- rep(0, length(TAB))
    names(n[[s]]) <- names(TAB)      
    n[[s]][names(dTAB)] <- dTAB

    Cn <- c(Cn, n[[s]])

    ##tmpZitZi[[s]] <- numeric(0)          ### REMOVED ON 20/10/2009
    
    ##### Computation of basis for initial values
    ##### ---------------------------------------------------------------------
    if (dist[s] %in% c("gaussian"))                             Cy_c <- c(Cy_c, dta[,"y"])
    else if (dist[s] %in% c("binomial(logit)", "poisson(log)")) Cy_d <- c(Cy_d, dta[,"y"])

    if (zempty[s] & !random.intercept[s]){             ## no random effects
      iEranef[[s]]  <- 0         
      iSDranef[[s]] <- 0
      ib[[s]]       <- 0
      z[[s]] <- "empty"
      
      if (xempty[s]){
        ##CXtX <- c(CXtX, nrow(dta))                     ## X'X = 1'1 = nrow(X)
                                                         ### REMOVED ON 21/10/2009
        if (do.init){
          FORM <- formula("y ~ 1")
          if (dist[s] %in% c("gaussian"))               ifit <- lm(FORM, data=dta)
          else if (dist[s] %in% c("binomial(logit)"))   ifit <- glm(FORM, family=binomial(link = "logit"), data=dta)
               else if (dist[s] %in% c("poisson(log)")) ifit <- glm(FORM, family=poisson(link = "log"), data=dta)
          sifit <- summary(ifit)

          ifixef[[s]] <- 0          
        }          
        
        x[[s]] <- matrix(1, nrow=nrow(dta), ncol=1)
      }  
      else{
        CX <- c(CX, t(dta[, -(1:2)]))
        ##### ----- CODE REMOVED ON 21/10/2009 ----- #####
        ##tmpX <- as.matrix(cbind(1, dta[, -(1:2)]))
        ##tmpXtX <- t(tmpX) %*% tmpX
        ##CXtX <- c(CXtX, tmpXtX[lower.tri(tmpXtX, diag=TRUE)])          
        ##### ----- END OF CODE REMOVED ON 21/10/2009 ----- #####
        
        if (do.init){        
          FORM <- formula(paste("y ~", paste(colnames(x[[s]]), collapse=" + ")))
          if (dist[s] %in% c("gaussian"))               ifit <- lm(FORM, data=dta)
          else if (dist[s] %in% c("binomial(logit)"))   ifit <- glm(FORM, family=binomial(link = "logit"), data=dta)
               else if (dist[s] %in% c("poisson(log)")) ifit <- glm(FORM, family=poisson(link = "log"), data=dta)          
          sifit <- summary(ifit)

          ifixef[[s]] <- data.frame(Est=coef(ifit)[-1], SE=sifit[["coefficients"]][-1, "Std. Error"])
          is.fixef[s] <- TRUE
        }  
        
        x[[s]] <- as.matrix(cbind(1, dta[, 3:(2+p[s])]))
      }
      if (do.init){
        iintcpt[s, "Est"] <- coef(ifit)["(Intercept)"]
        iintcpt[s, "SE"]  <- sifit[["coefficients"]][1, "Std. Error"]
        is.intcpt[s] <- TRUE

        if (dist[s] %in% c("gaussian")){
          isigma[s] <- sifit[["sigma"]]
          is.sigma[s] <- TRUE
        }
      }  
    }else{                                           ## there are some random effects
      #cat("Initial values for response ", s, " (", dist[s], "):\n", sep="")
      if (xempty[s]){
        if (!zempty[s] & !random.intercept[s]){
          #cat("xempty & !zempty & !random.intercept\n")
          ##CXtX <- c(CXtX, nrow(dta))                     ## X'X = 1'1 = nrow(X)
                                                           ### REMOVED ON 21/10/2009
          
          CZ <- c(CZ, t(dta[,-(1:2)]))
          ##### ----- CODE REMOVED ON 20/10/2009 ----- #####
          ##tmpZ <- as.matrix(z[[s]])
          ##tmpZ[is.NA,] <- 0
          ##for (i in 1:I){
          ##  tmpZi <- matrix(tmpZ[cumTAB[i]:(cumTAB[i+1]-1),], ncol=ncol(tmpZ))
          ##  tmpZtZ <- t(tmpZi) %*% tmpZi
          ##  tmpZitZi[[s]] <- c(tmpZitZi[[s]], tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])
          ##}  
          ##### ----- END OF CODE REMOVED ON 20/10/2009 ----- #####
          
          if (do.init){          
            FORM <- formula(paste("y ~ 1 +", paste(colnames(z[[s]]), collapse=" + "), " + (-1 +", paste(colnames(z[[s]]), collapse=" + "), " | id)"))
            OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)
            if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
            else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                 else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)
            options(OPT)                                   ### added on 20140521 (version 3.6)

            iintcpt[s, "Est"] <- lme4::fixef(ifit)["(Intercept)"]
            iintcpt[s, "SE"]  <- as.numeric(sqrt(vcov(ifit)[1, 1]))                                           ### lme4::vcov
            is.intcpt[s] <- TRUE        
         
            iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit)[-1], SE=sqrt(diag(as.matrix(vcov(ifit)))[-1]))   ### lme4::vcov
            is.ranef[s] <- TRUE
          }
            
          x[[s]] <- matrix(1, nrow=nrow(dta), ncol=1)
          z[[s]] <- as.matrix(dta[, (3+p[s]):(2+p[s]+q[s])])
        }else{
          if (zempty[s] & random.intercept[s]){
            #cat("xempty & zempty & random.intercept\n")            
            ##tmpZitZi[[s]] <- n[[s]]                    ## Z'Z = 1'1, where 1 = vector of length equal to number of obs. within cluster
                                                         ### REMOVED ON 20/10/2009

            if (do.init){            
              FORM <- formula("y ~ 1 + (1 | id)")
              OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)
              if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
              else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                   else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)
              options(OPT)                                   ### added on 20140521 (version 3.6)              

              iintcpt[s, "Est"] <- 0
              iintcpt[s, "SE"]  <- 0

              iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit), SE=sqrt(diag(as.matrix(vcov(ifit)))))   ### lme4::vcov
              is.ranef[s] <- TRUE
            }
              
            x[[s]] <- "empty"
            z[[s]] <- matrix(1, nrow=nrow(dta), ncol=1)
          }else{
            if (!zempty[s] & random.intercept[s]){
              #cat("xempty & !zempty & random.intercept\n")                          
              CZ <- c(CZ, t(dta[,-(1:2)]))
              ##### ----- CODE REMOVED ON 20/10/2009 ----- #####              
              ##tmpZ <- as.matrix(cbind(1, z[[s]]))
              ##tmpZ[is.NA,] <- 0
              ##for (i in 1:I){
              ##  tmpZi <- matrix(tmpZ[cumTAB[i]:(cumTAB[i+1]-1),], ncol=ncol(tmpZ))
              ##  tmpZtZ <- t(tmpZi) %*% tmpZi
              ##  tmpZitZi[[s]] <- c(tmpZitZi[[s]], tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])   ### REMOVED ON 20/10/2009
              ##}  
              ##### ----- END OF CODE REMOVED ON 20/10/2009 ----- #####              
              
              if (do.init){              
                FORM <- formula(paste("y ~ 1 +", paste(colnames(z[[s]]), collapse=" + "), " + (1 +", paste(colnames(z[[s]]), collapse=" + "), " | id)"))
                OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)                
                if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
                else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                     else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)
                options(OPT)                                   ### added on 20140521 (version 3.6)
                
                iintcpt[s, "Est"] <- 0
                iintcpt[s, "SE"]  <- 0

                iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit), SE=sqrt(diag(as.matrix(vcov(ifit)))))   ### lme4::vcov
                is.ranef[s] <- TRUE
              }
                
              x[[s]] <- "empty"
              z[[s]] <- as.matrix(cbind(1, dta[, (3+p[s]):(2+p[s]+q[s])]))
            }
          }  
        }
        if (do.init) ifixef[[s]] <- 0
      }                   ## end of if (xempty[s])
      else{               ## else of if (xempty[s])
        CX <- c(CX, t(dta[, 3:(2+ncol(x[[s]]))]))          
        if (!zempty[s] & !random.intercept[s]){
          #cat("!xempty & !zempty & !random.intercept\n")
          ##### ----- CODE REMOVED ON 21/10/2009 ----- #####
          ##tmpX <- as.matrix(cbind(1, dta[, 3:(2+ncol(x[[s]]))]))
          ##tmpXtX <- t(tmpX) %*% tmpX
          ##CXtX <- c(CXtX, tmpXtX[lower.tri(tmpXtX, diag=TRUE)])
          ##### ----- END OF CODE REMOVED ON 21/10/2009 ----- #####
          
          CZ <- c(CZ, t(dta[, (2+ncol(x[[s]])+1):(2+ncol(x[[s]])+ncol(z[[s]]))]))
          ##### ----- CODE REMOVED ON 20/10/2009 ----- #####                        
          ##tmpZ <- as.matrix(z[[s]])
          ##tmpZ[is.NA,] <- 0
          ##for (i in 1:I){
          ##  tmpZi <- matrix(tmpZ[cumTAB[i]:(cumTAB[i+1]-1),], ncol=ncol(tmpZ))
          ##  tmpZtZ <- t(tmpZi) %*% tmpZi
          ##  tmpZitZi[[s]] <- c(tmpZitZi[[s]], tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])
          ##}            
          ##### ----- END OF CODE REMOVED ON 20/10/2009 ----- #####              
          
          if (do.init){          
            FORM <- formula(paste("y ~ 1 +", paste(colnames(x[[s]]), collapse=" + "), " + ", paste(colnames(z[[s]]), collapse=" + "), " + (-1 +", paste(colnames(z[[s]]), collapse=" + "), " | id)"))
            OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)           
            if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
            else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                 else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)
            options(OPT)                                   ### added on 20140521 (version 3.6)            

            iintcpt[s, "Est"] <- lme4::fixef(ifit)["(Intercept)"]
            iintcpt[s, "SE"]  <- as.numeric(sqrt(vcov(ifit)[1, 1]))                                                                                  ### lme4::vcov
            is.intcpt[s] <- TRUE
      
            ifixef[[s]] <- data.frame(Est=lme4::fixef(ifit)[2:(1+ncol(x[[s]]))], SE=sqrt(diag(as.matrix(vcov(ifit)))[2:(1+ncol(x[[s]]))]))           ### lme4::vcov
            is.fixef[s] <- TRUE
      
            iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit)[-(1:(1+ncol(x[[s]])))], SE=sqrt(diag(as.matrix(vcov(ifit)))[-(1:(1+ncol(x[[s]])))]))    ### lme4::vcov
            is.ranef[s] <- TRUE
          }
            
          x[[s]] <- as.matrix(cbind(1, dta[, 3:(2+p[s])]))
          z[[s]] <- as.matrix(dta[, (3+p[s]):(2+p[s]+q[s])])
        }else{
          ##### ----- CODE REMOVED ON 21/10/2009 ----- #####
          ##tmpX <- as.matrix(dta[, 3:(2+ncol(x[[s]]))])
          ##tmpXtX <- t(tmpX) %*% tmpX
          ##CXtX <- c(CXtX, tmpXtX[lower.tri(tmpXtX, diag=TRUE)])
          ##### ----- END OF CODE REMOVED ON 21/10/2009 ----- #####
          
          if (zempty[s] & random.intercept[s]){
            #cat("!xempty & zempty & random.intercept\n")            
            ##tmpZitZi[[s]] <- n[[s]]                    ## Z'Z = 1'1, where 1 = vector of length equal to number of obs. within cluster
                                                         ### REMOVED ON 20/10/2009
            if (do.init){            
              FORM <- formula(paste("y ~ 1 +", paste(colnames(x[[s]]), collapse=" + "), " + (1 | id)"))
              OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)
              if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
              else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                   else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)
              options(OPT)                                   ### added on 20140521 (version 3.6)
              
              iintcpt[s, "Est"] <- 0
              iintcpt[s, "SE"]  <- 0

              ifixef[[s]] <- data.frame(Est=lme4::fixef(ifit)[-1], SE=sqrt(diag(as.matrix(vcov(ifit)))[-1]))           ### lme4::vcov
              is.fixef[s] <- TRUE
      
              iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit)[1], SE=sqrt(diag(as.matrix(vcov(ifit)))[1]))            ### lme4::vcov
              is.ranef[s] <- TRUE
            }
              
            x[[s]] <- as.matrix(dta[, 3:(2+p[s])])
            z[[s]] <- matrix(1, nrow=nrow(dta), ncol=1)            
          }else{
            if (!zempty[s] & random.intercept[s]){
              #cat("!xempty & !zempty & random.intercept\n")              
              CZ <- c(CZ, t(dta[, (2+ncol(x[[s]])+1):(2+ncol(x[[s]])+ncol(z[[s]]))]))
              ##### ----- CODE REMOVED ON 20/10/2009 ----- #####                            
              ##tmpZ <- as.matrix(cbind(1, z[[s]]))
              ##tmpZ[is.NA,] <- 0
              ##for (i in 1:I){
              ##  tmpZi <- matrix(tmpZ[cumTAB[i]:(cumTAB[i+1]-1),], ncol=ncol(tmpZ))
              ##  tmpZtZ <- t(tmpZi) %*% tmpZi
              ##  tmpZitZi[[s]] <- c(tmpZitZi[[s]], tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])
              ##}  
              ##### ----- END OF CODE REMOVED ON 20/10/2009 ----- #####              
              
              if (do.init){
                FORM <- formula(paste("y ~ 1 +", paste(colnames(x[[s]]), collapse=" + "), " + ", paste(colnames(z[[s]]), collapse=" + "), " + (1 +", paste(colnames(z[[s]]), collapse=" + "), " | id)"))
                OPT <- options(warn = -1)                      ### added on 20140521 (version 3.6)
                if (dist[s] %in% c("gaussian"))               ifit <- lme4::lmer(FORM, data=dta)
                else if (dist[s] %in% c("binomial(logit)"))   ifit <- lme4::glmer(FORM, family=binomial(link = "logit"), data=dta)
                     else if (dist[s] %in% c("poisson(log)")) ifit <- lme4::glmer(FORM, family=poisson(link = "log"), data=dta)                            
                options(OPT)                                   ### added on 20140521 (version 3.6)
                
                iintcpt[s, "Est"] <- 0
                iintcpt[s, "SE"]  <- 0

                ifixef[[s]] <- data.frame(Est=lme4::fixef(ifit)[2:(1+ncol(x[[s]]))], SE=sqrt(diag(as.matrix(vcov(ifit)))[2:(1+ncol(x[[s]]))]))     ### lme4::vcov
                is.fixef[s] <- TRUE
      
                iRAND <- c(1, (2+ncol(x[[s]])):(1+ncol(x[[s]])+ncol(z[[s]])))
                iEranef[[s]] <- data.frame(Est=lme4::fixef(ifit)[iRAND], SE=sqrt(diag(as.matrix(vcov(ifit)))[iRAND]))                              ### lme4::vcov
                is.ranef[s] <- TRUE
              }  
              
              x[[s]] <- as.matrix(dta[, 3:(2+p[s])])
              z[[s]] <- as.matrix(cbind(1, dta[, (3+p[s]):(2+p[s]+q[s])]))              
            }
          }  
        }  
      }                   ## end of else of if (xempty[s])

      if (do.init){      
        if (dist[s] %in% c("gaussian")){
          isigma[s] <- attr(lme4::VarCorr(ifit), "sc")
          is.sigma[s] <- TRUE
        }  

        iSDranef[[s]] <- attr(lme4::VarCorr(ifit)[["id"]], "stddev")

        if (length(iSDranef[[s]]) == 1) bb <- data.frame(b1=rnorm(I, iEranef[[s]][,"Est"], iSDranef[[s]]))
        else{
          Vb <- diag(iSDranef[[s]]^2)
          bb <- as.data.frame(rMVN(I, iEranef[[s]][,"Est"], Sigma=Vb)[["x"]])
        }
        rownames(bb) <- names(TAB)
        bb[names(dTAB),] <- lme4::ranef(ifit)[["id"]] + matrix(rep(iEranef[[s]][,"Est"], nrow(lme4::ranef(ifit)[["id"]])), ncol=ncol(lme4::ranef(ifit)[["id"]]), byrow=TRUE)
        ib[[s]] <- bb
      }  
    }                    ## end of else:  there are some random effects

    ##if (length(tmpZitZi[[s]]) != I * LT_q_ri[s]) stop(paste("BUG in the function (strange length of tmpZitZi[[", s, "]]), contact AK!", sep=""))
    ### REMOVED ON 20/10/2009
  }       ### end of for (s in 1:R)

  sumCn <- sum(Cn)

  ##### ----- CODE REMOVED ON 20/10/2009 ----- #####
  #CZitZi <- numeric(0)
  #if (sum(q_ri)){
  #  for (i in 1:I){
  #    for (s in 1:R){
  #      if (q_ri[s]){
  #        CZitZi <- c(CZitZi, tmpZitZi[[s]][((i - 1)*LT_q_ri[s] + 1):(i*LT_q_ri[s])])
  #      }  
  #    }  
  #  }      
  #}  
  #if (!length(CZitZi)) CZitZi <- 0
  #
  #if (length(CZitZi) != I*sum(LT_q_ri)) stop("BUG in the function (strange length of CZitZi), contact AK!")
  #if (any(abs(CZitZi) >= Inf)) stop("numerically not stable (ZitZi contains +-Inf)")
  #if (any(is.na(CZitZi))) stop("numerically not stable (ZitZi contains NaN/NA)")  
  ##### ----- END OF CODE REMOVED ON 20/10/2009 ----- #####
  
  names(n) <- colnames(y)
  if (!length(CX))     CX <- 0
  if (!length(CZ))     CZ <- 0

  ##### ----- CODE REMOVED ON 21/10/2009 ----- #####
  #if (!length(CXtX))   CXtX <- 0
  #
  #if (length(CXtX) != sum((p_fi * (1 + p_fi)) / 2)) stop("BUG in the function (strange length of CXtX), contact AK!")
  #if (any(abs(CXtX) >= Inf)) stop("numerically not stable (XtX contains +-Inf)")
  #if (any(is.na(CXtX))) stop("numerically not stable (XtX contains NaN/NA)")
  ##### ----- CODE REMOVED ON 21/10/2009 ----- #####
  
  RET <- list(Y      = Y,
              ID     = ID,
              time   = time,   
              x      = x,
              z      = z,
              I      = I,
              n      = n,
              Cn     = Cn,
              sumCn  = sumCn,
              Cy_c   = Cy_c,
              Cy_d   = Cy_d,
              CX     = CX,
              CZ     = CZ)
              #CXtX   = CXtX)      ## REMOVED ON 21/10/2009
              #CZitZi = CZitZi)    ## REMOVED ON 20/10/2009

  if (do.init){    
    names(ifixef) <- names(iEranef) <- names(iSDranef) <- names(ib) <- colnames(y)
    
    RET$iintcpt    <- iintcpt
    RET$ifixef     <- ifixef
    RET$isigma     <- isigma
    RET$iEranef    <- iEranef
    RET$iSDranef   <- iSDranef
    RET$ib         <- ib
    RET$is.intcpt  <- is.intcpt
    RET$is.fixef   <- is.fixef
    RET$is.ranef   <- is.ranef
    RET$is.sigma   <- is.sigma

    ########## ========== Unlist moments of random effects from ML fits                          ========= ##########
    ########## ========== Create a matrix with possible initial values of b                      ========= ##########  
    ########## ========== Create a vector of possible initial values of alpha                    ========= ##########
    ########## =========================================================================================== ##########
    if (dimb){
      iEranefVec <- iSEranefVec <- iSDranefVec <- numeric(0)
      ibempty <- TRUE
      for (s in 1:R){
        if (is.ranef[s]){
          iEranefVec <- c(iEranefVec, iEranef[[s]][,"Est"])
          iSEranefVec <- c(iSEranefVec, iEranef[[s]][,"SE"])
          iSDranefVec <- c(iSDranefVec, iSDranef[[s]])
          if (ibempty){
            ibMat <- ib[[s]]
            ibMat2 <- ib[[s]] + matrix(rnorm(nrow(ib[[s]])*ncol(ib[[s]]), mean=0, sd=rep(0.1*iSDranef[[s]], each=nrow(ib[[s]]))), ncol=ncol(ib[[s]]))
            ibempty <- FALSE
          }else{
            ibMat <- cbind(ibMat, ib[[s]])
            ibMat2 <- cbind(ibMat2, ib[[s]] + matrix(rnorm(nrow(ib[[s]])*ncol(ib[[s]]), mean=0, sd=rep(0.1*iSDranef[[s]], each=nrow(ib[[s]]))), ncol=ncol(ib[[s]])))
          }  
        }  
      }
      rownames(ibMat) <- rownames(ibMat2) <- paste(1:nrow(ibMat))
      colnames(ibMat) <- colnames(ibMat2) <- paste("b", 1:dimb, sep="")
    }else{
      iEranefVec <- iSEranefVec <- iSDranefVec <- ibMat <- ibMat2 <- 0
    }        
    
    if (lalpha){
      ialpha <- ialpha2 <- iSEalpha <- numeric(0)
      for (s in 1:R){
        if (is.intcpt[s]){
          ialpha  <- c(ialpha, iintcpt[s, "Est"])
          ialpha2 <- c(ialpha2, rnorm(1, mean=iintcpt[s, "Est"], sd=1*iintcpt[s, "SE"]))
          iSEalpha <- c(iSEalpha, iintcpt[s, "SE"])
        }  
        if (is.fixef[s]){
          ialpha  <- c(ialpha, ifixef[[s]][, "Est"])
          ialpha2 <- c(ialpha2, rnorm(nrow(ifixef[[s]]), mean=ifixef[[s]][, "Est"], sd=1*ifixef[[s]][, "SE"]))
          iSEalpha <- c(iSEalpha, ifixef[[s]][, "SE"])
        }  
      }  
    }else{
      ialpha <- ialpha2 <- iSEalpha <- 0
    }  

    RET$ibMat       <- ibMat
    RET$ibMat2      <- ibMat2
    RET$iEranefVec  <- iEranefVec
    RET$iSEranefVec <- iSEranefVec
    RET$iSDranefVec <- iSDranefVec
    RET$ialpha      <- ialpha
    RET$ialpha2     <- ialpha2
    RET$iSEalpha    <- iSEalpha
  }              

  return(RET)  
}  
