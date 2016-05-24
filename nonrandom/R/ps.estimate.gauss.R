ps.estimate.gauss <- function(data,
                              resp,
                              name.resp,
                              treat,
                              name.treat,
                              stratum.index,
                              adj.form,
                              adj,
                              name.adj,
                              weights,
                              levels.stratum.index,
                              family,
                              match.T,
                              match.id,
                              ...)
{
  ## #####################
  ## Set match.T indicator
  if (match.T){

    weights <- NULL
    weights.str <- NULL

  }
  
  ## #################################
  ## crude effect in original data set
  
  if (!match.T){ ## stratified

    lr     <- glm(resp~treat,
                  data=data,
                  family=family,
                  ...)

    crude  <- list(effect = summary(lr)$coeff[2,1],
                   se     = summary(lr)$coeff[2,2])

  }else{ ## matched

    data.o  <- data[stratum.index == 1,] ## strata == 1 ==> original data
    resp.o  <- resp[stratum.index == 1]
    treat.o <- treat[stratum.index == 1]

    lr    <- glm(resp.o ~ treat.o, 
                 data = data.o,  
                 family = family,
                 ...)

    crude <- list(effect = summary(lr)$coeff[2,1],
                  se     = (summary(lr)$coeff[2,2])) 
  }


  ## ###################################################################
  ## unadjusted analysis ==> stratum-specific parameter will be averaged
  ##                     ==> object$matched.data are used ##############
  
  if (!match.T){ ## stratified  

    effect.str <- q.str   <- vector(length=length(levels.stratum.index))
    se0.str    <- se1.str <- vector(length=length(levels.stratum.index))
    n0.str     <- n1.str  <- vector(length=length(levels.stratum.index))
    
    for (i in 1:length(levels.stratum.index)){

      if (any(c(length(resp[treat == max(treat, na.rm = TRUE) & stratum.index == i]),
                length(resp[treat == min(treat, na.rm = TRUE) & stratum.index == i])) == 0)){
        
        effect.str[i] <- q.str[i] <- NA
        se0.str[i] <- se1.str[i] <- NA
        n0.str[i] <- n1.str[i] <- NA
        
      }else{
        
        effect.str[i] <-
          mean(resp[treat == max(treat, na.rm = TRUE) & stratum.index == i], na.rm=TRUE) -
            mean(resp[treat == min(treat, na.rm = TRUE) & stratum.index == i], na.rm=TRUE)
        
        q.str[i] <-
          (dim(data[stratum.index == i & treat == min(treat, na.rm = TRUE),])[1] +
           dim(data[stratum.index == i & treat == max(treat, na.rm = TRUE),])[1]) /
             (dim(data[stratum.index == i & treat == min(treat, na.rm = TRUE),])[1] *
              dim(data[stratum.index == i & treat == max(treat, na.rm = TRUE),])[1])
        
        se0.str[i] <-
          sd(resp[stratum.index == i & treat == min(treat, na.rm = TRUE)], na.rm=TRUE)
        se1.str[i] <-
          sd(resp[stratum.index == i & treat == max(treat, na.rm = TRUE)], na.rm=TRUE)
        
        n0.str[i] <-
          dim(data[stratum.index == i & treat == min(treat, na.rm = TRUE),])[1]
        n1.str[i] <-
          dim(data[stratum.index == i & treat == max(treat, na.rm = TRUE),])[1]
      }
    }
    
    if (weights == "opt"){

      weights.str <- (1/sum(1/q.str, na.rm=TRUE))*(1/q.str)

    }else{

      weights.str <- as.numeric(table(stratum.index)/dim(data)[1])

    }
    
    unadj <-
      list(effect.str = effect.str,
           effect     = sum(weights.str*effect.str, na.rm=TRUE),
           se         = sqrt(sum(weights.str^2*q.str, na.rm=TRUE)*var(resp))#,
           ##se.unequal = sqrt(sum(weights.str^2*(se0.str^2/n0.str + se1.str^2/n1.str),na.rm=TRUE))
           )

  }else{ ## matched

    data.m     <- data[stratum.index == 2,] ## strata == 2 <==> matched data
    resp.m     <- resp[stratum.index == 2]
    treat.m    <- treat[stratum.index == 2]
    match.id.m <- match.id[match.id != 0]

    ## accounting for matched data structure
    ## require(lme4, quietly=TRUE)

    ## CHG 28.3.14: loaded by package depends.
#    require( "lme4", character.only=TRUE )

    lr.m <- glmer(resp.m ~ treat.m + (1|match.id.m),
                  data=data.m,
                  family=family,
                  ...)
    
    unadj <- list(effect = as.numeric(fixef(lr.m)[2]),
                  se     = sqrt(vcov(lr.m)[2,2]))
  }

  
  ## #############################
  ## adjusted analysis per stratum
  
  if (!is.null(adj)){   

    if (name.resp != names(model.frame(adj.form, data))[1]){

      stop(paste("Argument 'adj' does not include ", name.resp,
                 " as dependent variable.", sep=""))
    }else{

      if (name.treat != names(model.frame(adj.form, data))[2]){
        stop(paste("Argument 'adj' does not include ", name.treat,
                   " as treatment.", sep=""))

      }else{
        
        if (!match.T){ ## stratified
          
          effect.adj.str <-
            se.effect.adj.str <- vector(length=length(levels.stratum.index))
          
          for (i in 1:length(levels.stratum.index)){
            
            adj.model <- glm(adj.form,
                             data=data[stratum.index == i,], family=family, ...)
            
            if (any(rownames(summary(adj.model)$coeff) == name.treat)){
              
              effect.adj.str[i] <-
                summary(adj.model)$coeff[rownames(summary(adj.model)$coeff) ==
                                         name.treat,1]
              se.effect.adj.str[i] <-
                summary(adj.model)$coeff[rownames(summary(adj.model)$coeff) ==
                                         name.treat,2]
              
            }else{
              
              effect.adj.str[i] <- se.effect.adj.str[i] <- NA
              
            }
          }
          
          adj <- list(model      = adj.form,
                      effect.str = effect.adj.str,
                      ##se.str     = se.effect.adj.str,
                      effect     = sum(weights.str*effect.adj.str, na.rm=TRUE),
                      se         = sum(weights.str*se.effect.adj.str, na.rm=TRUE))
          
        }else{ ## matched

          data.m     <- data[stratum.index == 2,]   ## strata == 2 <==> matched data
          adj.form.m <- formula(paste(c(adj.form,
                                        "(1|match.id.m)"),
                                      collapse="+"))
          
          lr.m <- glmer(adj.form.m,
                        data=data.m,
                        family=family,
                        ...)

          adj <- list(model  = adj.form.m,
                      effect = as.numeric(fixef(lr.m)[2]),
                      se     = sqrt(vcov(lr.m)[2,2]))         
        }
      }
    }
  }else{
    adj <- "No adjustment"
  }

    ps.estimation <- list(crude       = crude,
                          unadj       = unadj,
                          adj         = adj,
                          weights     = weights,
                          weights.str = weights.str)

  return(ps.estimation)
  
}

