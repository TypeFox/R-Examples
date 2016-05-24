ps.estimate.bin <- function(data,           
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
  if ( match.T ){

    weights <- NULL
    weights.str <- NULL

  }

  ## ##################################
  ## crude exp(effect) in original data
  
  if ( !match.T ){ ## stratified

    lr <- glm(resp ~ treat,
              data = data,
              family = family,
              ...)
   
    crude <- list(effect = exp(summary(lr)$coeff[2,1]),
                  se     = summary(lr)$coeff[2,2])

  }else{ ## matched

 
    data.o  <- data[stratum.index == 1,]  ## stratum.index == 1 ==> original data
    resp.o  <- resp[stratum.index == 1]
    treat.o <- treat[stratum.index == 1]

    lr    <- glm(resp.o ~ treat.o, 
                 data = data.o,  
                 family = family,
                 ...)

    crude <- list(effect = exp(summary(lr)$coeff[2,1]),
                  se     = summary(lr)$coeff[2,2])
  }
  
  
  ## ###################################################################
  ## unadjusted analysis ==> stratum-specific parameter will be averaged
  ##                     ==> object$matched.data are used ##############
  
  if ( !match.T ){ ## stratified
    
    or.table <- table(resp, treat, stratum.index)
  
    p0.str <- (or.table[2,,]/colSums(or.table))[1,]
    p1.str <- (or.table[2,,]/colSums(or.table))[2,]

    ad.str <- or.table[1,1,]*or.table[2,2,]
    bc.str <- or.table[1,2,]*or.table[2,1,]

    n.str <- apply(or.table,3,sum)

    or.str <- (or.table[1,1,]*or.table[2,2,])/(or.table[1,2,]*or.table[2,1,])

    weights.str <- as.numeric(table(stratum.index)/dim(data)[1])

    ## response rates estimate
    effect.marg <-
      (sum(weights.str*p1.str, na.rm=TRUE)/(1-sum(weights.str*p1.str, na.rm=TRUE)))/
        (sum(weights.str*p0.str, na.rm=TRUE)/(1-sum(weights.str*p0.str, na.rm=TRUE)))

    se.marg <-
      sqrt(sum(weights.str^2*p1.str*(1-p1.str)/colSums(or.table)[2,], na.rm=TRUE)/
           ((sum(weights.str*p1.str, na.rm=TRUE)*(1-sum(weights.str*p1.str, na.rm=TRUE)))^2) +

           sum(weights.str^2*p0.str*(1-p0.str)/colSums(or.table)[1,], na.rm=TRUE)/
           ((sum(weights.str*p0.str, na.rm=TRUE)*(1-sum(weights.str*p0.str, na.rm=TRUE)))^2))

  
    ## mantel haenszel estimate
    effect.mh <- sum(ad.str/n.str, na.rm=TRUE)/sum(bc.str/n.str, na.rm=TRUE)

    se.mh <-
      sqrt(sum((or.table[1,1,]*or.table[2,2,]/table(stratum.index))*
               (or.table[1,1,]+or.table[2,2,])/table(stratum.index))/
           (2*sum(or.table[1,1,]*or.table[2,2,]/table(stratum.index), na.rm=TRUE)^2) +
           sum((or.table[1,2,]*or.table[2,1,]/table(stratum.index))*
               ((or.table[1,1,]+or.table[2,2,])/table(stratum.index)) +
               (or.table[1,1,]*or.table[2,2,]/table(stratum.index))*
               ((or.table[1,2,]+or.table[2,1,])/table(stratum.index))) /
           (2*sum(or.table[1,1,]*or.table[2,2,]/table(stratum.index), na.rm=TRUE)*
            sum(or.table[1,2,]*or.table[2,1,]/table(stratum.index), na.rm=TRUE)) + 
           sum((or.table[1,2,]*or.table[2,1,] /
                table(stratum.index))*(or.table[1,2,]+or.table[2,1,])/table(stratum.index))/
           (2*sum(or.table[1,2,]*or.table[2,1,]/table(stratum.index), na.rm=TRUE)^2))    
    
    unadj <- list(effect.mh = effect.mh,
                  odds.str  = or.str,
                  se.mh     = se.mh,
                  effect    = effect.marg,
                  se        = se.marg,
                  p1        = sum(weights.str*p1.str, na.rm=TRUE),
                  p0        = sum(weights.str*p0.str, na.rm=TRUE),
                  p1.str    = p1.str,
                  p0.str    = p0.str)

  }else{ ## matched

    data.m     <- data[stratum.index==2,]        ## strata == 2 <==> matched data
    resp.m     <- resp[stratum.index==2]
    treat.m    <- treat[stratum.index==2]
    match.id.m <- match.id[match.id != 0]

    ## accounting for matched data structure
    ## require(lme4, quietly=TRUE)
#    require( "lme4", character.only=TRUE )

    lr.m <- glmer(resp.m ~ treat.m + (1|match.id.m),
                  data = data.m,
                  family = family,
                  ...)
    unadj <- list(effect = exp(as.numeric(fixef(lr.m)[2])),
                  se     = sqrt(vcov(lr.m)[2,2]))

  }

  
  ## #############################
  ## adjusted analysis per stratum

  if ( !is.null(adj) ){

    if (name.resp != names(model.frame(adj.form, data))[1]){

      stop(paste("Argument 'adj' does not include ", name.resp,
                 " as dependent variable.", sep=""))
    }else{      

      if (name.treat != names(model.frame(adj.form, data))[2]){
        stop(paste("Argument 'adj' does not include ", name.treat,
                   " as treatment.", sep=""))
        
      }else{

        if ( !match.T ){## stratified   

          effect.adj.str <-
            se.effect.adj.str <- vector(length=length(levels.stratum.index))
        
          for (i in 1:length(levels.stratum.index)){

            adj.model <- glm(adj.form,
                             data = data[stratum.index == i,], family = family, ...)
            
            if (any(rownames(summary(adj.model)$coeff) == name.treat)){
              
              effect.adj.str[i]     <-
                exp(summary(adj.model)$coeff[rownames(summary(adj.model)$coeff) ==
                                             name.treat,1])
              se.effect.adj.str[i] <-
                summary(adj.model)$coeff[rownames(summary(adj.model)$coeff) ==
                                         name.treat,2]
            }else{
              
              effect.adj.str[i] <- se.effect.adj.str[i] <- NA
              
            }
          }
          
          adj <- list(model      = adj.form,
                      effect.str = effect.adj.str,
                      ## se.str  = se.effect.adj.str,
                      effect     = sum(weights.str*effect.adj.str, na.rm=TRUE),
                      se         = sum(weights.str*se.effect.adj.str, na.rm=TRUE))
          
        }else{ ## matched
          
          data.m     <- data[stratum.index == 2,]   ## strata == 2 <==> matched data 
          adj.form.m <- formula(paste(c(adj.form,
                                        "(1|match.id.m)"),
                                      collapse="+"))
          
          lr.m <- glmer(adj.form.m,
                        data = data.m,
                        family = family,
                        ...)
          
          adj <- list(model  = adj.form.m,
                      effect = exp(as.numeric(fixef(lr.m)[2])),
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

## effect estimates are exp, but se's are on log-scale !!!!!!
