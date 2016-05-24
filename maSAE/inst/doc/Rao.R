setOldClass(c('lme'))
setClass(Class = "sae"
         , representation = representation(
             lmeObj = "lme"
             , domain.df = "data.frame"
             )
         , validity = function(object){
           if(object@lmeObj$method != 'REML')
             return("lme method wasn't REML, can't handle this")
           if(! all(c(attr((getGroups(object@lmeObj)), 'label')
                      , dimnames(attr(terms(formula(object@lmeObj))
				      , "factors"))[[2]])
                    %in% names(object@domain.df)
                    )
              ) return('Names in domain.df do not match names in lmeObj')
           unitName <- attributes(getGroups(object@lmeObj))$label
           if(! all(eval(parse(text = paste(sep = ''
                                 , 'object@domain.df$'
                                 , unitName)))
                    %in% unique(getGroups(object@lmeObj))
                    )
              ) return('Unknown unit in domain.df')
         }
         )
setGeneric(name = "sae"
           , def = function(object, ...){standardGeneric("sae")}
           )


setMethod(f = "sae"
          , signature(object = "sae")
          , function(object){
            lmeObj  <- object@lmeObj
            domain.df <- object@domain.df

            eval(parse(text=paste(sep = ''
                         , 'domain.df$'
                         , setdiff(dimnames(attr(terms(formula(lmeObj))
						 , "factors"))[[1]]
                                   , dimnames(attr(terms(formula(lmeObj))
						   , "factors"))[[2]])
                         , '<- 1'))
                 )

            unitName <- attributes(getGroups(lmeObj))$label
            units <- unique(getGroups(lmeObj))
	    unitseq <- seq(along = units); names(unitseq) <- units
            n <- tapply(getGroups(lmeObj), getGroups(lmeObj), length)
	    ## the estimated Variances, see p.137, last paragraph
            varV <- as.numeric(VarCorr(lmeObj)[1,1])
            varE <- as.numeric(VarCorr(lmeObj)[2,1])
            ## Vector of the mean Residual per Unit
            mUnitRes <- tapply(resid(lmeObj, level = 0)
                               , names(resid(lmeObj, level = 0))
			       , mean)[levels(units)]
            ## list of design matrices, multipurpose
            lX <- lapply(X = unitseq
                         , FUN = function(i, lmeObj, unitName){
                           unit <- units[i]
                           X <- model.matrix(formula(lmeObj)
                                             , subset(lmeObj$data
                                                      , eval(parse(text
								   = unitName))
						      == unit))
                           return(X)
                         }
                         , lmeObj
                         , unitName
                         )
            ## list of R7.2.4, multipurpose
            lR7.2.4 <- lapply(X = unitseq
                              , FUN = function(i, n, varV, varE){
                                unit <- units[i]
                                n.i <- n[unit]
                                R7.2.4 <- varV / (varV + varE / n.i)
                                return(R7.2.4)
                              }
                              , n
                              , varV
                              , varE
                              )
            ## list of R7.2.7 to build the sum in R7.2.12
            lR7.2.7 <- lapply(X = unitseq
                              , FUN = function(i, n, lX, lR7.2.4, varV, varE){
                                unit <- units[i]
                                X <- lX[[unit]]
                                n.i <- n[unit]
                                R7.2.4 <- lR7.2.4[[unit]]
                                R7.2.2 <- 1 / varE *
                                  (diag(n.i)  - R7.2.4 / n.i * rep(1, n.i)
				   %*% t(rep(1, n.i)))
                                R7.2.7 <- t(X) %*% R7.2.2 %*% X
                                return(R7.2.7)
                              }
                              , n
                              , lX
                              , lR7.2.4
                              , varV
                              , varE
                              )
            sumAi <- Reduce('+',lR7.2.7)
            lR7.2.30 <- lapply(X = unitseq
                               , FUN = function(i, n, varV, varE){
                                 unit <- units[i]
                                 n.i <- n[unit]
                                 R7.2.30 <- varE + varV * n.i
                                 return(R7.2.30)
                               }
                               , n
                               , varV
                               , varE
                               )
            R7.2.27 <- 0.5 * sum(n^2 * unlist(lR7.2.30)^-2)
            R7.2.28 <- 0.5 * sum((n - 1) * sqrt(varE)^-4 + unlist(lR7.2.30)^-2 )
            R7.2.29 <- 0.5 * sum(n * unlist(lR7.2.30)^-2)
            informationMatrix <- matrix(data = c(R7.2.27, R7.2.29, R7.2.29, R7.2.28)
                                        , ncol = 2)
            asymtoticCovarianceMatrix <- solve(informationMatrix)
            R7.2.23 <- varE^2 * asymtoticCovarianceMatrix[1, 1] +
              varV^2 * asymtoticCovarianceMatrix[2, 2] -
                2 * varE * varV * asymtoticCovarianceMatrix[1, 2]
            units <- eval(parse(text = paste(sep = ''
					     , 'domain.df$'
					     , unitName)))
            unitseq <- seq(along = units)
            names(unitseq) <- units
            foo <- lapply(##  over observations in domain.df
                          X = unitseq
                          , FUN = function(i, n, lX, lR.7.2.4, varV
                              , varE, sumAi, R7.2.23){
                            unit <- units[i]
                            X <- lX[[unit]]
                            n.i <- n[unit]
                            R7.2.4 <- lR7.2.4[[unit]]
                            R7.2.11 <- R7.2.4 * varE / n.i
                            Xbar <- model.matrix(formula(lmeObj)
                                                 , subset(domain.df
							  , eval(parse(text =
								       unitName)) == unit
							  )
						 )
                            xbar <- colMeans(X)
                            tmp <- as.numeric(Xbar) - R7.2.4 * xbar
                            R7.2.12 <- t(tmp) %*% solve(sumAi) %*% tmp
                            R7.2.30 <- lR7.2.30[units[unit]]
                            R7.2.22 <- n.i^-2 * (varV + varE / n.i)^-3 * R7.2.23
                            R7.2.32 <- n.i^-2 * (varV + varE
						 / n.i)^-4 * R7.2.23 * mUnitRes[unit]^2
                            R7.2.33 <- R7.2.11 + R7.2.12 + 2 * R7.2.32
                            R7.2.34 <- R7.2.11 + R7.2.12 + R7.2.22 + R7.2.32
                            return(c(mse1 = R7.2.33
                                     , mse2 = R7.2.34
				     ))
                          }
                          , n
                          , lX
                          , lR7.2.4
                          , varV
                          , varE
                          , sumAi
                          , R7.2.23
                          )
	    ## transform list to transposed data.frame
            foo <- as.data.frame(t(as.data.frame(foo)))
            return(cbind(eblup = as.numeric(predict(lmeObj
						    , newdata = domain.df))
	    , foo))
          }
          )

