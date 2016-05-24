  bayes.regress <- function (data.values=NULL, 
                           beta.prior=list("flat"), 
                           sigmasq.prior=list("inverse.gamma", 1.0, 1.0, 1.0),
                           Tsamp.out=1000, zero.intercept=FALSE)

#   This function computes the samples for linear regression coefficients and the variance for normally distributed noise
#  
#   data.values=list(xtx = matrix kxk, 
#                    xty = vector of k values, 
#                    yty = scalar, 
#                    numsamp.data = 1000,
#                    xtx.inv (optional) = kxk matrix) 
#  
#   k = dim(xtx)[1] is the number of unknown regression coefficients (   which equals to the number of
#                                                                        predictors when zero.intercept=TRUE ).
#  
#   Possible beta priors with parameters in the order of numeration 
#   if user does not specifies the names of the list variables: 
#  
#   1) beta.prior = list(type="flat") 
#
#
#   2) beta.prior=list(type="mvnorm.known", 
#                      mean.mu   = matrix(k x 1),
#                      cov.C     = matrix(k x k),
#                      prec.Cinv = matrix(k x k) = optional, it is used in place of cov.C if provided)
#
#     -> match.list( list(type = NULL, mean.mu=NULL, cov.C = NULL, prec.Cinv = NULL ), beta.prior)
#
#   3) beta.prior=list("mvnorm.unknown",
#                       mu.hyper.mean.eta        = matrix(k x 1), # expectation of mean
#                       mu.hyper.prec.Dinv       = matrix(k x k), # inverse of the variation of mean
#                       Cinv.hyper.df.lambda     = value,         # Wishart degree of freedom
#                       Cinv.hyper.invscale.Vinv = matrix(k x k), # Wishart SPD matrix
#                       mu.init     = matrix(k x 1),              # initial value for mu vector
#                       Cinv.init   = matrix(k x k)))             # initial value for C^-1 matrix
# 
#      -> match.list(list( type = NULL,
#                          mu.hyper.mean.eta = NULL,
#                          mu.hyper.prec.Dinv = NULL,
#                          Cinv.hyper.df.lambda      = NULL,
#                          Cinv.hyper.invscale.Vinv = NULL,
#                          mu.init     = NULL,
#                          Cinv.init   = NULL
#                        ), beta.prior)
#
#
#    Possible sigmasq priors:
#
#    sigmasq.prior=list( type="inverse.gamma",   
#                        inverse.gamma.a = value, 
#                        inverse.gamma.b = value, 
#                        sigmasq.init = value)
  
#    -> match.list(list( type = NULL, 
#                        inverse.gamma.a = NULL, 
#                        inverse.gamma.b = NULL, 
#                        sigmasq.init = NULL), sigmasq.prior)
#
#    sigmasq.prior=list( type="sigmasq.inverse", sigmasq.init = value)
#  
#    -> match.list(list( type = NULL, sigmasq.init = NULL), sigmasq.prior)

  {
        
  #check data.values list
    
  if(!all(class(data.values)=="list",length(data.values)>2)) 
    stop("data.values argument must be a list with at least four elements. See help for details.")
  
  template.data.attr <- c("xtx", "xty", "yty", "numsamp.data", "xtx.inv", "")
  
  template.data <- list(xtx=NULL, xty=NULL, yty=NULL, numsamp.data=NULL, xtx.inv = NULL)
  
  data.attr <- attributes(data.values)[['names']]
  
  # check if there are no repeated fields (otherwise stop):
  
  if ( length(data.attr) > 0 )
  {
    i.dup <- anyDuplicated(data.attr[which(data.attr!="")])
    
    if (i.dup > 0)
      stop(paste('formal argument "', 
                 data.attr[i.dup], 
                 '" matched by multiple actual arguments.', sep=''))    
  }

  # check if fields with attributes in user provided priors coincide 
  # with default ones (that is, that there are no incorrect names)
  
  if (length(data.attr)>0)
  {    
    unused.attr <- setdiff(data.attr, template.data.attr)
    if (length(unused.attr)>0) 
      stop(paste('unused argument in data.values:',unused.attr,"\n  "))
  }              

  if (length(data.values)>length(template.data))
  {
    stop(paste('data.values must contain maximum of ', length(template.data),' parameters'))
  }
  
  # match/transform data.values to a template form:
  
  data.values <- match.list(template.data, data.values )
  
  # check that the data.values contains numerical values
  
  flag <- TRUE
  
  for(i in 1:length(data.values))
    if  ((!is.null(data.values[[i]])))
      flag <- (flag && (is.numeric(data.values[[i]])))  
  
  if (flag==FALSE)
    stop(paste('data.values elements must be numerical.\n'))  
    
  
  #check if xtx matrix is supplied. 
  
  if(!is.null(data.values[['xtx']]))
  {
    if(all(is.matrix(data.values[['xtx']]))&isSymmetric(data.values[['xtx']])) 
    {
      xtx <- data.values[['xtx']] 
    }      
    else stop("data.values argument \"xtx\" is not propery formatted");    
  }
  else
    stop("data.values argument \"xtx\" is not propery formatted");
    
  #xtx matrix is provided. Thus, we check the rest of the data.values argument
  
  matsize<-dim(xtx)
  
  k<-matsize[1]
  
  #check if xty matrix is supplied. note ytx is a transpose of xty 
    
  if(!is.null(data.values[['xty']]))
  {
    if (all(is.matrix(data.values[['xty']]),dim(data.values[['xty']])==c(k,1))||all(is.numeric(data.values[['xty']]),length(data.values[['xty']])==k))
    {
      xty <- data.values[['xty']]
      ytx <- t(xty)
    }
    else stop("data.values argument \"xty\" is not properly formatted")
  }
  else stop("data.values argument \"xty\" is not propery formatted")  
  
    
  #check if yty matrix is supplied, since it will most likely come out of a dot product of two matrix variables, we expect a matrix for it
  
  if(!is.null(data.values[['yty']]))
  {
    if(all(is.matrix(data.values[['yty']]),dim(data.values[['yty']])==c(1,1))||all(is.numeric(data.values[['yty']]),length(data.values[['yty']])==1)) 
    {
      yty <- data.values[['yty']]      
    }
    else stop("data.values argument \"yty\" is not properly formatted")
  }
  else stop("data.values argument \"yty\" is not propery formatted");
  
  
  # check numsamp.data entry in data.values
  
  if(!is.null(data.values[['numsamp.data']])) 
  {
    if ((length(data.values[['numsamp.data']])==1) && (is.numeric(data.values[['numsamp.data']])))
    {
      num.samp <- data.values[['numsamp.data']]  
    }
    else stop("data.values argument \"numsamp.data\" is not properly formatted")            
  }    
  else
  {
    warning("parameter numsamp.data size is not specified. 10^3 is used instead.")
    num.samp=1000    
  }
  
  # HERE REDUCE MATRIX XTX if zero.itercept = TRUE since matrices xtx and ytx have been extended
  #  
  if (zero.intercept==TRUE)# extended matrices xtx and xty are given, but someone wants intercept to be zero
  {
    if(k>1) # We assume that the matrices xtx and xty were extended, thus there should be at least two columns
    {
      xtx <- xtx[2:k,2:k]
      xty <- xty[2:k]
      k   <- k-1
    }
    else
      stop("Dimensions of the matrix xtx in data.values are too small for zero.intercept analysis")
  }
  
  # ***************************************************************************************************************************************
  # At this point zero.intersept is known and k has correct value, representing the number of unknown (and desired) regression coefficients
  # ***************************************************************************************************************************************
  
  
  # check beta.prior and sigmasq.prior lists:
  
  beta.prior.name <- c("flat",
                       "mvnorm.known",
                       "mvnorm.unknown")
  
  sigmasq.prior.name <- c("inverse.gamma",
                          "sigmasq.inverse")
  
  # Get attributes of user provided priors:
  
  beta.prior.attr    <- attributes(beta.prior)[['names']]
  
  sigmasq.prior.attr <- attributes(sigmasq.prior)[['names']]
  
  # Check if there are no repeated fields (otherwise stop):
  
  if ( length(beta.prior.attr) > 0 )
  {
    i.dup <- anyDuplicated(beta.prior.attr[which(beta.prior.attr!="")])
    
    if ( i.dup > 0 )
      stop(paste('formal argument "', 
                 beta.prior.attr[i.dup], 
                 '" matched by multiple actual arguments.',sep=''))    
  }
  
  if (length(sigmasq.prior.attr)>0)
  {
    i.dup <- anyDuplicated(sigmasq.prior.attr[which(sigmasq.prior.attr!="")])
    
    if ( i.dup > 0)
      stop(paste('formal argument "', 
                 sigmasq.prior.attr[i.dup], 
                 '" matched by multiple actual arguments.',sep=''))    
  }        
    
  # Check if "beta.prior[['type']]" exists and it is defined correctly:
  
  type <- ""
  
  if((is.null(beta.prior))||(class(beta.prior)!="list"))
    stop("beta.prior argument is not a list or missing")
  
  if (length(beta.prior)==0)
    stop("beta.prior list is empty")

  if (length(beta.prior.attr) == 0)
    type <- beta.prior[[1]][1]
  else
  {
    i.type  <- which(beta.prior.attr == "type")  
    
    if( length(i.type)>1 ) 
      stop(paste('formal argument "type" in beta.prior matched by multiple actual arguments.'))
    
    if (length(i.type)==1)
    {
      type <- beta.prior[i.type[1]][[1]]    
    }
    else
    {    
      i.empty <- which( beta.prior.attr == "" )
      
      if (length(i.empty)>0) type <- beta.prior[i.empty[1]][[1]]
      else type <- ""
    }      
  }

  if (!(type %in% beta.prior.name))  
    stop('"type" in beta.prior is either not set or has an unexpected value')
  

# Check if "sigmasq.prior[['type']]" exists and it is defined correctly:

  sqtype <- ""
  
  if((is.null(sigmasq.prior))||(class(sigmasq.prior)!="list"))
    stop("sigmasq.prior argument is not a list or missing")
  
  if (length(sigmasq.prior)==0)
    stop("sigmasq.prior list is empty")
  
  if (length(sigmasq.prior.attr) == 0)
    sqtype <- sigmasq.prior[[1]][1]
  else
  {
    i.type  <- which(sigmasq.prior.attr == "type")
    
    if( length(i.type)>1 ) 
      stop(paste('formal argument "type" is sigmasq.prior matched by multiple actual arguments.'))
    
    if (length(i.type)==1)
    {
      sqtype <- sigmasq.prior[i.type[1]][[1]]
    }
    else
    {    
      i.empty <- which( sigmasq.prior.attr == "" )
      
      if (length(i.empty)>0) sqtype <- sigmasq.prior[i.empty[1]][[1]]
      else sqtype <- ""
    }      
  }
  
  if (!(sqtype %in% sigmasq.prior.name))  
    stop('"type" in sigmasq.prior is either not set or has an unexpected value')

  # Now type and sqtype contain type of beta and sigmasq priors.

  
  # Assign the lists of priors 
  
  template.prior  <- list(  flat = list(type=NULL),
                            mvnorm.known = list(type = NULL, 
                                                mean.mu=NULL, 
                                                cov.C=NULL, 
                                                prec.Cinv=NULL),
                            mvnorm.unknown=list(type = NULL,
                                                mu.hyper.mean.eta        = NULL,
                                                mu.hyper.prec.Dinv       = NULL,
                                                Cinv.hyper.df.lambda     = NULL,
                                                Cinv.hyper.invscale.Vinv = NULL,
                                                mu.init     = NULL,                                 
                                                Cinv.init   = NULL),                                 
                            inverse.gamma = list(type = NULL, 
                                                 inverse.gamma.a=NULL, 
                                                 inverse.gamma.b=NULL, 
                                                 sigmasq.init = NULL),
                            sigmasq.inverse = list(type = NULL, 
                                                   sigmasq.init = NULL)
  )
  
  # Assign the lists of attributes for each prior ( which includes "" attribute):
  
  template.attr <- list( flat=c("type",""),
                         mvnorm.known = c("type", "mean.mu", "cov.C", "prec.Cinv",""),
                         mvnorm.unknown = c("type", 
                                            "mu.hyper.mean.eta", 
                                            "mu.hyper.prec.Dinv", 
                                            "Cinv.hyper.df.lambda",
                                            "Cinv.hyper.invscale.Vinv",
                                            "mu.init", 
                                            "Cinv.init",""),
                         inverse.gamma = c("type","inverse.gamma.a","inverse.gamma.b","sigmasq.init", ""),
                         sigmasq.inverse = c("type","sigmasq.init","" )
                        )

  # check if fields with attributes in user provided priors coincide 
  # with default ones (that is, that there are no incorrect names)
  
  if (length(beta.prior.attr)>0)
  {    
    unused.attr <- setdiff(beta.prior.attr, template.attr[[type]])
    if (length(unused.attr)>0) 
      stop(paste('unused argument in beta.prior:',unused.attr,"\n  "))
  }              
  
  if (length(sigmasq.prior.attr)>0)
  {
    unused.attr <- setdiff(sigmasq.prior.attr, template.attr[[sqtype]])
    if (length(unused.attr)>0) 
      stop(paste('unused argument in sigmasq.prior:',unused.attr,"\n  "))
  }  

  if (length(beta.prior)>length(template.prior[[type]]))
  {
    stop(paste('beta.prior of type "', type,'" must contain maximum of ', length(template.prior[[type]]),' parameters'))
  }

  if (length(sigmasq.prior)>length(template.prior[[sqtype]]))
  {
    stop(paste('sigmasq.prior of type "', type,'" must contain maximum of ', length(template.prior[[sqtype]]),' parameters'))
  }

  # match/transform proprs to a template form:

  beta.prior    <- match.list(template.prior[[type]],   beta.prior   )
      
  sigmasq.prior <- match.list(template.prior[[sqtype]], sigmasq.prior)
  

# here we check that prior's fields (except "type") are numerical
  
  beta.prior.attr    <- attributes(beta.prior)[['names']]
  
  sigmasq.prior.attr <- attributes(sigmasq.prior)[['names']]
  
  flag <- TRUE
  
  for(i in 1:length(beta.prior))
    if  ((beta.prior.attr[i] != "type") && (!is.null(beta.prior[[i]])))
      flag <- (flag && (is.numeric(beta.prior[[i]])))  
  
  if (flag==FALSE) 
    stop(paste('beta.prior elements except \"type\" must be numerical.\n'))
  
  for(i in 1:length(sigmasq.prior))
    if  ((sigmasq.prior.attr[i] != "type") && (!is.null(sigmasq.prior[[i]])))
      flag <- (flag && (is.numeric(sigmasq.prior[[i]])))  
  
  if (flag==FALSE) 
    stop(paste('sigmasq.prior elements except \"type\" must be numerical.\n'))


  # Now, when everything is in the template form and numerical, we can work:
  
  result<-switch(type,
                 # beta prior distribution is uniform
            flat={
              
              # check if xtx inverse is supplied
              
              if((!is.null(data.values[['xtx.inv']]))&&all(dim(data.values[['xtx.inv']])==dim(xtx))) 
                xtx.inv<-data.values[['xtx.inv']]
              else 
                xtx.inv<-chol2inv(chol(xtx))
              
              switch(sqtype,  
                     inverse.gamma=# sigma squared distribution is Gamma(a,b). we will need two parameters a and b, starting value for sigma squared and preferably the inverse of xtx
                       {
                                             
                         # set default values for inverse gamma prior parameters
                         
                         if (!is.null(sigmasq.prior[["inverse.gamma.a"]]))
                           inv.gamma.a <- as.numeric(sigmasq.prior[["inverse.gamma.a"]])
                         else inv.gamma.a <- 1.0

                         if (!is.null(sigmasq.prior[["inverse.gamma.b"]])) 
                           inv.gamma.b <- as.numeric(sigmasq.prior[["inverse.gamma.b"]])
                         else inv.gamma.b <- 1.0
                                                                          
                         if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                           sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                         else sigmasq.init <- 1.0
                         
                         
                         bayesregressB1S1(xtx,xtx.inv,xty,yty,
                                          numsamp.data = num.samp,
                                          inv.gamma.a = inv.gamma.a, 
                                          inv.gamma.b = inv.gamma.b, 
                                          sigmasq.init = sigmasq.init,                                    
                                          Tsamp.out = Tsamp.out)
                        },
                     sigmasq.inverse=
                        {
                          
                          if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                            sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                          else sigmasq.init=1.0
                          
                          bayesregressB1S2(xtx,xtx.inv,xty,yty,
                                           numsamp.data = num.samp,
                                           sigmasq.init = sigmasq.init,
                                           Tsamp.out = Tsamp.out)
                        },
                     stop("Prior distribution for sigma squared is not properly formatted")
              )
            },
            
            mvnorm.known={
              
              if (!is.null(beta.prior[['mean.mu']]))
              {
                mu <- as.matrix(beta.prior[['mean.mu']])
                
                # if one did not trim
                if ( (zero.intercept==TRUE) && ( length(mu)==(k+1) ) ) 
                  mu <- mu[2:(k+1)]
                
                if (length(mu)!=k) 
                  stop(paste('The dimension of mu-vector must be ',k,
                             ', \n which is the number of unknown regression coefficients.',sep=''))
                
              }
              else
              {
                warning("Mean vector is not supplied for beta prior. Zero vector is used instead.");
                mu <- as.matrix(rep(0.0,k))                
              }

              if(!is.null(beta.prior[['prec.Cinv']])&is.matrix(beta.prior[['prec.Cinv']])) # check if var inverse is supplied
              {
                
                # if one did not trim:
                if ((zero.intercept==TRUE) && (dim(beta.prior[['prec.Cinv']])[1]==(k+1)) && (dim(beta.prior[['prec.Cinv']])[2]==(k+1)))
                {
                  temp.C <- chol2inv(chol(beta.prior[['prec.Cinv']]))
                  temp.C <- temp.C[2:(k+1),2:(k+1)]
                  beta.prior[['prec.Cinv']] <- chol2inv(chol(temp.C))
                }                 
                
                if ((dim(beta.prior[['prec.Cinv']])[1]!=k) || (dim(beta.prior[['prec.Cinv']])[2]!=k))
                  stop(paste('The dimension of prec.Cinv matrix must be ',k,'x',k,'.',sep=''))
  
                var.inv <- beta.prior[['prec.Cinv']]
                if(!is.null(beta.prior[['cov.C']])) 
                  warning('beta.prior list contains both "prec.Cinv" and "cov.C" matrices. Only "prec.Cinv" will be used and "cov.C" will be ignored.');
              }
              else
              {
                if(!is.null(beta.prior[['cov.C']])&is.matrix(beta.prior[['cov.C']]))
                {
                  
                  # if one did not trim:                
                  if ((zero.intercept==TRUE) && (dim(beta.prior[['cov.C']])[1]==(k+1)) && (dim(beta.prior[['cov.C']])[2]==(k+1)))
                    beta.prior[['cov.C']] <- beta.prior[['cov.C']][2:(k+1),2:(k+1)]
                  
                  if ((dim(beta.prior[['cov.C']])[1]!=k) || (dim(beta.prior[['cov.C']])[2]!=k))
                    stop(paste('The dimension of cov.C matrix must be ',k,'x',k,'.',sep=''))
                  
                    var.inv <- chol2inv(chol(beta.prior[['cov.C']]))                  
                }                  
                else
                {
                  warning("Covariance matrix is not supplied for beta prior. Identity matrix is used instead.");
                  var.inv = diag(k);
                } 
              }
                
              switch(sqtype, 
                      inverse.gamma=# sigma squared distribution is Gamma(a,b). we will need two parameters a and b, starting value for sigma squared and preferably the inverse of xtx
                        {
                                                    
                          # set default values for inverse gamma prior parameters
                          
                          if (!is.null(sigmasq.prior[["inverse.gamma.a"]]))
                            inv.gamma.a <- as.numeric(sigmasq.prior[["inverse.gamma.a"]])
                          else inv.gamma.a <- 1.0
                          
                          if (!is.null(sigmasq.prior[["inverse.gamma.b"]])) 
                            inv.gamma.b <- as.numeric(sigmasq.prior[["inverse.gamma.b"]])
                          else inv.gamma.b <- 1.0
                          
                          if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                            sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                          else sigmasq.init <- 1.0
                          
                          
                          bayesregressB2S1(xtx,xty,yty,
                                           Tsamp.out = Tsamp.out,
                                           beta.prior.mean = mu,
                                           beta.prior.var=var,
                                           beta.prior.var.inv=var.inv,
                                           inv.gamma.a=inv.gamma.a, inv.gamma.b=inv.gamma.b,
                                           sigmasq.init=sigmasq.init,                                           
                                           numsamp.data=num.samp)
                        },
                      sigmasq.inverse=
                        {
                          
                          if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                            sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                          else sigmasq.init=1.0
                          
                          
                          bayesregressB2S2(xtx,xty,yty,numsamp.data=num.samp,
                                           sigmasq.init=sigmasq.init,
                                           beta.prior.mean = mu,
                                           beta.prior.var=var,
                                           beta.prior.var.inv=var.inv,                                           
                                           Tsamp.out=Tsamp.out)
                        },
                     stop("Prior distribution for sigma squared is not properly formatted")
              )
            },
            mvnorm.unknown={
              
              if (!is.null(beta.prior[['mu.hyper.mean.eta']]))
              {
                eta <- as.numeric(beta.prior[['mu.hyper.mean.eta']])
                
                # if one did not trim
                if ( (zero.intercept==TRUE) && ( length(eta)==(k+1) ) ) 
                  eta <- eta[2:(k+1)]
                
                if (length(eta)!=k) 
                  stop(paste('The dimension of mu.hyper.mean.eta-vector must be ',k,
                             ', \n which is the number of unknown regression coefficients.',sep=''))
                
              }                
              else
              {
                warning("Mean vector is not supplied for beta prior. Vector of 0's is used instead.");
                eta=as.numeric(rep(0.0,k));
              }
                
              
              if(!is.null(beta.prior[['mu.hyper.prec.Dinv']])&is.matrix(beta.prior[['mu.hyper.prec.Dinv']])) # check if var inverse is supplied
              {
                
                var.inv <- beta.prior[['mu.hyper.prec.Dinv']]
                
                # if one did not trim:
                if ((zero.intercept==TRUE) && (dim(var.inv)[1]==(k+1)) && (dim(var.inv)[2]==(k+1)))
                {
                  temp.D  <- chol2inv(chol(var.inv))
                  temp.D  <- temp.D[2:(k+1),2:(k+1)]
                  var.inv <- chol2inv(chol(temp.D))
                }                 
                
                if ((dim(var.inv)[1]!=k) || (dim(var.inv)[2]!=k))
                  stop(paste('The dimension of mu.hyper.prec.Dinv matrix must be ',k,'x',k,'.',sep=''))                
              }                
              else
              {
                warning("Variance matrix is not supplied for beta prior. Identity matrix is used instead.");
                var.inv=diag(k);
              }
               
              
              if (!is.null(beta.prior[['Cinv.hyper.df.lambda']]))
                wlambda<-as.numeric(beta.prior[['Cinv.hyper.df.lambda']])
              else
              {
                warning(paste("Wishart degrees freedom parameter is not supplied. Default value of lambda equal to the number of predictors (" ,k,") is used instead."));
                wlambda=k; 
              }

              
              if(!is.null(beta.prior[['Cinv.hyper.invscale.Vinv']])&is.matrix(beta.prior[['Cinv.hyper.invscale.Vinv']])) # check if var inverse is supplied
              {  
                wishart.V.inv <- beta.prior[['Cinv.hyper.invscale.Vinv']]
                
                # if one did not trim:
                if ((zero.intercept==TRUE) && (dim(wishart.V.inv)[1]==(k+1)) && (dim(wishart.V.inv)[2]==(k+1)))
                {
                  temp.V         <- chol2inv(chol(wishart.V.inv))
                  temp.V         <- temp.V[2:(k+1),2:(k+1)]
                  wishart.V.inv  <- chol2inv(chol(temp.V))
                }                 
                
                if ((dim(wishart.V.inv)[1]!=k) || (dim(wishart.V.inv)[2]!=k))
                  stop(paste('The dimension of Cinv.hyper.invscale.Vinv matrix must be ',k,'x',k,'.',sep=''))                
                
              }
              else 
              {
                warning("Wishart inverse variance matrix is not supplied for beta prior. Identity matrix is used instead.");
                wishart.V.inv=diag(k);
              }
              
              if (!is.null(beta.prior[['mu.init']]))
              {
                mu.init <- as.numeric(beta.prior[['mu.init']])
                
                # if one did not trim:
                if ( (zero.intercept==TRUE) && ( length(mu.init)==(k+1) ) ) 
                  mu.init <- mu.init[2:(k+1)]
                
                if (length(mu.init)!=k) 
                  stop(paste('The dimension of mu.init-vector must be ',k,
                             ', \n which is the number of unknown regression coefficients.',sep=''))
                                
              }                
              else
              {
                warning("Initial value for mu vector is not supplied for beta prior. Vector of 1's is used instead.");
                mu.init=as.numeric(rep(1.0,k));
              }
              
              if(!is.null(beta.prior[['Cinv.init']] )& is.matrix(beta.prior[['Cinv.init']])) # check if var inverse is supplied
              {
                Cinv.init<-beta.prior[['Cinv.init']] 
                
                # if one did not trim:                
                if ((zero.intercept==TRUE) && (dim(Cinv.init)[1]==(k+1)) && (dim(Cinv.init)[2]==(k+1)))
                {
                  temp.C     <- chol2inv(chol(Cinv.init))
                  temp.C     <- temp.C[2:(k+1),2:(k+1)]
                  Cinv.init  <- chol2inv(chol(temp.C))
                }
                                  
                if ((dim(Cinv.init)[1]!=k) || (dim(Cinv.init)[2]!=k))
                  stop(paste('The dimension of Cinv.init matrix must be ',k,'x',k,'.',sep=''))
                
              }                
              else{
                warning("Initial value for variance matrix inverse is not supplied for beta prior. Identity matrix is used instead.")
                Cinv.init=diag(k);
              }
              
              switch(sqtype, 
                     inverse.gamma=# sigma squared distribution is Gamma(a,b). we will need two parameters a and b, starting value for sigma squared and preferably the inverse of xtx
                          {                            
                                                                                    
                            # set default values for inverse gamma prior parameters
                            
                            if (!is.null(sigmasq.prior[["inverse.gamma.a"]]))
                              inv.gamma.a <- as.numeric(sigmasq.prior[["inverse.gamma.a"]])
                            else inv.gamma.a <- 1.0
                            
                            if (!is.null(sigmasq.prior[["inverse.gamma.b"]])) 
                              inv.gamma.b <- as.numeric(sigmasq.prior[["inverse.gamma.b"]])
                            else inv.gamma.b <- 1.0
                            
                            if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                              sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                            else sigmasq.init <- 1.0
                            
                            
                            bayesregressB3S1(xtx,xty,yty,numsamp.data=num.samp,                                             
                                             eta       = eta,
                                             Rinv      = var.inv, # Dinv
                                             mu.init   = mu.init,
                                             lambda    = wlambda,
                                             Vinv      = wishart.V.inv, # Cinv                                            
                                             Cinv.init = Cinv.init,
                                             inv.gamma.a = inv.gamma.a, 
                                             inv.gamma.b = inv.gamma.b,
                                             sigmasq.init=sigmasq.init,                                             
                                             Tsamp.out = Tsamp.out)
                          },
                          sigmasq.inverse=
                          {
                            
                            if(!is.null(sigmasq.prior[["sigmasq.init"]])) 
                              sigmasq.init <- as.numeric(sigmasq.prior[["sigmasq.init"]])
                            else sigmasq.init = 1.0                            
                            
                            bayesregressB3S2(xtx, xty, yty, numsamp.data=num.samp,                                             
                                             eta=eta,
                                             Rinv=var.inv,
                                             mu.init=mu.init,
                                             lambda=wlambda,
                                             Vinv=wishart.V.inv,
                                             Cinv.init=Cinv.init,
                                             sigmasq.init=sigmasq.init,
                                             Tsamp.out=Tsamp.out)
                          },
                          stop("Prior distribution for sigma squared is not properly formatted")
              )
            }

  )
  return(result)
}
