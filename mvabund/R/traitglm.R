traitglm = function( L, R, Q=NULL, family="negative.binomial", formula = NULL, method="manyglm", composition=FALSE, col.intercepts = TRUE, ...  )
{

#subfunctions get.design and get.polys defined below.
  
  # extract any arguments that work with cv.glm1path and save separately so they are not passed to glm1path.
  L = as.data.frame(L)
  allargs <- match.call(expand.dots = FALSE)
  dots <- allargs$...
  if( "best" %in% names(dots) )
    best <- dots$best    
  else
    best = "1se"
  if( "plot" %in% names(dots) )
    plot <- dots$plot    
  else
    plot=TRUE
  if( "prop.test" %in% names(dots) )
    prop.test <- dots$prop.test
  else
    prop.test = 0.2
  if( "n.split" %in% names(dots) )
    n.split <- dots$n.split
  else
    n.split=10
  if( "seed" %in% names(dots) )
    seed <- dots$seed
  else
    seed=NULL
  if( "show.progress" %in% names(dots) )
    show.progress <- dots$show.progress
  else
    show.progress = FALSE
  if( "get.fourth" %in% names(dots) )
    get.fourth <- dots$get.fourth
  else
    get.fourth = TRUE
  
  deactive <- c("best", "plot", "prop.test", "n.split", "seed", "show.progress", "get.fourth") 
  deactivate <- (1:length(dots))[names(dots) %in% deactive ]  
  for (i in length(deactivate):1) 
    dots[ deactivate[i] ]<-NULL
  
  dots <- lapply( dots, eval, parent.frame() )
  
  
  n.sites = dim(L)[1] #use R.des gives number of sites for prediction / model fitting
  n.spp   = dim(L)[2]

  # get standardised R, Q, orthogonal polys, and poly coeffs
  if(is.null(formula))
    R.des = get.polys(R)
  else
    R.des = list(X=R)
    
  if(is.null(Q))
  {
      cat(paste("No traits matrix entered, so will fit SDMs with different env response for each spp","\n"))
      Q.des = list(X=data.frame(names(L)),X.squ=NULL,var.type="factor")
  }
  else
  {
    if(is.null(formula))
      Q.des = get.polys(Q)
    else
      Q.des = list(X=Q)
  }  

  any.penalty = method=="cv.glm1path" || method=="glm1path"
  # get mega-matrix of design for regression against vectorised l
  marg.penalty = TRUE
  X.des = get.design( R.des, Q.des, names(L), formula=formula, marg.penalty=marg.penalty, composition = composition, col.intercepts = col.intercepts, any.penalty=any.penalty, get.fourth=get.fourth )
  # setting marg.penalty=TRUE means that spp are always in the model, and penalised if any.penalty=TRUE 
  X = X.des$X
  l <- as.vector(as.matrix(L))
  
  if( method=="cv.glm1path" || method=="glm1path" )
  {
    if( "block" %in% names(dots) )
      blockID <- dots$block    
    else
      blockID = 1:n.sites
    block = factor(rep(blockID,n.spp))
    ft = do.call(glm1path,c(list(y=l, X=X, family=family, penalty = X.des$penalty), k=log(n.sites), dots))
    if( method=="cv.glm1path" )
      ft = do.call(cv.glm1path,c(list(object=ft, block=block, best=best, plot=plot, prop.test=prop.test, n.split=n.split, 
                                    seed=seed, show.progress=show.progress), dots))
    id.use = which(ft$lambdas==ft$lambda)
    ft$deviance = -2*ft$logL[id.use]
    ft$phi = ft$glm1$phi
    if(ft$df[1]==1)
      null.deviance = ft$logL[1]
#    ft$family=ft$glm1$family
  }
  else
    ft = do.call( method, c(list(formula=l~., family=family, data=data.frame(X)), dots) )
  # note user-entered formula was not used in fitting, because X is the model matrix generated from this formula.

  # report fourth corner matrix for "best" model
  ft$fourth.corner = matrix( coef(ft)[X.des$is.4th.corner], length(X.des$names.Q), length(X.des$names.R) )
  ft$fourth.corner = provideDimnames(ft$fourth,base=list(X.des$names.Q,X.des$names.R))

#    ft$fourth.corner = Matrix(fourth.corner, sparse=T)
  # cat("\n")
  # cat("Fourth corner matrix:    ")
  # cat("\n")
  # print(round(ft$fourth.corner,dec.pl))


  ### Plot results

  # do a lattice "levelplot" of the coefficients in the matrix of fourth corner interactions
#  if( fourthPlot == TRUE & get.fourth == TRUE )
#  {
#      library(lattice)
#      a        = max( abs(ft$fourth.corner) )
#      if(a < 1.e-8)
#          warning("No fourth corner interactions were included in final model")
#      else
#      {
#          colort   = colorRampPalette(c("blue","white","red")) 
#          plot.4th = levelplot(t(as.matrix(ft$fourth.corner)), xlab="Environmental Variables", ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100), scales = list( x= list(rot = 45)))
#          print(plot.4th)
#      }
#  }
    
  ft$R.des = R.des
  ft$Q.des = Q.des
  ft$any.penalty = any.penalty
  ft$L = L
  ft$scaling = X.des$scaling
  ft$call=match.call()

  if(is.null(formula)==FALSE & get.fourth==FALSE) # if formula was user-entered but not we exclude 4th corner, redefine it
    ft$formula = X.des$formula
  else
    ft$formula = formula

  class(ft)=c("traitglm",class(ft))
  return( ft )
    
}



################ get.polys for getting orthogonal polynomials ###################
     
get.polys = function( X, X.des.train=NULL)
{
# get.polys will take a matrix of env vars (or trait vars), and standardise the quant ones
# as well as return orthogonal poly's. Importantly, if training matrices are given as input,
# these will be used in matrix construction.

    
    if(is.null(dim(X)))
      X = data.frame(X)
    n.sites = dim(X)[1]
    n.params  = dim(X)[2]
    if(is.null(X.des.train))
        n.train.sites = n.sites
    else
        n.train.sites <- dim(X.des.train$X)[1]
    if(is.null(X.des.train$var.type))
        var.type = rep("quantitative",n.params)
    else
        var.type = X.des.train$var.type
    for (i in 1:n.params)
    {

        # test if binary quantitative, if so, change to a factor to avoid poly error.  But only if training data
        if(is.null(X.des.train$var.type) & is.factor(X[,i])==FALSE)
        {
          testBinary = try(duplicated(X[,i],nmax=2), silent=TRUE)
          if(class(testBinary)=="logical")
          {
            X[,i] = factor(X[,i])
            warning(paste0("Binary variable '", names(X)[i], "' found and changed to a factor"))
          }
        }
        
        if(is.factor(X[,i]))
        {
            n.levels    = length(levels(X[,i]))
            if(n.levels==2)
            {
                var.type[i]="binary" #treat as quantitative but don't find its square
                #change variable name to indicate what it is doing
                dimnames(X)[[2]][i]=paste(names(X)[i],levels(X[,i])[2],sep="")
                #change entry in X to numerical
                X[,i]  = as.numeric(as.factor(X[,i]))*2-3
            }
            else
            {
                var.type[i]="factor"
                contrasts(X[,i])[contrasts(X[,i])==0] = -1
            }
        }
    }

    # to return standardised values of quant vars, with coeff, where needed:
    is.quant = which(var.type=="quantitative")
    n.quant = length(is.quant)
    if( n.quant>0 )
    {
        X.squ = X[,0]
        degs = c()
        names.X.squ = c()
        poly.coefs = as.list( rep( NA, n.quant ) )
        for(i.quant in is.quant)
        {
            poly.i = poly( X[,i.quant], degree=2, coefs=X.des.train$coefs[[i.quant]] )
            X.squ = cbind( X.squ, poly.i )
            degs  = c( degs, attr(poly.i, "degree") )
            poly.coefs[[ i.quant ]] = attr(poly.i, "coefs")   
            names.X.squ = c( names.X.squ, dimnames(X)[[2]][i.quant], paste( dimnames(X)[[2]][i.quant], ".squ", sep="") )
        }
        X.squ = X.squ * sqrt(n.train.sites)
        dimnames(X.squ)[[2]] = names.X.squ
        X[,var.type=="quantitative"] = X.squ[,degs==1]
        #get rid of the linear terms:
        X.squ = data.frame(X.squ[,degs==2])
    }
    else
    {
        X.squ=NULL
        poly.coefs=NULL
    }
    # to return orthogonal poly values of quant vars (with coeff):
    return( list( X=X, X.squ=X.squ, var.type=var.type, coefs=poly.coefs ) )
}


################ get.design for getting the design matrix ###################
get.design = function( R.des, Q.des, L.names, formula = formula, marg.penalty=TRUE, composition = FALSE, col.intercepts = TRUE, any.penalty=TRUE, scaling=NULL, get.fourth=TRUE )
{

# get.design will take matrices of linear env and trait terms, and orthogonal quadratic terms,
# and return a mega-matrix that can be regressed against vectorised abundance.
# also returns logical for fourth corner terms, and their row and column names

# unless of course a formula is specified, in which case it will use that to build model matrix.
  
    #How many site are there? This will be used later
  n.sites = dim(R.des$X)[1] #use R.des gives number of sites for prediction / model fitting
  n.spp   = length(L.names)

  is.scaling.given = is.null(scaling)==F #remember if scaling was specified - if not, build it up
  #get spp indicator
  if(col.intercepts==TRUE)
  {
      spp     = rep(L.names,each=n.sites)
      spp     = as.factor(spp)
      mod     = as.formula("~spp-1")
      X.spp   = model.matrix(mod)
      if(marg.penalty==FALSE || any.penalty==FALSE) #if species not penalised need to remove species 1 indicator
        X.spp   = X.spp[,-1]
      X.spp[X.spp==0] = -1
      if(is.scaling.given==F) #rescale X.spp so ALL variables have variance 1.
      {
        scaling = list()
        X.spp = scale(X.spp)
        scaling$spp$center = attr(X.spp,"scaled:center")
        scaling$spp$scale  = attr(X.spp,"scaled:scale")      
      }
      if(is.scaling.given)
        X.spp = scale(X.spp,center=scaling$spp$center, scale=scaling$spp$scale)
      X.spp   = cbind(1,X.spp) #add intercept since my code doesn't do that automatically
  }
  else
  {
    X.spp = as.matrix( rep(1,n.sites*n.spp) )
      if(is.scaling.given==F) #if required, initiate scaling object (to be built later) 
        scaling = list() 
  }
  if(composition==TRUE)
  {
    site    = rep(dimnames(R.des$X)[[1]],n.spp)
    site    = as.factor(site)
    mod     = as.formula("~site-1")
    X.site  = model.matrix(mod)
    if(marg.penalty==FALSE || any.penalty==FALSE) #if site not penalised need to remove site 1 indicator
      X.site = X.site[,-1]
    X.site[X.site==0] = -1
    if(is.scaling.given==F) #rescale X.site so ALL variables have variance 1.
    {
      X.site  = scale(X.site)
      scaling$site$center = attr(X.site,"scaled:center")
      scaling$site$scale  = attr(X.site,"scaled:scale")      
    }
    if(is.scaling.given)
      X.site = scale(X.site,center=scaling$site$center, scale=scaling$site$scale)
  }
  else #otherwise make it empty to avoid errors later
    X.site = X.spp[,0]

  if(is.null(formula)) #if no formula provided, build up a design matrix by default with quad terms etc
  {
    
    # R terms
    X.R     = X.spp[,0]
    if( is.null(R.des$X.squ) )
    {
        R.small = R.des$X
        var.type= R.des$var.type
        is.lin.small = rep( TRUE,NCOL(R.des$X) )
    }
    else
    {
        R.small = cbind( R.des$X, R.des$X.squ )
        var.type= c( R.des$var.type, rep("quantitative",dim(R.des$X.squ)[2]) )
        is.lin.small = c( rep( TRUE,NCOL(R.des$X) ), rep( FALSE, NCOL(R.des$X.squ) ) )
    }    
    names.R = c()
    is.lin.R= c()
    for( iR in 1:NCOL(R.small) )
    {
        R.i   = rep( R.small[,iR], times=n.spp )
        mod     = as.formula("~0+R.i")
        mm      = model.matrix(mod)
        names.i = dimnames(R.small)[[2]][iR]
        if(var.type[iR]=="factor")
        {
            if(any.penalty==TRUE)
              names.i = paste(dimnames(R.small)[[2]][iR], levels(R.small[,iR]),sep="")
            else #if unpenalised, remove the first level
            {
              mm = mm[,-1]
              names.i = paste(dimnames(R.small)[[2]][iR], dimnames(contrasts(R.small[,iR]))[[2]],sep="")
            }
        }
        names.R = c( names.R, names.i )
        is.lin.R= c( is.lin.R, rep( is.lin.small[iR] , dim(mm)[2] ) )
        X.R     = cbind( X.R, mm )
    }
    dimnames(X.R)[[2]]=names.R
    if(is.scaling.given==F) #rescale X.R so ALL variables have variance 1.
    {
      X.R = scale(X.R)
      scaling$R$center = attr(X.R,"scaled:center")
      scaling$R$scale  = attr(X.R,"scaled:scale")      
    }
    if(is.scaling.given)
      X.R = scale(X.R,center=scaling$R$center, scale=scaling$R$scale)
    
    # Q terms
    X.Q     = X.spp[,0]
    if( is.null(Q.des$X.squ) )
    {
        Q.small = Q.des$X
        var.type= Q.des$var.type
        is.lin.small = rep( TRUE,NCOL(Q.des$X) )
    }
    else
    {
        Q.small = cbind( Q.des$X, Q.des$X.squ )
        var.type= c( Q.des$var.type, rep("quantitative",dim(Q.des$X.squ)[2]) )
        is.lin.small = c( rep( TRUE,NCOL(Q.des$X) ), rep( FALSE, NCOL(Q.des$X.squ) ) )
    }    
    names.Q = c()
    is.lin.Q= c()
    for( iQ in 1:NCOL(Q.small) )
    {
        Q.i  = rep( Q.small[,iQ], each=n.sites )
        mod     = as.formula("~0+Q.i")
        mm      = model.matrix(mod)
        names.i = dimnames(Q.small)[[2]][iQ]
        if(var.type[iQ]=="factor")
        {
          if(any.penalty==TRUE)
            names.i = paste(dimnames(Q.small)[[2]][iQ], levels(Q.small[,iQ]),sep="")
          else #if unpenalised, remove the first level
          {
            mm = mm[,-1]
            names.i = paste(dimnames(Q.small)[[2]][iQ], dimnames(contrasts(Q.small[,iQ]))[[2]],sep="")
          }
        }
        names.Q = c( names.Q, names.i )
        is.lin.Q= c( is.lin.Q, rep( is.lin.small[iQ] , dim(mm)[2] ) )
        X.Q     = cbind( X.Q, mm )
    }
    dimnames(X.Q)[[2]]=names.Q
    if(is.scaling.given==F)
    {
      X.Q = scale(X.Q)
      scaling$Q$center = attr(X.Q,"scaled:center")
      scaling$Q$scale  = attr(X.Q,"scaled:scale")      
    }
    if(is.scaling.given)
      X.Q = scale(X.Q, center=scaling$Q$center, scale=scaling$Q$scale)
    
    if(get.fourth==TRUE)
    {
      # R*Q interaction
      X.RQ    = X.spp[,0]
      n.R = sum(is.lin.R)
      n.Q = sum(is.lin.Q)
      ref.R = rep(1:n.R,each=n.Q)
      ref.Q = rep(1:n.Q,n.R)
      X.RQ = as.matrix( X.R[,ref.R] * X.Q[,ref.Q] )
      dimnames(X.RQ)[[2]] = paste(dimnames(X.R)[[2]][ref.R], dimnames(X.Q)[[2]][ref.Q], sep=":")
      if(is.scaling.given==F)
      {
        X.RQ = scale(X.RQ)
        scaling$RQ$center = attr(X.RQ,"scaled:center")
        scaling$RQ$scale  = attr(X.RQ,"scaled:scale")      
      }
      if(is.scaling.given)
        X.RQ = scale(X.RQ,center=scaling$RQ$center, scale=scaling$RQ$scale)
    }
    else
      X.RQ=X.R[,0]
    
    if(any.penalty)
    {
      if(composition==TRUE)
        X             = cbind(X.spp,X.site,X.R,X.Q,X.RQ)
      else
        X             = cbind(X.spp,X.R,X.Q,X.RQ)
    }
    else
    {  
      if(composition==TRUE)
        X             = cbind(X.spp,X.site,X.RQ) #no trait nor env main effects if there is no penalty on trait params, covered by spp and site.
      else
        X             = cbind(X.spp,X.R,X.RQ) #no trait main effects if there is no penalty on trait params, covered by spp.
    }
    
    n.X           = dim(X)[2]
    if(get.fourth==TRUE)
      is.4th.corner = c( rep(F,n.X-n.R*n.Q), rep(T,n.R*n.Q) )
    else
      is.4th.corner = rep(F,n.X)
    names.R = dimnames(X.R)[[2]][is.lin.R]
    names.Q = dimnames(X.Q)[[2]][is.lin.Q]
  }
  else # if formula has been provided, use it to build up design matrix
  {
    # construct design matrix X (minus any spp and site terms)
    rowReps = rep(1:n.sites, times=n.spp)
    R.i     = data.frame(R.des$X[rowReps,])
    names(R.i) = names(R.des$X)
    colReps = rep(1:n.spp, each=n.sites)
    Q.i     = data.frame(Q.des$X[colReps,])
    names(Q.i) = names(Q.des$X)
    X       = model.matrix(formula,cbind(R.i,Q.i))
    
    # now extract fourth corner as a matrix with just one column
    tt       = terms(formula)
    facts    = attr(tt,"factors")
    which.env   = charmatch(names(R.des$X),dimnames(facts)[[1]])
    which.env   = which.env[is.na(which.env)==FALSE]
    which.trait = charmatch(names(Q.des$X),dimnames(facts)[[1]])
    which.trait = which.trait[is.na(which.trait)==FALSE]
    if(length(which.env)==1) #get a vector saying which terms involve env variables
      is.env = facts[which.env,]>0
    else
    {
      if(length(which.env)==0)
         is.env = logical(0)
      else
         is.env = apply(facts[which.env,],2,sum)>0
    }
    if(length(which.trait)==1) #get a vector saying which terms involve traits
      is.trait = facts[which.trait,]>0
    else
    {
      if(length(which.trait)==0)
        is.trait = logical(0)
      else
        is.trait = apply(facts[which.trait,],2,sum)>0
    }
    is.4th.term = is.env & is.trait # 4th corner terms (involve both env and traits)
    if(get.fourth==TRUE)
      is.4th.corner = is.4th.term[attr(X,"assign")] #logical for whether columns of X are 4th corner, missing a FALSE for intercept though
    else
      is.4th.corner = rep( FALSE, length(attr(X,"assign")) ) #set them all to FALSE if no 4th corner terms to be included in model

    if(attr(tt,"intercept")==1) # get rid of intercept in X so everything matches up
    {
      names.Q = dimnames(X)[[2]][c(FALSE,is.4th.corner)]
      X = as.matrix(X[,-1]) # as.matrix include in case it is reduced to a vector
    }
    else
      names.Q = dimnames(X)[[2]][is.4th.corner]
    names.R = "coef"
    
    # if no 4th corner terms in model, kick them out of X, make is.4th.corner empty, reconstruct formula to ditch 4th terms
    if(get.fourth==FALSE)
    {
      X = as.matrix(X[,-which(is.4th.term==TRUE)])
      is.4th.corner = c()
      if(sum(is.4th.term==FALSE)>0)
      {
        formula = paste0(dimnames(facts)[[2]][is.4th.term==FALSE],collapse="+")
        formula = as.formula(paste0("~",formula))
      }
      else
        formula=as.formula("~1")
    }
    # add on spp and site terms, and add to is.4th.corner too:
    X = cbind(X.spp,X.site,X)
    is.4th.corner = c(rep(FALSE, dim(X.spp)[2]+dim(X.site)[2]), is.4th.corner)
  }
  if(any.penalty) #set the value of penalty, if using glm1path
  {
    penalty = c( 0, rep(1,dim(X)[2]-1) )
    if (col.intercepts==TRUE & marg.penalty==FALSE) #change if spp intercept terms are to be left unpenalised:
      penalty  = c( rep( 0,dim(X.spp)[2]+dim(X.site)[2] ), rep( 1, dim(X)[2]-dim(X.spp)[2]-dim(X.site)[2] ) )
  }
  else #otherwise penalty can be empty
    penalty = NULL

  return(list(X=X, is.4th.corner=is.4th.corner, names.R=names.R, names.Q=names.Q, penalty=penalty, any.penalty=any.penalty, scaling=scaling, formula=formula) )  
}
