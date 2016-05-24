lcmixed <- function ( formula = .~. , continuous, discrete, ppdim,
                     diagonal = TRUE, pred.ordinal=FALSE, printlik=FALSE )
{
#  require ("mvtnorm")
#  require(flexmix)
  retval <- new ("FLXMC", weighted = TRUE ,
                 formula = formula , dist = "mixnormmulti",
                 name = "latent class for normal/multinomial mixed data")
  retval@defineComponent <- expression ({
    logLik <- function (x, y) {
#      print("LogLik")
      if (continuous!=0){
        out <- dmvnorm (as.matrix(y)[,1:continuous,drop=FALSE], mean = center ,
               sigma = cov , log = TRUE )
      }
      else
        out <- 0
      if (discrete>0){
        for (k in 1:discrete)
          out <- out+log(pp[[k]][y[,continuous+k]])
      }
      if (printlik){
        cat("LogLikelihood= ",sum(out), "\n")
        cat("pp= ",pp[[k]],"\n")
      }
#      print("end Loglik")
      out
    }
    predict <- function (x) {
#      print("predict")
      if (!is.null(center))
        out <- matrix (center , nrow = nrow (x),
                     ncol = length ( center ), byrow = TRUE )
      else
        out <- matrix(0,ncol=2,nrow=nrow(x))
      if (discrete>0){
        for (k in 1:discrete){
          if (pred.ordinal)
            out[,continuous+k] <- sum((1:ppdim[k])*pp[[k]])
          else  
            out[,continuous+k] <- which.max(pp[[k]])
        }
      }
#      print("end predict")
      out
    }
    if (continuous==0){
      center <- NULL
      cov <- NULL
    }
    new ("FLXcomponent",
      parameters = list ( center = center , cov = cov , pp=pp),
      df = df , logLik = logLik , predict = predict )
  })
  retval@fit <- function (x, y, w) {
#    print("fit") 
    n <- nrow(x)
    sw <- sum(w)
#    print(sum(w[y[,3]==1]))
#    print(sw)
    if (continuous!=0){
      para <- cov.wt(as.matrix(y)[,1:continuous,drop=FALSE], wt = w)[c("center", "cov")]
    }
    else{
      para <- list(center=NULL,cov=NULL)
    }
    para$pp <- list()
    if (discrete>0){
      for (k in 1:discrete){
        para$pp[[k]] <- numeric(0)
#        print(table(as.matrix(y)[,continuous+k]))
        for (l in 1:ppdim[k])
          para$pp[[k]][l] <- sum(w*(as.matrix(y)[,continuous+k]==l))/sw
      }
    }
    df <- (3 * continuous + continuous^2)/2
    if (discrete>0)
      df <- df + sum(ppdim-1)
    if (continuous>0){
      if ( diagonal ) {
        if (ncol(para$cov)>1)
          para$cov <- diag ( diag ( para$cov ))
        df <- 2 * continuous + sum(ppdim-1)
      }
    }
#    print("end fit")
    with (para , eval ( retval@defineComponent ))
  }
  retval
}


# Recodes discrete/categorical variables to standard numbering
discrete.recode <- function(x,xvarsorted=TRUE,continuous=0,discrete){
  xold <- x
  x <- data.matrix(x)
#  str(x)
  if (!xvarsorted){
    x <- x[,c(continuous,discrete)]
    xold <- xold[,c(continuous,discrete)]
    continuous <- length(continuous)
    discrete <- length(discrete)
  }
  x.recode <- matrix(0,ncol=ncol(x),nrow=nrow(x))
  if (continuous>0)
    x.recode[,1:continuous] <- x[,1:continuous]
#  str(x.recode)
  discretelevels <- list()
  ppdim <- numeric(0)
  for (i in 1:discrete){
      discretelevels[[i]] <- levels(as.factor(xold[,continuous+i]))
      ppdim[i] <- length(discretelevels[[i]])
#      print(i)
#      print(discretelevels[[i]])
      for (j in 1:ppdim[i])
        x.recode[as.factor(xold[,continuous+i])==discretelevels[[i]][j],
                 continuous+i] <- j
  }
  out <- list(data=x.recode,ppdim=ppdim,discretelevels=discretelevels,
              continuous=continuous,discrete=discrete)
  out
}
  
# xvarsorted: TRUE if continuous variables come first and then the discrete
# ones. In this case continuous and discrete give numbers of these variables.
# Otherwise they are vectors indicating sets of variable numbers.

flexmixedruns <- function(x,diagonal=TRUE,xvarsorted=TRUE,
                          continuous,discrete,ppdim=NULL,initial.cluster=NULL,
                          simruns=20,n.cluster=1:20,verbose=TRUE,recode=TRUE,
                          allout=TRUE,control=list(minprior=0.001),silent=TRUE){
  x <- as.matrix(x)
#  require(flexmix)
  if (recode){
    drx <- discrete.recode(x,xvarsorted=xvarsorted,
                           continuous=continuous,discrete=discrete)
    continuous <- drx$continuous
    discrete <- drx$discrete
    x <- drx$data
    ppdim <- drx$ppdim
    discretelevels <- drx$discretelevels
#  print(str(x))
#  print(ppdim)      
  }
  else{
    if (!xvarsorted){
      x <- x[,c(continuous,discrete)]
      continuous <- length(continuous)
      discrete <- length(discrete)
    }
    discretelevels <- list()
    if (!identical(discrete,0)){
      for (i in 1:discrete)
        discretelevels[[i]] <- 1:ppdim[i]
    }
  }
  minBIC <- Inf
  flexout <- list()
  bicvals <- rep(NA,max(n.cluster))
  optimalk <- 0
  errcount <- rep(NA,max(n.cluster))
  for (k in n.cluster){
    errcount[k] <- 0
    maxlik <- -Inf
    flexout[[k]] <- NA
    for (i in 1:simruns){
      fmo <- try(flexmix(x~1,k=k,
               model=lcmixed(continuous=continuous,discrete=discrete,
                 ppdim=ppdim,diagonal=diagonal),cluster=initial.cluster,
               control=control),silent=silent)
      if (class(fmo)=="try-error"){
        if (verbose) cat("k= ",k," flexmix error in run ",i," \n")
        errcount[k] <- errcount[k]+1
      }
      else{
        if(fmo@logLik>maxlik){
          maxlik <- fmo@logLik
          flexout[[k]] <- fmo
          if (verbose) cat("k= ",k," new best fit found in run ",i,"\n")
        }
        else
          if (verbose) cat("Nonoptimal or repeated fit found in run ",i,"\n")
      }
      if (k==1) break
    }
    if (!identical(flexout[[k]],NA)){
      bicvals[k] <- BICk <- BIC(flexout[[k]])
      if (verbose) cat("k= ",k," BIC= ",BICk,"\n")
      if (BICk<minBIC){
        optimalk <- k
        minBIC <- BICk
      }
    }
  }
  optsummary <- summary(flexout[[optimalk]])
  if(!allout)
    flexout <- flexout[[optimalk]]
  out <- list(optsummary=optsummary,optimalk=optimalk,
              errcount=errcount,flexout=flexout,bicvals=bicvals,
              ppdim=ppdim,discretelevels=discretelevels)
  out      
}

# recodes categorical variables, positions given in categorical,
# as binary variables (one for every category).
cat2bin <- function(x,categorical=NULL){
  x <- as.matrix(x)
  variableinfo <- list()
  dimnam <- c()
  px <- ncol(x)
  vn <- 1
  x.con <- matrix(0,nrow=nrow(x),ncol=2)
  for (i in 1:px){
    if (i %in% categorical){      
      variableinfo[[i]] <- list()
      variableinfo[[i]]$type <- "categorical recoded"
      variableinfo[[i]]$levels <- levels(as.factor(x[,i]))
      variableinfo[[i]]$ncat <- length(variableinfo[[i]]$levels)
      variableinfo[[i]]$varnum <- vn:(vn+variableinfo[[i]]$ncat-1)
      for (j in 1:(variableinfo[[i]]$ncat)){
        dimnam[vn+j-1] <- paste("var",i,"cat",j,sep="")
        if (vn==1 & j==1)
          x.con <- as.integer(x[,i]==variableinfo[[i]]$levels[j])
        else
          x.con <- cbind(x.con,as.integer(x[,i]==variableinfo[[i]]$levels[j]))
      }
      vn <- vn+(variableinfo[[i]]$ncat)
    }
    else{
      dimnam[vn] <- paste("var",i,sep="")
      if (vn==1)
        x.con <- x[,i]
      else
        x.con <- cbind(x.con,x[,i])
      variableinfo[[i]] <- list()
      variableinfo[[i]]$type <- "not recoded"
      variableinfo[[i]]$varnum <- vn
      vn <- vn+1
    }
  }
#  print(x.con)
#  print(dimnam)
#  print(variableinfo)
  dimnames(x.con) <- list(NULL,dimnam)
  out <- list(data=x.con,variableinfo=variableinfo)
}

# Gives out data matrix with variables jitterv jittered
jittervar <- function(x,jitterv=NULL,factor=1){
  x.out <- x
  for (i in jitterv)
    x.out[,i] <- jitter(x[,i],factor=factor)
  x.out
}

# Computes the factor by which dummy or rank variables for type="categorical"
# or type="ordinal" variables should be multiplied when standardising
# them for comparison with continuous variables with variance 1.
# cat: number of categories.
# n: number of points. catsizes: number of observations per category.
# One of n and catsizes must be supplied.
# The normfactor 2 is E(X_1-X_2)^2 for X N(0,1).
# The qfactor is for relative weighting of average squared dissimilarities
# in cat/ordinal variables vs. continuous ones to account for discretisation.
# In the computation of categorical values there is a factor 1/sqrt(2),
# which accounts for the fact that dummy variables code categorial
# differences twice, which becomes sqrt(2) in Euklidean distance calculation.
# The idea is that the effective squared distance between the two levels of a
# binary variable, if catsizes=(n/2,n/2), is E(X_1-X_2)^2 above. Binaries are
# categorical and ordinal (though with categoricals "Euklidisation" of
# two dummies needs to be taken into account) and this is generalised in
# different ways for ordinals and categoricals.

# If distances shouldn't be determined data-driven, it makes sense not to
# specify catsizes and take n/cat by default.
distancefactor <- function(cat,n=NULL,
                           catsizes=NULL,type="categorical",
                           normfactor=2,
                           qfactor=ifelse(type=="categorical",1/2,
                             1/(1+1/(cat-1)))){
  if (is.null(catsizes))
    catsizes <- rep(round(n/cat),cat)
  n <- sum(catsizes)
  overall <- (n+1)*n/2
  if (type=="categorical"){
    wit <- 0
    for (i in 1:cat)
      wit <- wit+(catsizes[i]+1)*catsizes[i]/2
    bet <- overall-wit
    out <- sqrt(qfactor*normfactor*overall/(2*bet))
#    print(qfactor)
#    print(overall/bet)
  }
  if (type=="ordinal"){
    bd <- 0
    for (i in 1:(cat-1))
      for(j in (i+1):cat)
        bd <- bd+(j-i)^2*catsizes[i]*catsizes[j]
    out <- sqrt(qfactor*normfactor*overall/bd)
  }
  out
}


# clustering is assumed to run from 1 to max(clustering)
# vardata: data with variables to be summarised
# contdata: data on which projections are based
# tablevar: numbers of variables for which tables are used (categorical)
# catvar: variables that should be categorised for variablewise tables output
# quantvar: variables for which mean, sd and some quantiles are given 
# projmethod: projection method for discrproj
# rangefactor: how far outside the cluster range should the projection
# range go to both sides ("1" means that the cluster is in the middle 1/3.)
cluster.varstats <- function(clustering,vardata,contdata=vardata,
                             clusterwise=TRUE,
                            tablevar=NULL,catvar=NULL,
                             quantvar=NULL, catvarcats=10,
                             proportions=FALSE,
                            projmethod="none",minsize=ncol(contdata)+2,
                             ask=TRUE,
                            rangefactor=1){
  p <- ncol(vardata)
  vnames <- names(vardata)
  if (!(length(vnames)==p)){
    vnames <- c()
    for (i in 1:p)
      vnames[i] <- paste("Variable",i)
  }
  n <- length(clustering)
  if (clusterwise){
    if (ask)
      par(ask=TRUE)
    for (i in 1:max(clustering)){
      clusteri <- (clustering==i)
      cat("\nCluster ",i," ",sum(clusteri)," out of ",n," points.\n\n")
      for (j in 1:p){
        cat("Cluster ",i," ",vnames[j],"\n")
        if (j %in% tablevar){
          tt <- table(clusteri,vardata[,j],
                dnn=c(paste("In cluster ",i),vnames[j]))
          if (proportions)
            print(prop.table(tt,1))
          else
            print(tt)
          cat("\n")
        } # if tablevar
        else{
          par(mfrow=c(2,1))
          h1 <- hist(vardata[,j],main=paste("All obs.", vnames[j]),
                     xlab=vnames[j])
          h2 <- hist(vardata[clusteri,j],main=paste("Cluster ",i,vnames[j]),
                     breaks=h1$breaks,
                     xlab=vnames[j])
        } # else (!tablevar)
        if (j %in% quantvar){
          cat("  Mean=", mean(vardata[clusteri,j])," all obs.=",
              mean(vardata[,j]),"\n")
          cat("  Standard deviation=", sd(vardata[clusteri,j])," all obs.=",
              sd(vardata[,j]),"\n")
          print(quantile(vardata[clusteri,j], probs = seq(0, 1, 0.25)))
          print("All obs.:")
          print(quantile(vardata[,j], probs = seq(0, 1, 0.25)))
        } # if quantvar         
        par(mfrow=c(1,1))
      } # for j
      if (!(projmethod=="none" | sum(clusteri)<minsize)){
# print(j)
# print("proj")
        dp <- discrproj(contdata,as.integer(clusteri),clnum=1,method=projmethod)
        rx <- range(dp$proj[clusteri,1])
  #    print(rx)
        rrx <- rx[2]-rx[1]
        ry <- range(dp$proj[clusteri,2])
  #    print(ry)
        rry <- ry[2]-ry[1]
        plot(dp$proj[,1:2],xlab=paste(projmethod,1),ylab=paste(projmethod,2),
             xlim=c(rx[1]-rangefactor*rrx,rx[2]+rangefactor*rrx),
             ylim=c(ry[1]-rangefactor*rry,ry[2]+rangefactor*rry),
             main= paste("Cluster ",i," projection by ",projmethod),
             col=1+clusteri,pch=1+18*clusteri)
      } # if projection
    } # for i
    if (ask)
      par(ask=FALSE)
  } # if clusterwise
#  print("varwisetables")
  varwisetables <- list()
  for (i in 1:p){
    if (i %in% catvar){
#       print("tables cat")
     qx <- quantile(vardata[,i],seq(1/catvarcats,1,1/catvarcats))
      vx <- rep(1,n)
      for (j in 2:catvarcats)
        vx[vardata[,i]>qx[j-1]] <- j
      varwisetables[[i]] <- table(clustering,vx,
                                  dnn=c("Cluster",
                                    paste("Categorised ",vnames[i])))
      varwisetables[[i]] <- addmargins(varwisetables[[i]],1)
      if (proportions)
        varwisetables[[i]] <- prop.table(varwisetables[[i]],1)
    }
    else{
#      print("tables !cat")
      varwisetables[[i]] <- table(clustering,vardata[,i],dnn=c("Cluster",
                                    vnames[i]))
      varwisetables[[i]] <- addmargins(varwisetables[[i]],1)
      if (proportions)
        varwisetables[[i]] <- prop.table(varwisetables[[i]],1)
    }
  }
  class(varwisetables) <- "varwisetables"
  attr(varwisetables,"prop") <- "counts"
  if (proportions)
  attr(varwisetables,"prop") <- "proportions"
  varwisetables
}

print.varwisetables <- function(x,digits=3,...){
  l <- length(x)
  td <- 10^digits
  for (i in 1:l){
    if(attr(x,"prop")=="counts")
      print(x[[i]])
    else
      print(round(x[[i]]*td)/td)
    cat(" \n")
  }
  invisible(x)
}
          
clucols <- function(i, seed=NULL){
  if (!is.null(seed))
    set.seed(seed)
  out <- sample(colours()[-1],i)
  out
}

clugrey <- function(i,max=0.9)
  grey(seq(0,max,by=max/(i-1)))

clusym <- c(sapply(c(1:9,0),toString),intToUtf8(c(97:122,65:90),multiple=TRUE))






                   
  
  
      
    
        
        
        
          
         
        
        
                              
    
    
    
    
    
      
      
    
