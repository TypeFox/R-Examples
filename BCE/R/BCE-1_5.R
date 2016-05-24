## Main developers: Karel Van den Meersche, Karline Soetaert
## 
#####################################################################
#####################################################################


rescaleRows <- function (A,                 # matrix or dataframe to be row-rescaled: rowSums(A[rescale])=1
                         columns=1:ncol(A)) # vector containing indices of the columns that should be included in the normalisation
  {
    if (is.null(columns)) return(A)
    if (nrow(A)==1)
      {
        R <- sum(A[,columns])
      } else {
        R <- rowSums(A[,columns])
      }
    A[R>0,columns]<-A[R>0,columns]/R[R>0]
    A
  } # rescaled dataframe

####################################################################################

rdirichlet2 <- function(alpha)          # input and output are matrices; each row is a point in a simplex
  {
    l <- length(alpha)
    n <- ncol(alpha)
    x <- matrix(rgamma(l, alpha),ncol=n)
    x/rowSums(x)
  }

## sum of logprobabilities of the rowvectors of x given the rowvectors of alfa; log = TRUE!!!
logddirichlet2 <- function(x,alpha) sum((alpha - 1) * log(x)) - sum(lgamma(alpha)) + sum(lgamma(rowSums(alpha)))


####################################################################################

lsei1 <- function(A,                     # search x for which min||Ax-B||
                 B,                     # 
                 E,                     # Ex=F
                 F,                     # Ex=F
                 G,                     # Gx>H
                 H)                     # Gx>H
  {
    ## require(quadprog,quietly=TRUE)

    ## dvec  <- t(A) %*% B
    ## Dmat  <- t(A) %*% A
    ## Amat  <- t(rbind(E,G))
    ## bvec  <- c(F,H)
    ## solve.QP(Dmat ,dvec, Amat , bvec, meq=1)$solution

    ## require(limSolve,quietly=TRUE)
    lsei(A,B,E,F,G,H)$X
  }


lsei2 <- function(Rat,                   # min(X%*%Rat-Dat)
                 Dat,                   # min(X%*%Rat-Dat)
                 sddat,                 # weighting
                 G=diag(1,nrow(Rat)),   # G*x>= h 
                 H=rep(0,nrow(Rat)),    # G*x>= h
                 E=matrix(1,1,nrow(Rat)), # E*x=F
                 F=1                      # E*x=F
                 )
  {
    X <- matrix(NA,nrow(Dat),nrow(Rat))
    for (i in 1:nrow(Dat))
      {
        select <- !is.na(Dat[i,])
        A <- t(Rat[,select])/sddat[i,select]
        B <- Dat[i,select]/sddat[i,select]
        X[i,] <- lsei1(A,B,E,F,G,H)
      }
    dimnames(X) <- list(rownames(Dat),rownames(Rat))
    return(X)
  } # lsei2


#########################################################################################################

panel.cor <- function(x, y, digits=2, prefix="", cex.cor,...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
  }

###########################################################################################

BCE <- function(
                ## parameters
                Rat,                        # initial ratio matrix
                Dat,                        # initial data matrix

                relsdRat = 0,               # relative standard deviation on ratio matrix: a number or a matrix
                abssdRat = 0,               # absolute standard deviation on ratio matrix: a number or a matrix
                minRat   = 0,               # minimum values of ratio matrix: a number or a matrix
                maxRat   = +Inf,            # maximum values of ratio matrix: a number or a matrix
                
                relsdDat = 0,               # relative standard deviation on data matrix: a number or a matrix
                abssdDat = 0,               # absolute standard deviation on data matrix: a number or a matrix

                tol      = 1e-4,            # minimum standard deviation for data matrix
                tolX     = 1e-4,            # minimum x values for MCMC initiation
                positive = 1:ncol(Rat),     # which columns contain strictly positive data; other columns are not rescaled, and can become negative)
                iter     = 100,             # number of iterations for MCMC
                outputlength = 1000,        # number of iterations kept in the output
                burninlength = 0,           # number of initial iterations to be removed from analysis
                jmpRat   = 0.01,            # jump length of ratio matrix (in normal space): a number, vector or matrix
                jmpX     = 0.01,            # jump lenth of composition matrix (in a simplex): a number, vector or matrix
                unif     = FALSE,           # do we take uniform distro's for ratio matrix? (as in chemtax)
                verbose  = TRUE,            # if TRUE, extra information is provided during the run of the function, such as extra warnings, elapsed time and expected time until the end of the MCMC 
                initRat  = Rat,             # ratio matrix used to start the markov chain: default the initial ratio matrix
                initX    = NULL,            # composition matrix used to start the markov chain: default the LSEI solution of Ax=B
                userProb = NULL,            # posterior probability for a given ratio matrix and composition matrix: should be a function with 2 arguments RAT and X, and as returned value a number giving the -log posterior probability of ratio matrix RAT and composition matrix X. Dependence of the probability on the data should be incorporated in the function. 
                confInt  = 2/3,             # confidence interval in output; because the distributions are not symmetrical, standard deviations are not a useful measure; instead, upper and lower boundaries of the given confidence interval are given. Default is 2/3 (equivalent to standard deviation), but a more or less stringent criterion can be used. 
                export   = FALSE,           # if true, a list of variables and plots are exported to the specified filename in a folder "out". If a valid path, the list is exported to the location indicated by the path. 
                file     = "BCE"            # objects are saved to this file. 
                )
  
                                        # the mcmc function BCE assesses probability distributions of
                                        # ratio matrix Rat giving biomarker composition of a
                                        # number of taxa
                                        # and a composition matrix x giving taxonomical composition of a number of
                                        # stations
                                        # with Dat = the biomarker composition of the samples.
                                        # the probability of outcome x %*% Rat ~= Dat is evaluated with a mcmc. 
  
  {
    input.list <- list(
                       Rat               = Rat,
                       Dat               = Dat,
                       relsdRat          = relsdRat,
                       abssdRat          = abssdRat,
                       relsdDat          = relsdDat,
                       abssdDat          = abssdDat,
                       minRat            = minRat,
                       maxRat            = maxRat,
                       tol               = tol,     
                       tolX              = tolX,           
                       positive          = positive,
                       userProb          = userProb,
                       unif              = unif,
                       verbose           = verbose,
                       jmpX              = jmpX,
                       jmpRat            = jmpRat,
                       initRat           = initRat,
                       initX             = initX,
                       confInt           = confInt
                       )

    init.list <- init(input.list)

    with(init.list,{

      ##==========================================================##
      ## initial posterior probability of x, Rat and x%*%Rat = y given data
      ##==========================================================##

      rat2 <- rat1
      x2 <- x1
      
      logp1 <- logProbabilityBCE(rat1,x1,init.list)

      ##=======================================##
      ## initialise mcmc objects
      ##=======================================##
      if (burninlength>=iter) stop("The burninlength
defines the number of iterations to be removed from the analysis.
This burninlength cannot be smaller than iter, the total number of iterations!")
      ou <- ceiling((iter-burninlength)/outputlength)
      iter <- iter-(iter-burninlength)%%ou
      outputlength <- (iter-burninlength)%/%ou
      ou1 <- burninlength+ou
      i1 <- 1
      
      mcmc.Rat <- array(dim=c(nalg,npig,outputlength),dimnames=list(algnames,pignames,NULL))
      mcmc.X   <- array(dim=c(nst,nalg,outputlength),dimnames=list(stnames,algnames,NULL)) 
      mcmc.logp   <- vector(length=outputlength)         
      naccepted <- 0
      init.time <- proc.time()

      
      ##==========##
      ## mcmc loop
      ##==========##

      for (i in 1:iter)
        {
          ## new parameters
          x2 <- rdirichlet2(x1*(alfajmp-nalg)+1)
          rat2[] <- rnorm(lr,rat1,jmpRat.matrix)

          ## new posterior probability p2 (-log)
          logp2  <-  logProbabilityBCE(rat2,x2,init.list)
          r <- exp(logp2-logp1 + logddirichlet2(x1,x2*(alfajmp-nalg)+1) - logddirichlet2(x2,x1*(alfajmp-nalg)+1))
          
          ## METROPOLIS algorithm: select the new point? or stick to the old one? 
          if (r>=runif(1))
            {
              ## update x1, rat1, logp1, naccepted
              x1 <- x2
              rat1 <- rat2
              logp1 <- logp2
              naccepted <- naccepted+1
            }
          
          ## update mcmc objects
          if (i==ou1)
            {
              mcmc.Rat[,,i1] <- rat1
              mcmc.X[,,i1] <- x1
              mcmc.logp[i1] <- logp1
              i1 <- i1+1
              ou1 <- ou1+ou
              
              ## give some process feedback
              if (i %in% c(100,1000,1:10*10000,iter))
                {
                  if (verbose)
                    {
                      present.time <- proc.time()-init.time
                      print(cat("runs:",i,"; elapsed time:",present.time[3],"s ; estimated time left:",present.time[3]*(iter-i)/i,"s ; speed:",i/present.time[3]," runs/s ")) 
                      flush.console()
                    }
                }
            }
          
        } # end mcmc loop

      if (dim(mcmc.X)[1]==1) {dim(mcmc.X) <- dim(mcmc.X)[2:3] ; rownames(mcmc.X) <- algnames}
      if (verbose) print(cat("number of accepted runs: ",naccepted," out of ",iter," (",100*naccepted/iter,"%) ",sep=""))

      mcmc <- list(Rat=mcmc.Rat,X=mcmc.X,logp=mcmc.logp,naccepted=naccepted)
      class(mcmc) <- c("bce","list")
      if (export) export(mcmc,file,input.list)
      return(mcmc)           # bce object: a list containing 4 elements:
                                        # - mcmc.Rat: array with dimension c(nrow(Rat),ncol(Rat),iter) containing the random walk values of the ratio matrix
                                        # - mcmc.X: array with dimension c(nrow(x),ncol(x),iter) containing the random walk values of the composition matrix
                                        # - mcmc.logp: vector with length iter containing the random walk values of the posterior probability
                                        # - naccepted: integer indicating the number of runs that were accepted
    })}





####################################################

init <- function(input.list)
  {
    with(input.list,{

      ## warnings and error messages
      if (ncol(Rat)!=ncol(Dat)) stop("ratio matrix and data matrix must have same number of columns")
2
      ##===============================================##
      ## initialisations
      ##===============================================##

      if (is.vector(Dat)) Dat <- t(Dat)                          # Dat has to be a matrix
      if (is.data.frame(Dat)) Dat <- as.matrix(Dat)
      if (is.data.frame(Rat)) Rat <- as.matrix(Rat)
      if (is.data.frame(initRat))   initRat <- as.matrix(initRat)
      if (is.data.frame(relsdDat)) relsdDat <- as.matrix(relsdDat)
      if (is.data.frame(relsdRat)) relsdRat <- as.matrix(relsdRat)
      if (is.data.frame(abssdDat)) abssdDat <- as.matrix(abssdDat)
      if (is.data.frame(abssdRat)) abssdRat <- as.matrix(abssdRat)

      ## useful numbers
      nalg <- nrow(Rat)    ;    algnames   <- rownames(Rat)       # number & names of taxonomic groups
      nst  <- nrow(Dat)    ;    stnames    <- rownames(Dat)       # number & names of stations or samples
      npig <- ncol(Rat)    ;    pignames   <- colnames(Rat)       # number & names of biomarkers
      lx <- nst*nalg
      lr <- nalg*npig
      ld <- nst*npig
      
      ## standard deviations  for Rat and Dat

      if (length(relsdRat)==1) relsdRat <- rep(relsdRat,ncol(Rat))
      if (is.vector(relsdRat)&length(relsdRat)==ncol(Rat)) relsdRat <- t(matrix(relsdRat,nrow=ncol(Rat),ncol=nrow(Rat)))
      if (length(relsdRat)!=length(Rat)) stop("invalid dimensions of relsdRat")

      if (length(abssdRat)==1) abssdRat <- rep(abssdRat,ncol(Rat))
      if (is.vector(abssdRat)&length(abssdRat)==ncol(Rat)) abssdRat <- t(matrix(abssdRat,nrow=ncol(Rat),ncol=nrow(Rat)))
      if (length(abssdRat)!=length(Rat)) stop("invalid dimensions of abssdRat")

      sdrat=as.matrix(relsdRat*Rat+abssdRat)

      if (length(minRat)==1) minRat <- rep(minRat,ncol(Rat))
      if (is.vector(minRat)&length(minRat)==ncol(Rat)) minRat=t(matrix(minRat,nrow=ncol(Rat),ncol=nrow(Rat)))
      if (length(minRat)!=length(Rat)) stop("invalid dimensions of minRat")

      if (length(maxRat)==1) maxRat <- rep(maxRat,ncol(Rat))
      if (is.vector(maxRat)&length(maxRat)==ncol(Rat)) maxRat=t(matrix(maxRat,nrow=ncol(Rat),ncol=nrow(Rat)))
      if (length(maxRat)!=length(Rat)) stop("invalid dimensions of maxRat")
      
      if (length(relsdDat)==1) relsdDat <- rep(relsdDat,ncol(Dat))
      if (is.vector(relsdDat)&length(relsdDat)==ncol(Dat)) relsdDat <- t(matrix(relsdDat,nrow=ncol(Dat),ncol=nrow(Dat)))
      if (length(relsdDat)!=length(Dat)) stop("invalid dimensions of relsdDat")

      if (length(abssdDat)==1) abssdDat <- rep(abssdDat,ncol(Dat))
      if (is.vector(abssdDat)&length(abssdDat)==ncol(Dat)) abssdDat <- t(matrix(abssdDat,nrow=ncol(Dat),ncol=nrow(Dat)))
      if (length(abssdDat)!=length(Dat)) stop("invalid dimensions of abssdDat")
      
      sddat=as.matrix(relsdDat*Dat+abssdDat)
      
      if (any(sddat==0,na.rm=TRUE))
        {
          sddat[sddat==0] <- tol
          warning("Some elements in the data matrix have standard deviation = 0. They are set to a minimum value (tol)")
          if (verbose) flush.console()
        }
                                        # for elements of B that have sd=0, we want them to be changed very little when determining the
                                        # optimal posterior distribution. Instead of excluding them from analysis, which would be fairly
                                        # complex to implement, we give them a standard deviation tol.

      ## lamda and k for gamma distro Rat and Dat

      krat <- Rat^2/sdrat^2
      lrat <- Rat/sdrat^2
      krat[Rat==0&sdrat!=0] <- 1
      lrat[Rat==0&sdrat!=0&!is.na(Rat)] <- 1/sdrat[Rat==0&sdrat!=0&!is.na(Rat)]

      kdat <- Dat^2/sddat^2
      ldat <- Dat/sddat^2
      kdat[Dat==0] <- 1
      ldat[Dat==0&!is.na(Dat)] <- 1/sddat[Dat==0&!is.na(Dat)]


      select <- sdrat>0; if (unif) select <- Rat>0

      wholeranged <- !(1:npig%in%positive)
      select.pos <- select ; select.pos[,wholeranged] <- FALSE
      ind.pos <- which(select.pos)
      select.r <- select ; select.r[,positive] <- FALSE
      ind.r <- which(select.r)

      ##==========================================================##
      ## initialisation x with LSEI
      ##==========================================================##

      if (is.null(initRat)) rat1 <- Rat else rat1 <- initRat
      rat1[is.na(rat1)] <- 0
      if (is.null(initX))
        {
          x <- lsei2(rat1,Dat,sddat,H=rep(tolX,nalg))
          x1 <- x
        } else x1 <- x <- initX

      ##=========================================================##
      ## initialisation jump lengths
      ##=========================================================##
      alfajmp <- 1/(4*jmpX^2)-1
      
      if (length(jmpRat)==1)
        {
          jmpRat.matrix <- matrix(jmpRat,nrow=nalg,ncol=npig)
        } else {

          if (is.vector(jmpRat)&length(jmpRat)==npig) {
            jmpRat.matrix <- t(matrix(jmpRat,ncol=nalg,nrow=npig))
          } else {
            if (all(dim(jmpRat)==dim(Rat))) {
              jmpRat.matrix <- as.matrix(jmpRat)
            } else {
              stop("The jump length of the ratio matrix should be either a single value, or specified for each biomarker separately")
            }
          }
        }
      
      jmpRat.matrix[sdrat==0] <- 0             # only jump when standard deviation >0
      if (any(is.na(jmpRat.matrix))) stop("missing values in jump ratio matrix; please specify a valid jump ratio matrix.")

      
      return(list(
                  Rat               = Rat,
                  Dat               = Dat,
                  sdrat             = sdrat,
                  sddat             = sddat,
                  minRat            = minRat,
                  maxRat            = maxRat,
                  tol               = tol,     
                  tolX              = tolX,           
                  positive          = positive,
                  wholeranged       = wholeranged,
                  userProb          = userProb,
                  unif              = unif,
                  verbose      = verbose,
                  
                  krat         = krat,
                  lrat         = lrat,
                  kdat         = kdat,
                  ldat         = ldat,
                  algnames     = algnames,
                  stnames      = stnames,
                  pignames     = pignames,
                  nalg         = nalg,
                  nst          = nst,
                  npig         = npig,
                  lx           = lx,
                  lr           = lr,
                  ld           = ld,
                  x            = x,
                  ind.r        = ind.r,
                  ind.pos      = ind.pos,
                  x1           = x1,
                  rat1         = rat1,
                  alfajmp      = alfajmp,
                  jmpRat.matrix = jmpRat.matrix,
                  confInt       = confInt
                  ))
    })
  } # end initializations


#############################################################################################

logProbabilityBCE <- function(RAT,        # ratio matrix
                            X,          # composition matrix
                            init.list)  # list with variables, output of function init(input.list)
  {
    with(init.list,{

      if (!is.null(userProb)) logp <- log(userProb(RAT,X)) else {
        
        y <- X%*%RAT

        if (!is.null(positive))
          {
            dA.p <- dgamma(RAT[ind.pos],krat[ind.pos],lrat[ind.pos],log=TRUE)
            kdat <- y^2/sddat^2
            ldat <- y/sddat^2
            kdat[Dat==0] <- 1
            ldat[Dat==0&!is.na(Dat)] <- 1/y[Dat==0&!is.na(Dat)]
            dB.p <- dgamma(Dat[,positive],kdat[,positive],ldat[,positive],log=TRUE)
          } else {
            dA.p <- 0
            dB.p <- 0
          }
        
        if (any(wholeranged))
          {
            dA.r <- dnorm(RAT[ind.r],Rat[ind.r],sdrat[ind.r],log=TRUE)
            dB.r <- dnorm(y[,wholeranged],Dat[,wholeranged],sddat[,wholeranged],log=TRUE)
          } else {
            dA.r <- 0
            dB.r <- 0
          }
        
        dA.p[dA.p < -1e8] <- -1e8
        dA.p[dA.p > +1e8] <- +1e8
        dA.r[dA.r < -1e8] <- -1e8
        dA.r[dA.r > +1e8] <- +1e8
        
        drat <- sum(dA.p,dA.r,na.rm=TRUE)
        
        if (unif) drat <- 0
        if (any(RAT<minRat|RAT>maxRat,na.rm=TRUE)) drat <- -Inf
        
        
        dB.p[dB.p < -1e8] <- -1e8
        dB.p[dB.p > +1e8] <- +1e8
        dB.r[dB.r < -1e8] <- -1e8
        dB.r[dB.r > +1e8] <- +1e8

        ddat <- sum(dB.p,dB.r,na.rm=TRUE) /nst

        logp  <-  drat + ddat
      }
      
      return(logp)
    })
  }


##########################################################################################

export <- function(x,...) UseMethod("export")
export.bce <- function(x,             # a bce object, output of the function bce()
                       file="BCE",  # the bce object is written to this file
                       input.list=NULL, # a list of the arguments in bce() can be provided and saved as well.
                       ...)             # additional arguments
  {
    if (!is.null(attributes(x)$A_not_null)) {
      return("no export function is available for output of the function bce(); see ?BCE for the use of export.bce")
    } else {
      
      BCE <- x
      save(BCE,input.list,file=file)
      
      BCEsummary <- summary(BCE)
      
      with(c(BCE,BCEsummary),{
      write.csv(firstX,paste(file,"-firstX.csv",sep=""))
      write.csv(bestRat,paste(file,"-bestRat.csv",sep=""))
      write.csv(bestX,paste(file,"-bestX.csv",sep=""))
      write.csv(bestDat,paste(file,"-bestDat.csv",sep=""))
      write.csv(meanRat,paste(file,"-meanRat.csv",sep=""))
      write.csv(lbRat,paste(file,"-lbRat.csv",sep=""))
      write.csv(ubRat,paste(file,"-ubRat.csv",sep=""))
      write.csv(covRat,paste(file,"-covRat.csv",sep=""))
      write.csv(meanX,paste(file,"-meanX.csv",sep=""))
      write.csv(lbX,paste(file,"-lbX.csv",sep=""))
      write.csv(ubX,paste(file,"-ubX.csv",sep=""))
      write.csv(covX,paste(file,"-covX.csv",sep=""))

      png(paste(file,"%03d.png",sep=""),width=1903,height=1345,pointsize=10)
      par(mfrow=c(4,6))
      
      nalg <- nrow(Rat)
      npig <- ncol(Rat)
      nst <- nrow(bestDat)
      algnames <- rownames(Rat)
      pignames <- colnames(Rat)
      stnames <- rownames(bestDat)
      
      for (i in 1:nalg)
        {
          for (j in (1:npig)[meanRat[i,]!=0])
            {
              plot(Rat[i,j,],type="l",main=paste("trace of",algnames[i],pignames[j]),xlab="",ylab="")
              hist(Rat[i,j,],100,main=paste("histogram of",algnames[i],pignames[j]),xlab="")
            }
          a <- aperm(Rat)[,meanRat[i,]!=0,i]
          if (!is.null(dim(a))) pairs(a,upper.panel=panel.cor,pch=".")
        }
      if (is.matrix(X))
        {
          for(j in 1:nalg)
            {
              plot(X[j,],type="l",main=paste("trace of",algnames[j]),xlab="",ylab="")
              hist(X[j,],100,main=paste("histogram of",algnames[j]),xlab="")
            }
          a <- aperm(X[1:nalg,])
          if (!is.null(dim(a))) pairs(a,upper.panel=panel.cor,pch=".")
          
        } else {
          for (i in 1:nst)
            {
              for(j in 1:nalg)
                {
                  plot(X[i,j,],type="l",main=paste("trace of",algnames[j],stnames[i]),xlab="",ylab="")
                  hist(X[i,j,],100,main=paste("histogram of",algnames[j],stnames[i]),xlab="")
                }
              a <- aperm(X[i,1:nalg,])
              if (!is.null(dim(a))) pairs(a,upper.panel=panel.cor,pch=".")
            }
        }
      
      par(mfrow=c(1,1))
      barplot(t(bestX),legend.text=algnames)
      
      dev.off()
    })
    } 
  }                                    #end function export.bce()

#############################################################################################

plot.bce <- function(x,...)               # bce object
  if (!is.null(attributes(x)$A_not_null)) {
    NextMethod("modMCMC")
  } else {
    
    with(x,{
      
      nalg <- nrow(Rat)
      npig <- ncol(Rat)
      algnames <- rownames(Rat); if (is.null(algnames)) algnames <- paste("taxon",1:nalg)
      pignames <- colnames(Rat); if (is.null(pignames)) pignames <- paste("biomarker",1:npig)
    if (is.matrix(X)) {nst <- 1; stnames <- NULL} else { nst <- nrow(X); stnames <- rownames(X)}
      if (is.null(stnames)) stnames <- paste("sample",1:nst)
      
      oldpar <- par(no.readonly=TRUE)
      par(mfcol=c(nalg,5),mar=c(0,0,0,0),oma=c(0,0,1,0),ask=TRUE)
      
      for (j in 1:npig)
        {
          for (i in 1:nalg)
            {
              R <- Rat[i,j,]
              Rr <- range(R)
              if (all(Rr==0)) plot.new() else {
                ylim <- matrix(c(1,-.2,0,1.2),2)%*%Rr
                plot(R,type="l",xlab="",ylab="",xaxt="n",yaxt="n",ylim=ylim)
                text(0,ylim[2],paste(algnames[i],pignames[j]),adj=0)
              }
            }
        }
      mtext("ratio matrix traces",outer=TRUE)
      
      
      par(mfcol=c(nalg,5),mar=c(0,0,0,0),oma=c(0,0,1,0),ask=TRUE)
      for (i in 1:nst)
        {
          for(j in 1:nalg)
            {
            ifelse(is.matrix(X),x <- X[j,],x <- X[i,j,])
            xr <- range(x)
            ylim <- matrix(c(1,-.2,0,1.2),2)%*%xr
            plot(x,type="l",xlab="",ylab="",xaxt="n",yaxt="n",ylim=ylim)
            text(0,ylim[2],paste(stnames[i],algnames[j]),adj=0)
          }
        }
      par(oldpar)
    })
  }
