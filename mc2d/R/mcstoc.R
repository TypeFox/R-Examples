#<<BEGIN>>
mcstoc <- function(func=runif, type=c("V","U","VU","0"), ..., nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each", nsample="n", seed=NULL, rtrunc=FALSE, linf=-Inf, lsup=Inf, lhs=FALSE)
#TITLE Creates Stochastic mcnode Objects
#DESCRIPTION
# Creates a \code{\link{mcnode}} object using a random generating function.
#KEYWORDS methods
#INPUTS
#{func}<<A function providing random data or its name as character.>>
#[INPUTS]
#{type}<<The type of \samp{mcnode} to be built. By default, a \samp{"V"} node. see \code{\link{mcnode}} for details.>>
#{\dots}<<All other arguments but the size of the sample to be passed to \samp{func}. These arguments
#should be vectors or \samp{mcnode}s (arrays prohibited).>>
#{nsv}<<The number of simulations in the variability dimension.>>
#{nsu}<<The number of simulations in the uncertainty dimension.>>
#{nvariates}<<The number of variates of the output.>>
#{outm}<<The  output of the \samp{mcnode} for multivariates nodes. May be "each" (default)
#if an output should be provided for each variates considered independently, "none" for no output
#or a vector of functions (as a character string) that will be applied on the variates dimension
#before any output (ex: \samp{"mean"}, \samp{"median"}, \samp{c("min","max")}). Each function should return 1
#value when applied to 1 value (ex. do not use \samp{"range"}).
#Note that the \samp{outm} attribute may be changed further using the \code{\link{outm}} function.>>
#{nsample}<<The name of the parameter of the function giving the size of the vector.
#By default, \samp{n}, as in most of the random sampling distributions
# of the \samp{stats} library (with the exceptions of \samp{rhyper} and \samp{rwilcox} where \samp{nsample="nn"} should be used).>>
#{seed}<<The random seed used for the evaluation. If \samp{NULL} the \samp{seed} is unchanged.>>
#{rtrunc}<<Should the distribution be truncated? See \code{\link{rtrunc}}.>>
#{linf}<<If truncated: lower limit. May be a scalar, an array or a mcnode.>>
#{lsup}<<If truncated: upper limit. May be a scalar, an array or a mcnode. \samp{lsup} should be pairwise strictly greater then \samp{linf}>>
#{lhs}<<Should a Random Latin Hypercube Sampling be used? see \code{\link{lhs}}>>
#VALUE
#An \samp{mcnode} object.
#DETAILS
#Note that arguments after \dots must match exactly.
#
#Any function who accepts vectors/matrix as arguments may be used (notably: all current random generator of the \samp{stats} package).
#The arguments may be sent classically but it is STRONGLY recommended to use consistant \samp{mcnode}s
#if arguments should be recycled, since a very complex recycling is handled for \samp{mcnode} and not for vectors.
#The rules for compliance of \samp{mcnode} arguments are as following (see below for special functions):
#{type="V"}<<accepts \samp{"0" mcnode} of dimension \samp{(1 x 1 x nvariates)} or of dimension \samp{(1 x 1 x 1)} (recycled)
# and \samp{"V" mcnode} of dimension \samp{(nsv x 1 x nvariates)} or \samp{(nsv x 1 x 1)} (recycled).>>
#{type="U"}<<accepts \samp{"0" mcnode} of dimension \samp{(1 x 1 x nvariates)} or of dimension \samp{(1 x 1 x 1)} (recycled)
#and \samp{"U" mcnode} of dimension \samp{(1 x nsu x nvariates)} or of dimension \samp{(1 x nsu x 1)} (recycled).>>
#{type="VU"}<<accepts \samp{"0" mcnode} of dimension \samp{(1 x 1 x nvariates)} or of dimension \samp{(1 x 1 x 1)} (recycled),
#\samp{"V" mcnode} of dimension \samp{(nsv x 1 x nvariates)} (recycled classicaly) or \samp{(nsv x 1 x 1)} (recycled classically),
#\samp{"U" mcnode} of dimension \samp{(1 x nsu x nvariates)} (recycled by rows) or \samp{(1 x nsu x 1)}
#(recycled by row on the uncertainty dimension and classicaly on variates),
#\samp{"VU" mcnode} of dimension \samp{(nsv x nsu x nvariates)} or of dimension \samp{(nsv x nsu x 1)} (recycled).>>
#{type="0"}<<accepts \samp{"0" mcnode} of dimension \samp{(1 x 1 x nvariates)} or \samp{(1 x 1 x 1)}
#(recycled).>>
#
#Multivariate nodes and multivariate distributions:
#
#The number of variates should be provided (not guesses by the function).
#A multivariates node may be built using a univariate distribution and
#\samp{nvariates!=1}. See examples.
#
#\code{\link{rdirichlet}} needs for \samp{alpha} a vector or a multivariates nodes and returns a multivariate node.
#\code{\link{rmultinomial}} needs for \samp{size} and \samp{prob} vectors and/or multivariate nodes and return a univariate or a multivariate node.
#\code{\link{rmultinormal}} needs for \samp{mean} and \samp{sigma} vectors and/or multivariate nodes and return a multivariate node.
#\code{\link{rempiricalD}} needs for \samp{values} and \samp{prob} vectors and/or multivariate nodes and return a a univariate or a multivariate node.
#See examples.
#
#\samp{trunc=TRUE} is valid for univariates distributions only.
#The distribution will be truncated on \samp{(linf, lsup]}.
#The function 'func' should have a 'q' form (with first argument 'p') and a 'p' form, as
#all current random generator of the \samp{stats} library.
#Example : 'rnorm' (has a 'qnorm' and a 'pnorm' form), 'rbeta', 'rbinom', 'rgamma', ...
#
#If \samp{lhs=TRUE}, a Random Hypercube Sampling will be used on \samp{nsv} and \samp{nsu}
#The function 'func' should have a 'q' form (with argument 'p').
#\samp{lhs=TRUE} is thus not allowed on multivariates distributions.
#
#
#SEE ALSO
#\code{\link{mcnode}} for a description of \samp{mcnode} object, methods and functions on \samp{mcnode} objects.</>
#\code{\link{Ops.mcnode}} for operations on \samp{mcnode} objects.
#\code{\link{rtrunc}} for important warnings on the use of the \samp{trunc} option.
#EXAMPLE
#Oldnvar <- ndvar()
#Oldnunc <- ndunc()
#ndvar(5)
#ndunc(4)
#
### compatibility with mcdata as arguments
#x0 <- mcstoc(runif,type="0")
#xV <- mcstoc(runif,type="V")
#xU <- mcstoc(runif,type="U")
#xVU <- mcstoc(runif,type="VU")
#
### "0" accepts mcdata "0"
#mcstoc(runif,type="0",min=-10,max=x0)
#
### "V" accepts "0" mcdata and "V" mcdata
#mcstoc(rnorm,type="V",mean=x0,sd=xV)
#
### "U" accepts "0" mcdata and "U" mcdata
#mcstoc(rnorm,type="U",mean=x0,sd=xU)
#
### "VU" accepts "0" mcdata, "U" mcdata
### "V" mcdata and "U" mcdata with correct recycling
#mcstoc(rnorm,type="VU",mean=x0,sd=xVU)
#mcstoc(rnorm,type="VU",mean=xV,sd=xU)
#
### any function giving a set (vector/matrix) of value of length 'size' works
#f <- function(popi) 1:popi
#mcstoc(f,type="V",nsample="popi")
#
###Multivariates
#
#ndvar(2)
#ndunc(5)
###Build a multivariate node with univariate distribution
#mcstoc(rnorm,"0",nvariates=3)
#mcstoc(rnorm,"V",nvariates=3)
#mcstoc(rnorm,"U",nvariates=3)
#mcstoc(rnorm,"VU",nvariates=3)
#
###Build a multivariate node with multivariates distribution
#alpha <- mcdata(c(1,1000,10,100,100,10,1000,1),"V",nvariates=4)
#(p <- mcstoc(rdirichlet,"V",alpha=alpha,nvariates=4))
#mcstoc(rmultinomial,"VU",size=10,p,nvariates=4)
#
###Build a univariates node with "multivariates" distribution
#size <- mcdata(c(1:5),"U")
#mcstoc(rmultinomial,"VU",size,p,nvariates=1) #since a multinomial return one value
#
###Build a multivariates node with "multivariates" distribution
#mcstoc(rmultinomial,"VU",size,p,nvariates=4) #sent 4 times to fill the array
#
###Use of rempiricalD with nodes
###A bootstrap
#ndunc(5)
#ndvar(5)
#dataset <- c(1:9)
#(b <- mcstoc(rempiricalD,"U",nvariates=9,values=dataset))
#unclass(b)
###Then we build a VU node by sampling in each set of bootstrap
##in the uncertainty dimensions
#(node <- mcstoc(rempiricalD,"VU",values=b))
#unclass(node)
#
### truncated
#ndvar(2)
#ndunc(5)
#linf <- mcdata(-1:3,"U")
#x <- mcstoc(rnorm,"VU",rtrunc=TRUE,linf=linf)
#unclass(round(x))
## lhs and truncated with linf as mcnode
#linf <- mcdata(1:5,"U")
#mcstoc(rnorm,"VU",nsv=100,rtrunc=TRUE,linf=linf,lhs=TRUE)
#
#ndvar(Oldnvar)
#ndunc(Oldnunc)

#CREATED 08-01-25
#--------------------------------------------
{

    func <- match.fun(func)
    if(!is.null(seed)) set.seed(seed)

    if(!is.character(outm)  || (outm != "none" && outm != "each" && !all(sapply(outm,exists,mode="function"))))
      stop("outm should be 'none','each' or a vector of name(s) of valid function(s)")

    type <- match.arg(type)
    argsd <- list(...)
    dimf <- switch(type, "V"=c(nsv,1,nvariates),"U"=c(1,nsu,nvariates),"VU"=c(nsv,nsu,nvariates),"0"=c(1,1,nvariates))
    nsv <- dimf[1]
    nsu <- dimf[2]
    nva <- dimf[3]

    if(rtrunc) argsd <- c(argsd,list(linf=linf),list(lsup=lsup))          # launch linf and lsup in the process

    largsd <- length(argsd)

#### A function to deal mcnodes (including linf and lsup) as arguments

    LAFUNC <- function(argsd,typethismc){
          if(!is.null(typethismc)){                                    #mcnode as arguments

            if(!(type=="VU" || typethismc=="0" || typethismc==type)) stop("Incompatible type of nodes") # incompatible node

            dimm <- dim(argsd)
            if ((typethismc == "V" && dimm[1] != nsv) || 
                (typethismc == "U" && dimm[2] != nsu) || 
                (typethismc == "VU" && (dimm[1] != nsv || dimm[2] != nsu))) 
                stop("Nodes of incompatible dimensions")
 # incompatible dimension

            if(maxdim3 > 1){  #at least one multivariate node as parameter, need recycling on the third dimension
              if(typethismc=="U") argsd <- apply(argsd, 3, matrix, nrow=maxdim1, ncol=maxdim2, byrow=TRUE)  # recycling U as matrix (maxdim1*maxdim2) x nvariates
                else {
                  if(maxdim1 ==1 && maxdim2 ==1) argsd <- matrix(argsd, nrow=1) # Very special case to be added 
                  else argsd <- apply(argsd, 3, matrix, nrow = maxdim1, ncol = maxdim2)
                  }        # recycling 0, V, VU as matrix (maxdim1*maxdim2) x nvariates
            }                           

            else { dim(argsd) <- NULL    # as vector
                   if(typethismc == "U" && maxdim1!=1) argsd <- rep(argsd, each = maxdim1)     #recycling U as vector nsv*nsu
                    }
          }                                                           
          else if(is.array(argsd)) stop("Array prohibited in mcstoc as parameter. Use an mcnode instead")
      return(unclass(argsd))
    }

#### 

      typemc <- lapply(argsd, attr, which = "type")
      yamc <-   !is.null(unlist(typemc))                                        # At least one mcnode
      if(yamc){
            # evaluate the minimal common to build the minimal recycling level...
            maxdim1 <- unlist(lapply(argsd, function(x) dim(x)[1]))
            maxdim1 <- ifelse(is.null(maxdim1), 1, max(maxdim1))
            maxdim2 <- unlist(lapply(argsd, function(x) dim(x)[2]))
            maxdim2 <- ifelse(is.null(maxdim2), 1, max(maxdim2))
            maxdim3 <- unlist(lapply(argsd, function(x) dim(x)[3]))
            maxdim3 <- ifelse(is.null(maxdim3), 1, max(maxdim3))
      }

    if(largsd != 0){
      argsd <- mapply(LAFUNC, argsd, typemc, SIMPLIFY=FALSE)
      #print(argsd)
      }
#####################################
# If lhs or rtrunc, redefine the function to draw random variables 
#########
  #keep a copy of original function
  funcorigin <- func
  
  if(lhs || rtrunc){                                             #define good function for the random sampling
    distr <- as.character(match.call()$func)                     #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                              #remove the r
    qfun <- paste("q",distr,sep="")                              #define "qfunc"

  if(rtrunc){
       pfun <- paste("p",distr,sep="")                            #define pfunc

       func <- function(...){
          argsd <- list(...)
          nnfin <- argsd[[nsample]]
          linf <- if(length(argsd$linf) <= nnfin) as.vector(argsd$linf) else rep(argsd$linf, length.out=nnfin)      # solve a problem when linf was multivariate
          lsup <- if(length(argsd$lsup) <= nnfin) as.vector(argsd$lsup) else rep(argsd$lsup, length.out=nnfin)      # solve a problem when linf was multivariate
          lmax <- max(length(linf),length(lsup))
          if(any(rep(linf, length.out = lmax) >= rep(lsup, length.out = lmax))) stop("linf should be < lsup")  #recycle vectors
          argsd$linf <- argsd$lsup <- argsd[[nsample]] <- NULL
          #find the p of the limit
          pinf <- as.vector(do.call(pfun,c(list(q=linf),argsd),quote=TRUE))
          psup <- as.vector(do.call(pfun,c(list(q=lsup),argsd),quote=TRUE))
          #sample uniformely between the limits
          if(!lhs) lesp <- runif(nnfin,min=pinf,max=psup)
          else     lesp <- lhs(distr="runif", nsv=dimf[1], nsu=dimf[2], nvariates=dimf[3], min=pinf, max=psup) 
          #get the q
          data <- (do.call(qfun,c(list(p=lesp),argsd)))[1:nnfin]
          data[pinf==0 & data > lsup] <- NaN          #ex: rtrunc("lnorm",10,linf=-2,lsup=-1)
          data[psup==1 & data < linf] <- NaN          #ex: rtrunc("unif",10,linf=2,lsup=4,max=1)
          data[is.na(linf) | is.na(lsup)] <- NaN      #ex: rtrunc("norm",10,sd=-2)

          #Two tests for extreme situations. None Catch all possibilities. THe error is first to avoid the warning
          if(any(data <= linf | data > lsup, na.rm=TRUE)) stop("Error in rtrunc: some values are not in the expected range (maybe due to rounding errors)")
          if(isTRUE(all.equal(pinf,1)) | isTRUE(all.equal(psup,0)) ) warning("Warning: check the results from rtrunc. It may have reached rounding errors")
          return(data)
       }#end redefinition func
    } 
    else func <- function(...) {                      # LHS only
          argsd <- list(...)
          argsd[[nsample]] <- NULL
          lesp <- lhs(distr="runif", nsv=dimf[1], nsu=dimf[2], nvariates=dimf[3], min=0, max=1)
          return(do.call(qfun,c(list(p=lesp),argsd)))}
  }
  
  # do a try to test the length if nvariates != 1
  if(nvariates != 1){                                                 
    if(largsd != 0) argsdtest <- mapply(function(x,typemc){
                                          if(is.null(typemc)) return(unclass(x))
                                          if(is.matrix(x))    return(x[1,,drop=FALSE])     # mc (they have been unclassed)
                                          return(x[1])}, 
                                        argsd, typemc, SIMPLIFY=FALSE)
      else argsdtest <- vector(mode="list",length=0)

    argsdtest[[nsample]] <- 1
    if(rtrunc) argsdtest$linf <- argsdtest$lsup <- NULL
    dimf <- c(1,1,1)
    data <- do.call(funcorigin,argsdtest,quote=TRUE)
    l <- length(data)
    if(l == nvariates) {
        if(rtrunc | lhs) stop("mcstoc does not handle rtrunc and lhs for multivariate distributions")
          dimf <- c(nsv,nsu,1)}      # If it returns a vector 
     else if(l == 1) dimf <- c(nsv,nsu,nvariates)      # if it returns a number
          else stop("the function should return a vector of size 1 or nvariates if",nsample,"=1")
        argsd[[nsample]] <- prod(dimf)
        data <- do.call(func, argsd, quote = TRUE)
        
       #Post Production, multivariate
       if (yamc){
          if(l==1){ # univariate distribution
                if(maxdim1 == 1 && maxdim2 == 1) data <- aperm(array(data, dim = c(nvariates, nsv, nsu)), c(2, 3, 1))
           else if(maxdim1 == 1 && nsv!=1)       data <- aperm(array(data, dim = c(nsu, nvariates, nsv)), c(3, 1, 2))
           else if(maxdim2 == 1 && nsu!=1)       data <- aperm(array(data, dim = c(nsv, nvariates, nsu)), c(1, 3, 2))
           else data <- array(data, dim = c(nsv, nsu, nvariates))
          }
          else {     # l != 1 : multivariate 
                if(maxdim1 == 1 && nsv != 1)     data <- aperm(array(data, dim = c(nsu, nsv, nvariates)), c(2, 1, 3))
           else data <- array(data, dim = c(nsv, nsu, nvariates))
           }
        }
        else data <- array(data, dim = c(nsv, nsu, nvariates))
      } #end multivariates
   else{  # univariate
        argsd[[nsample]] <- prod(dimf)
        data <- do.call(func, argsd, quote = TRUE)

        if (yamc && maxdim1 == 1 && nsv != 1) data <- aperm(array(data, dim = c(nsu, nsv, nvariates)), c(2, 1, 3))
        else data <- array(data, dim = c(nsv, nsu, nvariates))
   } 
      
      class(data) <- "mcnode"
      attr(data,"type") <- type
      attr(data,"outm") <- outm
      return(data)
    }
