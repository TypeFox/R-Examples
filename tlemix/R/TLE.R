###########################################################################################
###   version info removed, have a look at the svn info
###     TLE procedure
###
###             kStar   = k*- size of the initial random subsample;
###             kTrim   = trimming parameter: size of the C-steps random subsample;
###             data    = data  frame containing the x and y variables;
###     .................................................................................
###             nit     = number of iterations;
###             msglvl  = level of messages;
###             result  = restart/continuation information;
###             cit     = number of iteration in refinement steps;
###             test    = expected true loglikelihood of the model. stop procedure if reached;
###     .................................................................................
###             Density = density function of type - function(data,solution);
###             Estimate = specific estimation procedure interface: function(data,ind,...);
###             ...     = specific parameters for Estimate procedure.
###
#########################################################################################
###
###
### Changelog:
### 2008-12-16:  - added sensible default parameters
###              - changed number of iterations from 100 to 10
###
### 2010-12-02:  - set proposed default parameters kTrim and kStar as integer
###              - repaired executing message "TLE C-step converged"
###              - set nonzero message level for final print and improving the print
###
##########################################################################################
##### setClass TLE
##########################################################################################
setClass("TLE",
                representation(estimate="ANY",
                                                iterbest="numeric",
                                                it="numeric",
                                                maxloglik="numeric",
                                                indbest="integer",
                                                indout="integer",
                                                tleweights="numeric",
                                                tlelogliks="numeric",
                                                tleclusters="ANY",
                                                kTrim="numeric",
                                                kStar="numeric",
                                                mcomp="ANY",
                                                nobs="integer",
                                                stop="logical",
                                                call="call"))

##########################################################################################
##### setMethod show for TLE object
##########################################################################################
if(!isGeneric("show")) setGeneric("show")
setMethod("show", "TLE",
          function(object){
            cat("\nCall:", deparse(object@call,0.75*getOption("width")),
                sep="\n")
            cat("\nEstimation:\n")
            print(object@estimate)

          })

##########################################################################################
##### setClass summary.TLE - class for summary object
##########################################################################################
setClass("summary.TLE",
         representation(call="call",
                                estimate="ANY",
                        kTrim="numeric",
                        kStar="numeric",
                        maxloglik="numeric",
                        nobs="integer",
                        nout="integer"))

##########################################################################################
##### setMethod for show() of summary.TLE objects
##########################################################################################
setMethod("show", "summary.TLE",
          function(object){
            cat("\nCall:", deparse(object@call,0.75*getOption("width")),
                sep="\n")
            cat("\n")
            print(object@estimate)
            cat("\n")
            cat("kTrim:", object@kTrim, "   Number of Observations:", object@nobs, "   Number of Outliers:",object@nout,"\n")
            cat("\n")
          })

#######################################
### plot
#######################################
### In vielen Fallen mochte man auch plot-Methoden fur eine S4-Klasse erstellen.
### Dies funktioniert z.B. folgenderma_en:


##if(!isGeneric("plot")) setGeneric("plot")
setGeneric("plot")
setMethod("plot", signature("TLE","missing"),
        #function(x, y = NULL, ...){
        function(x,y=NULL,...) {
        plot(x@estimate)
})



##########################################################################################
##### setGeneric function TLE
##########################################################################################
setGeneric("TLE",
        function(formula,family,data, kStar=NULL,kTrim=NULL, nit=10, msglvl=0, result=NULL, cit=9, test=NULL, nc=1, Density, Estimate, ...)
        standardGeneric("TLE"))

### just to make sure that some S3 generics are available in S4
setGeneric("summary")



### #######################################################################################
### #######################################################################################
### #######################################################################################
### ## setMethod TLE
### #######################################################################################
### if called with existing TLE object (as result=), then restart estimate...
### #######################################################################################
setMethod("TLE",
          signature(formula="formula",family="character",data="ANY",kStar="ANY",
                    kTrim="ANY", nit="ANY", msglvl="ANY", cit="ANY",
                    result="ANY",test="ANY",nc="numeric", Density="function", Estimate="function"),
          function(formula,family,data, kStar=NULL,kTrim=NULL, nit=10, msglvl=0, result=NULL,
                   cit=9, test=NULL,nc=1, Density, Estimate, ...)
          {
print(proc.time())
            mycall = match.call()
            nobs = nrow(data)            # number of observations
            dimension = ncol(data)       # number of dimensions

            ## proposed default parameter
            k = trunc((nobs + nc * (dimension + 1) + 1)/2)

            if(is.null(kTrim)) {
              kTrim = k
            } else if(kTrim<1) {
                kTrim = round(kTrim*nobs)
            } else kTrim=trunc(kTrim)

            if(kTrim<nobs/2) {
              kTrim=nobs-kTrim
            } else if(kTrim>nobs) {
                kTrim=nobs
            }

            if(kTrim<k) {
              warning("kTrim should be greater than (nobs+nc*(dimension+1)+1)/2")
            }

            if(is.null(kStar)) {
              ## set kStar to the average of kTrim and g(p+1).
              ## round to nearest integer value.
              kStar = round((kTrim + (nc*(dimension+1)))/2)
            } else kStar = round(kStar)
            ## cat(paste("k* and ktrim chosen: ",kStar,kTrim))

            ## check restart object
            ## scan if current run is restart run
            if(      is.null(result) ||
               is.null(result@kTrim) ||
               result@kTrim!=kTrim ||
               is.null(result@nobs) ||
               result@nobs!=nobs ||
               is.null(result@maxloglik) ||
               is.null(result@estimate)   ) {
              if(msglvl>0) print("New estimate start")
              estimate  = NULL                       # current model estimate
              iterbest  = 0                          # iteration of the best estimate
              it        = 0                          # current iteration number
              maxloglik = NULL                       # maximal loglikelihood
              mcomp     = 1                          # initial number of clusters if any
            } else {                                 # no, it's a restart...
              if(msglvl>0) print("Restart estimate!")
              ## save previous estimate elements
              estimate  = result@estimate            # iteration of the best estimate
              iterbest  = result@iterbest
              it                = result@it          # current iteration number
              maxloglik         = result@maxloglik   # maximal loglikelihood
              indbest           = result@indbest     # selected best subsample
              indout            = result@indout      # outliers
              tleweights        = result@tleweights  # tle-weights
              tlelogliks        = result@tlelogliks  # tle-logliks
              tleclusters       = result@tleclusters # tle-clusters
              mcomp             = result@mcomp       # number of components
            }
            Qnew        = NULL                       # loglikelihood of selected data

            ## if model based clustering, we have to create a list
            ## this is some kind of workaround
            if(family=="mvtnormal") {
              var_name <- as.character(substitute(data))
              data <- list(data)
              names(data) <- var_name
            }

            ## #######################################################################################
            ## #        main TLE loop
            ## #######################################################################################
            oldit=it    # set oldit to it (0 if new estimate, else result@it)
            for(it in (1:nit+oldit)) # 'nit' new iteratins, (1+oldit):(nit+oldit)
              {
                ##      set the initial random subsample of rank k*:
                ind <-sort(sample(nobs,kStar))

                ##      control prints:
                if(msglvl==1) { if((it-1)%%80==0) cat("\n.") else cat(".") }
                        if(msglvl>1) cat("\n iteration:",it,"\n")


                ## Initial Step. Estimate using random subsample of size k*
                est     = Estimate(data,ind,nc=nc,model=formula,family=family, ...)
                if(is.null(est)) {
                  next() # if we don't get an estimate, go to the next iteration
                }

                ## get estimated densities for all data; update: added model and family
                wrk     = Density(data,est,model=formula,family=family)
                if(is.null(wrk)) {
                  next() # if we don't get an estimate, go to the next iteration
                }
                loglik  = log(wrk$lik)          # estimated Loglikelihood vector
                ## the difference between $c and $cc is that $cc is in matrix form whereas $c is just a vector
                ## ie $cc has one column for every component;
                cluster = wrk$c                 # estimated clusters
                ccomp   = wrk$cc                # estimated clusters
                ncomp   = dim(wrk$cc)[2]        # number of components

                qtmp    = loglik[ind]           # loglikelihood of the data points in the subsample
                if(any(is.na(qtmp))) next()     # incorrect estimate
                Qold    = sum(qtmp)             # LogLikelihood of selected data

                if(msglvl>1) cat("\n Initial sample of size   ",kStar,
                                 "\n Estimated LogLik=        ",Qold,
                                 "\n sample observations:     ",ind,
                                 "\n")
                ## ###################################################################################
                ## C-steps:     ##########################################################################
                for(cstep in 1:cit) {
                  ## if(msglvl>1)progress.bar((it-oldit-1)*9+cstep-1,nit*9+1,title="TLE progress")
                  ## Improvement:
                  ## Get kTrim best cases.
                  ## Replace all NA values with machine minimal value,
                  ## in order to obtain correct kTrim indeces:
                  ##    if(any(is.na(loglik))) loglik[is.na(loglik)]=.Machine$double.xmin

                  kSort = sort(loglik,index.return=T,decreasing=TRUE)$ix[1:kTrim]

                  ## Use appropriate estimation function interface
                  est   = Estimate(data,kSort,cluster=ccomp, nc=nc,model=formula,family=family,...) # added model
                  if(is.null(est)) break()
                  ##

                  wrk   = Density(data,est,model=formula,family=family) # get estimate densities for all data, added model and family
                  if(is.null(wrk)) break()
                  loglik        = log(wrk$lik)           # estimated Loglikelihood vector
                  cluster       = wrk$c                  # estimated clusters
                  ccomp = wrk$cc                         # estimated clusters
                  ncomp = dim(wrk$cc)[2]                 # number of components

                  qtmp  = loglik[kSort]
                  if(any(is.na(qtmp))) break()           # incorrect estimate
                  Qnew  = sum(qtmp)                      # LogLikelihood of selected data

                  if(is.na(Qnew)) break()                # If not all logliks are incorrect (Binary model!)

                  totlik        = sum(loglik,na.rm=TRUE) # LogLikelihood of all data
                  weights       = rep(0,nobs)            # set all weights to 0
                  weights[kSort] = 1                     # set weights of objects in the subsample to 1; TLE weights
                                                         #
                  if(msglvl>1) cat("\n C-step",cstep," of size ",length(kSort),
                                   "\n Estimates LogLik=                ",Qnew,
                                   "\n sample observations:       ",sort(kSort),
                                   "\n")

                  if( is.null(maxloglik) || Qnew > maxloglik ) {           # save new "best" attributes
                                        maxloglik       = Qnew             # LogLikelihood of the best estimate
                                        iterbest        = it               # iteration of the best subsample
                                        totloglik       = totlik           # LogLikelihood of all data
                                        tleweights      = weights          # TLE weights
                                        tlelogliks      = loglik           # TLE logliks
                                        tleclusters     = cluster          # TLE clusters
                                        indbest         = sort(kSort)      # ordered indexes of best subsample
                                        indout          = (1:nobs)[-kSort] # indices of outliers
                                        indout          = sort(indout)     # sorted indices outliers
                                        estimate        = est              # estimate & coefficients
                                        if(is.null(ncomp)) ncomp=1
                                        mcomp           = ncomp            # number of clusters/components

                                        if(msglvl>0)   cat("\n =========================================",
                                       "\n loglik updated!",maxloglik,iterbest,mcomp,
                                       if(is.null(test)) "" else paste("(",format(maxloglik-test,digits=3),")",sep=""),
                                       "\n =========================================",
                                       "\n")
                  }

                  if(length(ind)==kTrim && all(ind==sort(kSort))) {
                    if(msglvl>1) print("TLE C-step converged")
                    break ()
                  }
                  ind<-sort(kSort)
                  Qold = Qnew
                }
                ## ########################
                ## End of the Improvement loop (C-steps)
                ## ########################
                stop=!is.null(test)&&!is.null(Qnew)&&test<=Qnew
                if(stop) {
                  print("target OK")
                  break()
                }
              }
            ## End of the main  iteration loop
            ##
            ##  Final print:
            if(msglvl>0) {
              if(!is.null(estimate))
                {
                  li=length(indbest)
                  if(!is.null(maxloglik))
                    {
                      cat("\n best subsample found:",
                          "\n ==========================================================================",
                          "\n TLE parameters        ",
                          "\n   kTrim             = ",kTrim,
                          "\n   kStar             = ",kStar,
                          "\n   nobs              = ",nobs,
                          "\n TLE estimates         ",
                          "\n   loglik            = ",maxloglik,
                          "\n   num of components = ",mcomp,
                          "\n   iter of best      = ",iterbest,
                          "\n   total num.iter.   = ",it,
                          "\n   sorted subsample  = ",indbest[1:(li%/%4)],
                          "\n                     = ",indbest[(  li%/%4+1):(2*li%/%4)],
                          "\n                     = ",indbest[(2*li%/%4+1):(3*li%/%4)],
                          "\n                     = ",indbest[(3*li%/%4+1):li],
                          "\n   outliers          = ",indout,
                          "\n ==========================================================================",
                          "\n")
                    }
                } else
              { cat(    "\n No improvement \n")
              }
            }
            cat("\n")
            ## if(msglvl>1)progress.bar(1,0)
            if(!is.null(estimate))                                # there is an estimate
              {
                result = new("TLE", estimate=estimate,
                  iterbest=iterbest,it=it, maxloglik=maxloglik, indbest=indbest, indout=indout,
                  tleweights=tleweights, tlelogliks=tlelogliks, tleclusters=tleclusters,
                  kTrim=kTrim, kStar=kStar, mcomp=mcomp, nobs=nobs, stop=stop)
                ## result=list(estimate=estimate,
                ## iterbest=iterbest,it=it,maxloglik=maxloglik,indbest=indbest,indout=indout,tleweights=tleweights,
                ## kTrim=kTrim,kStar=kStar,mcomp=mcomp,nobs=nobs)
                ## attr(result,"stop")=stop
                result@call <- mycall
                result
              }
            else return(new("TLE"))
          })


### #################################################################################################
### ####### summary method for TLE object
### #################################################################################################
setMethod("summary", "TLE",
          function(object, ...) {
            z <- new("summary.TLE",
                     call = object@call,
                     estimate = object@estimate,
                     kTrim = object@kTrim,
                     kStar = object@kStar,
                     ## glik = object@maxloglik,
                     nobs = object@nobs,
                     nout = length(object@indout))
            z
          } )


if(!isGeneric("estimate"))
        setGeneric("estimate", function(object) standardGeneric("estimate"))

### # Set method for estimate() function used to return the estimate
setMethod("estimate", signature="TLE",
          function(object)
          {
                        object@estimate
          } )
