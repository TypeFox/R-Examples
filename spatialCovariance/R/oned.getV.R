
if(!exists("aniso")) aniso <- 1

collapse.sum <- function(vect,positions)
      {
        res <- NULL
        pos <- 1
        for(count.i in 1:length(positions))
          {
            res <- c(res,sum(vect[pos:(pos+positions[count.i]-2)]))
            pos <- pos+positions[count.i]-1
          }
        res
      }

computeV <- function(info,class="ldt",params,rel.tol = .Machine$double.eps^0.25,abs.tol = rel.tol,cat.level=0, K=NULL)
  {
    getVTime <- proc.time()
    
    n <- info$nrows*info$ncols
    area <- info$rowwidth*info$colwidth

    ## density for distance between rectangles
    ptm <- proc.time()
    
    if(class != "ldt" && class != "power")
      {
        KandCov <- defineK(class,K)
        if(class != "misc") K <- KandCov$K  ## covariance function
        cov.f <- KandCov$cov.f  ## product of covariance function and density
      }

    if(cat.level>=1) cat("loading in functions",(proc.time()-ptm)[3],"seconds\n")
    
    ## this is what needs to be computed
    ptm <- proc.time()
    if(class=="ldt") {
      results <- apply(info$indices,1,f.anal.ldt,info=info)
      message <- "evaluating ldt analytic results"
    } else if(class=="power") {
      results <- apply(info$indices,1,f.anal.power,h=params[1],info=info)
      message <- "evaluating power analytic results"
    } else {
      results <- apply(info$indices,1,f.NI,params=params,rel.tol=rel.tol,abs.tol=abs.tol,K=K,cov.f=cov.f,info=info)
      message <- "evaluating numerical integrals"
    }
    
    ## collapse results back down
    results <- collapse.sum(results,info$lengths)
    if(cat.level>=1) cat(message,(proc.time()-ptm)[3],"seconds\n")

    ptm <- proc.time()
    results <- rep(results,info$rowReps)
    V <- matrix(0,n,n)
    
    ## now stick the results into V
    V[info$locations] <- results
    V[info$locations[,c(2,1)]] <- results
    diag(V) <- V[1,1]
    if(cat.level>=1) cat("Inserting values in V",(proc.time()-ptm)[3],"seconds\n")
    if(cat.level) cat("computing V takes",(proc.time()-getVTime)[3],"seconds\n")
    V
  }

if(FALSE) {
  ## compute info when we have at least one row or one column
  ## code fails when we only have one plot
  if(!exists("info") && nrows+ncols>2) {
    getV.prec <- getV.precompute(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level)
    info <- getV.prec$info
    rm(getV.prec)
  }
  
  ## this is a hack to get V when we have one row and one column
  ## this is hardly ever used, but just for a complete setup we have included it here
  ## it is used when checking the accuracy and timing of evaluating the power jobs
  
  if(!exists("info")) {
    getV.prec <- getV.precompute(nrows+1,ncols+1,rowwidth,colwidth,rowsep,colsep,cat.level)
    info <- getV.prec$info
    source("11case.R",local=T)
    rm(getV.prec)
  }
  
  if(info$aniso != aniso) {
    getV.prec.update <- getV.precompute.update(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level)
    info <- getV.prec.update$info
    rm(getV.prec.update)
  }
  
  if(exists("getV.precompute")) rm(getV.precompute)
  
  getVresult <- getV(info,class,rel.tol,abs.tol,unitarea,cat.level,params)
  V <- getVresult$V
  area <- getVresult$area
  rm(getV,getVresult)
}
