#' @export
#' @import stats
calibSI <- function(dat, n=NULL, family, conf.level=.95,B=2000, gridlo=4, gridhi=20, gridix=2, tol=.03){
  #make sure the cov prob is measured with at least some accuracy
  if(B<1000){
    message("B not large enough and taken to be 1000")
    B <- 1000
  }
  
  #Argument error checking
  if(!family%in% c("binomial", "poisson", "pois", "gaussian", "normal")){stop("family is not supported")}
  if(family=="binomial") {if (length(n)>1|length(dat)>1){stop("dat and n argument must be scalars for binomial family")}}
  if(class(conf.level)!="numeric" | 0>conf.level|1<conf.level) {stop("conf.level must be between 0 and 1")}
  if(gridhi<=gridlo  | gridhi-gridlo < tol){stop("search grid misspecified; please specify gridlo<gridlo+tol<gridhi")}
  
  #initialize initial grid search likelihood values (st.levels) and the variables where the 
  #estimated coverage will be stored(cov.st.levels)
  st.levels <- seq(gridlo, gridhi, gridix)
  cov.st.levels <- rep(0,length(st.levels))  
  
  #check some likelihood levels to figure out where the nominal coverage likely is
  #Binomail family calculations
  if(family=="binomial"){
    
    for(k in 1:length(st.levels)){
      b.int <- binLikSI(dat, n, st.levels[k], conf=T, B=5000)
      cov.st.levels[k] <- b.int$conf.equiv
    }
  }
  
  #Poisson family calculations
  if(family=="poisson"| family=="pois"){
    for(k in 1:length(st.levels)){
      b.int <- poisLikSI(dat, st.levels[k], conf=T, B=5000)
      cov.st.levels[k] <- b.int$conf.equiv
    }
  }
  
  #Gaussian family calculations
  if(family=="gaussian"|family=="normal"){
    out <- list(normLikSI(dat,  1/exp(-1/2*qchisq(conf.level, 1))), 1/exp(-1/2*qchisq(conf.level,1)))
    names(out) <- c("si", "support.level")
    return(out)
  }
  
  #Error check to make sure estimated coverages are in the right range to 
  #cover conf.level specified.
  if(sum(cov.st.levels>=conf.level)==0|sum(cov.st.levels<= conf.level)==0 | cov.st.levels[1] >= conf.level){
    print(cbind(st.levels, cov.st.levels))
    stop("Support levels attempted were out of confidence range")}
  
  #Generate appropriate output for binomial and pois families
  if(family=="binomial") {
    #Coverage vs likelihood level in Bin is a step function so need to find front edge of step
    getSmallbinLev <- function(lev, tar, steps) {
      low.end <- lev-steps
      high.end <- lev
      golden.ratio = 2/(sqrt(5) + 1)
      while(high.end-low.end > tol){
        check.val <- high.end-(golden.ratio*steps)
        cov.check <- binLikSI(dat, n, check.val, conf=T, B=B)
        if( cov.check$conf.equiv < tar){
          low.end <- check.val
          steps <- high.end-check.val
        }
        if(cov.check$conf.equiv >= tar){
          high.end <- check.val
          steps <- high.end-low.end
        }
      }
      if(cov.check$conf.equiv > tar) { return(check.val)}
      if(cov.check$conf.equiv <= tar) { return(high.end)}
    }
    level <- getSmallbinLev(min(st.levels[cov.st.levels >= conf.level]), conf.level, gridix)
    out <- list(binLikSI(dat, n, level), level, cbind(st.levels, cov.st.levels))
    names(out) <- c("si", "support.level", "init.grid")
    return(out)
  } 
  
  
  if(family=="poisson"|family=="pois") {
    getSmallpoisLev <- function(lev, tar, steps) {
      low.end <- lev-steps
      high.end <- lev
      golden.ratio = 2/(sqrt(5) + 1)
      while(high.end-low.end > tol){
        check.val <- high.end-(golden.ratio*steps)
        cov.check <- poisLikSI(dat, check.val, conf=T, B=B)
        if( cov.check$conf.equiv < tar){
          low.end <- check.val
          steps <- high.end-check.val
        }
        if(cov.check$conf.equiv >= tar){
          high.end <- check.val
          steps <- high.end-low.end
        }
      }
      if(cov.check$conf.equiv > tar) { return(check.val)}
      if(cov.check$conf.equiv <= tar) { return(high.end)}
    }
    level <- getSmallpoisLev(min(st.levels[cov.st.levels >= conf.level]), conf.level, gridix)
    out <- list(poisLikSI(dat, level), level, cbind(st.levels, cov.st.levels))
    names(out) <- c("si", "support.level", "init.grid")
    return(out)
  } else {print("whoops"); return(NULL)}
}