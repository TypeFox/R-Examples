## make a null group for comparison
hzar.make.LLfunc.null <- function (obsData,
    model.LL=obsData$model.LL, LLrejectedModel = -1e+08) 
{
    model.req <- function(pVal) { return(all(pVal<1,pVal>0));}
    model.gen <- function(pVal) { res <- substitute(function(x) rep(p,length(x)), list(p=pVal));
                                  return(eval(res));}
    eval.clineLL <- model.LL
    myRejectionLL <- LLrejectedModel
    param.free.names="pVal";
    param.fixed=list();
    old.formals <- formals(model.gen)
    if (length(old.formals) != (length(param.free.names) + length(param.fixed))) {
        warning("The length of the method formals does not match the length of the parameters supplied.")
    }
    ttt.formals <- old.formals[param.free.names]
    names(ttt.formals) <- param.free.names
    new.formals <- c(ttt.formals, param.fixed)
    formals(model.req) <- new.formals
    formals(model.gen) <- new.formals
    llFunc <- function(theta) {
        if (!do.call(model.req, as.list(theta))) 
            return(myRejectionLL)
        model = do.call(model.gen, as.list(theta))
        result <- eval.clineLL(model)
        if (identical(is.finite(result), TRUE)) 
            return(result)
        return(myRejectionLL)
    }
    return(llFunc)
}

hzar.dataGroup.null <-function(obsData){
  
  clineLLfunc<-hzar.make.LLfunc.null(obsData)
  
  nPS <- 100;
  pSeries <- (1:(nPS-1))/(nPS)
  pExpected <- sum(obsData$frame$n*obsData$frame$obsFreq)/sum(obsData$frame$n)
  pSeries <- c(pSeries,pExpected);
  data.mcmc <- mcmc(data=data.frame(pVal=pSeries));
  return(hzar.make.dataGroup(data.mcmc=data.mcmc,
                             llFunc=clineLLfunc,
                             ML.cline=hzar.gen.cline(list(pVal=pExpected),clineLLfunc)
                             ));
}
