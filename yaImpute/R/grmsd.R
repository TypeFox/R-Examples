grmsd <- 
function (...,ancillaryData=NULL,vars=NULL,wts=NULL,rtnVectors=FALSE)
{
  if (missing(...)) stop ("... required")

  args <- list(...)  
  argLabs <- as.list(substitute(list(...)))[-1]
  names(args) <- if (is.null(names(argLabs))) unlist(argLabs) else
  {
    fixNames <- names(argLabs) == ""
    names(argLabs)[fixNames] <- argLabs[fixNames]
    names(argLabs)
  }
  
  okClasses <- c("yai","impute.yai","data.frame","matrix","lm")
  if (!is.null(wts)) 
  {
    if (any(wts < 0) || sum(wts) <= 0) stop("wts must be positive and sum > 0")
  }

  mgd <- list()
  for (objName in names(args))
  {
    object <- args[[objName]]
    if (!inherits(object,okClasses))
    {
      warning("object ",objName," class is not one of ",paste(okClasses,collapse=", "))
      next
    }
    if (inherits(object,"yai")) object <- 
           impute.yai(object,ancillaryData=ancillaryData,vars=vars,observed=TRUE)               
    # try to allow "lm" objects. This code may fail as there are many
    # methods in R that inherit from "lm". 
    if (inherits(object,"lm")) 
    {
      pr <- predict(object)
      ob <- pr + resid(object)
      # only one column?
      if (is.null(dim(pr))) 
      {
        object <- cbind(pr,ob)
        colnames(object) = c(objName,paste0(objName,".o"))
      } else {
        colnames(ob) = paste0(colnames(ob),".o")
        object <- cbind(pr,ob)
      }
    }
    object <- na.omit(object)
    if (nrow(object) == 0) 
    {
      warning("argument ",objName," has no rows.")
      next
    }    
    if (inherits(object,"matrix") & mode(object) != "numeric") 
    {
      warning("argument ",objName," must be numeric.")
      next
    }
    facts = if (inherits(object,"matrix")) FALSE else unlist(lapply(object,is.factor))
    if (any(facts)) 
    {
      if (all(facts)) 
      {
        warning("all variables are factors in ",objName) 
        next
      } else {
        nams <- names(facts)[facts]
        nams <- nams[-grep("[.]o$",nams)]
        warning("factor(s) have been removed from ",objName,": ",paste0(nams,collapse=", "))
        object <- object[,!facts,drop=FALSE]
      }
    }
    useVars <- if (is.null(vars)) colnames(object) else 
      {
        ov <- grep ("[.]o$",vars)
        ov <- if (length(ov) == 0) unique(c(vars,paste0(vars,".o"))) else vars
        intersect(ov,colnames(object))
      }
    if (length(useVars) == 0)
    {
      warning ("needed variables not found in ",objName)
      next
    }
    ov = useVars[grep ("[.]o$",useVars)]
    if (length(ov) == 0) 
    {
      warning ("no observed variables found in ",objName)
      next
    }
    pv <- unique(sub("[.]o$","",ov))
    pv <- intersect(pv,useVars)
    if (length(pv) == 0) 
    {
      warning("nothing to compute in ",objName)
      next
    }
 
    ob <- as.matrix(object[,ov,drop=FALSE])  
    pr <- as.matrix(object[,pv,drop=FALSE]) 
    qr <- qr(ob)

    uvars <- qr$pivot[1:qr$rank]
    if (length(uvars) < length(ov))
      warning("rank deficiency in ",objName," was addressed by removing: ",
             paste0(c(colnames(ob)[qr$pivot[(qr$rank+1):length(qr$pivot)]],
                      colnames(pr)[qr$pivot[(qr$rank+1):length(qr$pivot)]]),collapse=", "))

    p <- solve(chol(cov(ob[,uvars,drop=FALSE])))
    ob <- as.matrix(ob[,uvars]) %*% p
    pr <- as.matrix(pr[,uvars]) %*% p

    wt <- wts
    wt <- if (is.null(wt)) rep(1,ncol(pr)) else 
      {
        if (length(names(wt)) > 0) 
        {
          names(wt)  <- sub("[.]o$","",names(wt))
          wt <- na.omit(wt[names(pr)])
        }
        if (length(wt) != ncol(pr)) 
        {
          warning ("weights do not match variables in ",objName," and were ignored.")
          wt <- rep(1,ncol(pr))
        }
        wt
      }
    wt <- wt/sum(wt)
    md <- apply((pr-ob),1,function (x,wt) sum((x^2)*wt), wt)
    mgd[[objName]] <- if (rtnVectors) sqrt(md) else sqrt(mean(md))
  }
  if (rtnVectors) 
  { 
    idx <- sort(unlist(lapply(mgd,function (x) sqrt(mean(x)))),index.return=TRUE)$ix
    mgd[idx]
  } else sort(unlist(mgd))
}
