#### function for testing for the presence of groupwide signal in the multivariate model
#### four options are possible: exact and approximate likelihood ratio tests, F test and marginal screening prototype tests
#### user has to specify a sigma value (error standard deviation)
#### user has choice to specify intercept or not (mu parameter) -- if not, mu is estimated and tests performed to account for this
#### user may also specify selected columns. this is done in a list of length equal to the number of groups
#### if no columns are specified for a given group, they are estimated using either the lasso or the marginal screening procedure, depending on the type specified
####  input:
####      - x = Matrix of all predictors
####      - y = Response VECTOR. Same number of rows as x
####      - groups = Vector of length equal to the number of columns in x. Labels the group membership of the corresponding columns of x.
####      - test.group = Label of group (as in 'groups') the nullity of which is to be tested.
####      - type = Type of test to perform. Choices are "ELR", "ALR", "F" or "MS"
####      - selected.col = Vector of indices into columns of x that represent pre-selected columns (assumed to be done unsupervised). Default is NULL.
####      - lambda = Regularisation parameter used for the selection of group columns using lasso for the unspecified groups
####      - mu, sigma = Intercept and error standard deviation. Sigma default is 1 and cannot be NULL. Mu default is NULL, in which case it is estimated during the testing process
####      - hr.iter, hr.burn.in = Number of hit-and-run samples and burn-in samples required when generating the reference distribution
####      - verbose, tol = Extra parameters
prototest.multivariate <-
function(x, y, groups, test.group, type=c("ELR", "ALR", "F", "MS"), selected.col=NULL, lambda, mu=NULL, sigma=1, hr.iter=50000, hr.burn.in=5000, verbose=FALSE, tol=10^-8){
  if (!is.matrix (x)){stop("x should be a matrix")}
  ### sort out the data
  y = as.matrix (y)
  x = apply (x, 2, function(col){col-mean(col)})
  x = apply (x, 2, function(col){col/sqrt(sum(col^2))})
  
  #######################
  ##### INPUT CHECK #####
  #######################
  
  # dimensions
  
  if (length(y) != nrow (x)){stop("x and y should have the same number of rows")}
  n = length(y)
  p = ncol(x)
  
  # groups
  unique.groups = unique (groups)
  K = length (unique.groups)
  groups.of.selected = groups[selected.col]
  if (!(test.group %in% unique.groups)){stop("'test.group' should be contained in 'groups'")}
  test.group.index = which (unique.groups == test.group)
  
  # type
  type = type[1]
  if (!(type %in% c("ELR", "ALR", "F", "MS"))){stop("type should be one of ELR, ALR, F or MS")}
  
  # selected cols
  if (!is.null(selected.col)){
    selected.col = sort (unique(selected.col))
    if (length(selected.col) == 0){stop("must select non-zero number of columns")}
    if (length(selected.col) > p | min(selected.col) <= 0 | max(selected.col) > p){stop("invalid columns selected")}
    if ((length(selected.col) + !is.null(mu)) > n){stop("too many columns selected")}
    if (type == "MS" & any(table (groups.of.selected) > 1)){stop("must select maximum one column per group for the marginal screening prototype method")}
  }
  
  ##########################
  ##### SELECT COLUMNS #####
  ##########################
  ## check for NULL groups and select columns, generating constraint matrices
  if (verbose){print ("Finalising column selection...", quote=FALSE)}
  selection.obj = lapply (unique.groups, function(group){
    selected.in.group = selected.col[groups.of.selected == group]
    if (length(selected.in.group) == 0){ # none specified, so select
      if (verbose){print(paste("   Selecting columns in group ", group, sep=""), quote=FALSE)}
      if (type == "MS"){ # marginal screening
        Ab.obj = mc.selection.A.b (y=y, X=x[,groups==group, drop=FALSE], mu=mu)
      }else{
        Ab.obj = enet.selection.A.b (y=y, X=x[,groups==group, drop=FALSE], lambda=lambda, mu=mu)
      }
      return (list(selected.col=(1:p)[groups==group][Ab.obj$which.col], selected.group=rep(group, length(Ab.obj$which.col)), A=Ab.obj$A, b=Ab.obj$b))
    }else{ # specified
      return (list(selected.col=selected.in.group, selected.group=rep(group, length(selected.in.group)), A=NULL, b=NULL))
    }
  })
  if (length(selection.obj[[test.group.index]]$selected.col) == 0){ # check whether we have any selected columns in the interest group
    warning ("no columns selected in the test.group: returning p-value 1")
    out = list(ts=0, p.val=1, selected.col=NULL, y.hr=NULL)
    class (out) = "prototest"
    return (out)
  }

  ## combine constraint matrices
  if (verbose){print ("Finding constraint matrices...", quote=FALSE)}
  A = do.call (rbind, lapply(selection.obj, function(obj){obj$A}))
  b = do.call (c, lapply(selection.obj, function(obj){obj$b}))
  selected.col = do.call (c, lapply(selection.obj, function(obj){obj$selected.col}))
  groups.of.selected = do.call (c, lapply(selection.obj, function(obj){obj$selected.group}))
  
  ## check whether the selected columns outnumber the observations
  if ((length (selected.col) + !is.null(mu)) > n){
    stop (paste("number of selected columns (", length(selected.col) + 1, ") exceeds number of observations (", n, ")", sep=""))
  }
  
  
  ################
  ##### TEST #####
  ################
  
  if (is.null(A)){ ### if A is null, do non-selective test
    y.hr = NULL
    if (verbose){print ("Non-selective test.", quote=FALSE)}
    if (type == "F"){
      if (verbose){print("   Performing F test...", quote=FALSE)}
      F.test.obj = nonselective.multivariate.F.test (x=x[,selected.col, drop=FALSE], y=as.matrix(y), groups=groups.of.selected, test.group=test.group, mu=mu)
      ts = F.test.obj$ts
      p.val = F.test.obj$p.val
    }else if(type == "MS"){
      if (verbose){print("   Performing marginal screening prototype test...", quote=FALSE)}
      mc.test.obj = nonselective.mc.test (x=x[,selected.col,drop=FALSE], y=y, test.index=test.group.index, mu=mu, sigma=sigma)
      ts = mc.test.obj$ts
      p.val = mc.test.obj$p.val
    }else if (type == "ELR"){
      if (verbose){print("   Performing exact likelihood ratio test...", quote=FALSE)}
      ts = compute.lr.stat.multi (x=x[,selected.col, drop=FALSE], y=as.matrix(y), groups=groups[selected.col], test.group=test.group, mu=mu, sigma=sigma, verbose=verbose, tol=tol)
      p.val = pchisq(ts, df=1, lower.tail=FALSE)
    }else if (type == "ALR"){
      if (verbose){print("   Performing approximate likelihood ratio test...", quote=FALSE)}
      ts = compute.lr.stat.multi (x=x[,selected.col, drop=FALSE], y=as.matrix(y), groups=groups[selected.col], test.group=test.group, mu=mu, sigma=sigma, exact=FALSE, verbose=verbose, tol=tol)
      p.val = pchisq(ts, df=1, lower.tail=FALSE)
    }
  }else{
    ### if A is not null, do selective test
    ### F and MS can have both hit-and-run and non-sampling tests
    ### ELR and ALR can only be sampling tests
    if (verbose){print("Selective test.", quote=FALSE)}
    
    if (type == "F" | type == "MS"){ # F and MC test have same structure in multivariate case
      
      if (hr.iter == 0){ # we want a non-sampling test
        y.hr = NULL
        if (verbose){print("   Performing non-sampling F or MS test...")}
        QRS.obj = compute.QRS.vectors (y=y, X.E=x[,selected.col,drop=FALSE], groups=groups[selected.col], test.group=test.group, A=A, b=b, mu=mu, type="multivariate") # test stat and per line intervals
        F.trunc.interval = find.overall.truncation.interval(QRS.obj$Q, QRS.obj$R, QRS.obj$S, verbose=FALSE)/QRS.obj$c # overall interval
        ts = QRS.obj$F.stat
        p.val = compute.trunc.F.test.p.value (ts, F.trunc.interval, df1=QRS.obj$df1, df2=QRS.obj$df2) # p-value
      }else{ # hit-and-run
        if (verbose){print("   Performing hit-and-run F or MS test...")}
        y.hr = hit.and.run.samples.multivariate.model (x=x[,selected.col,drop=FALSE], y=y, groups=groups[selected.col], test.group=test.group, A=A, b=b, mu=mu, sigma=sigma, hr.iter=hr.iter, hr.burn.in=hr.burn.in)
        ts.hr = nonselective.multivariate.F.test (x=x[,selected.col,drop=FALSE], y=cbind(y, y.hr), groups=groups[selected.col], test.group=test.group, mu=mu)$ts
        ts = ts.hr[1]
        p.val = mean (ts.hr[-1] >= ts.hr[1])
      }
      
    }else if (type == "ELR"){
      
      if (hr.iter == 0){ # not implemented
        stop ("non-sampling selective tests not implemented for ELR and ALR")
      }else{ # hit-and-run
        if (verbose){print("   Performing hit-and-run ELR test...")}
        y.hr = hit.and.run.samples.multivariate.model (x=x[,selected.col,drop=FALSE], y=y, groups=groups[selected.col], test.group=test.group, A=A, b=b, mu=mu, sigma=sigma, hr.iter=hr.iter, hr.burn.in=hr.burn.in)
        ts.hr = compute.lr.stat.multi (x=x[,selected.col,drop=FALSE], y=cbind(y, y.hr), groups=groups[selected.col], test.group=test.group, mu=mu, sigma=sigma, verbose=verbose, tol=tol)
        ts = ts.hr[1]
        p.val = mean (ts.hr[-1] >= ts.hr[1])
      }
      
    }else if (type == "ALR"){
      if (hr.iter == 0){ # not implemented
        stop ("non-sampling selective tests not implemented for ELR and ALR")
      }else{ # hit-and-run
        if (verbose){print("   Performing hit-and-run ALR test...")}
        y.hr = hit.and.run.samples.multivariate.model (x=x[,selected.col,drop=FALSE], y=y, groups=groups[selected.col], test.group=test.group, A=A, b=b, mu=mu, sigma=sigma, hr.iter=hr.iter, hr.burn.in=hr.burn.in)
        ts.hr = compute.lr.stat.multi (x=x[,selected.col,drop=FALSE], y=cbind(y, y.hr), groups=groups[selected.col], test.group=test.group, mu=mu, sigma=sigma, exact=FALSE, verbose=verbose, tol=tol)
        ts = ts.hr[1]
        p.val = mean (ts.hr[-1] >= ts.hr[1])
      }
      
    }
  }
  
  out = list(ts=ts, p.val=p.val, selected.col=selected.col, y.hr=y.hr)
  class (out) = "prototest"
  return (out)
}
