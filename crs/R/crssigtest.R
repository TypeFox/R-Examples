## This function performs an individual bootstrap test of significance
## for a crs regression model. If the index `i' is provided it
## computes the test for predictor `i', otherwise it computes the test
## for all predictors, one at a time.

## (C) Jeffrey S. Racine (2011)

## Two methods for imposing the null are implemented, a) a residual
## bootstrap and b) a reordering of the predictor being tested (in
## place) that breaks any systematic relationship between the
## predictor and outcome.

crssigtest <- function(model = NULL,
                       index = NULL,
                       boot.num = 399,
                       boot.type = c("residual","reorder"),
                       random.seed = 42,
                       boot = TRUE) {

  ## Save any existing seed prior to setting the seed for this routine
  ## (we set the seed to ensure the same outcome when the test is run
  ## on the same data, otherwise resampling error would result in
  ## different P values each time the test is run which could be
  ## unsettling for some users - this can naturally be overridden by
  ## setting a unique seed each time the test is run)

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed <- TRUE
  } else {
    exists.seed <- FALSE
  }

  set.seed(random.seed)

  ## Some basic checking for invalid (unintentional) use of the
  ## function

  if(is.null(model)) stop(" you must provide a crs model")
  if(is.null(index)) index <- 1:NCOL(model$xz)
  if(index < 1 || index > NCOL(model$xz)) stop(" you must provide a valid index")
  boot.type <- match.arg(boot.type)

  ## Some storage vectors/matrices

  df1.vec <- numeric(length(index))
  df2.vec <- numeric(length(index))    
  uss.vec <- numeric(length(index))
  rss.vec <- numeric(length(index))    
  F.vec <- numeric(length(index))  
  P.vec.boot <- numeric(length(index))
  P.vec.asy <- numeric(length(index))  
  F.boot <- numeric(length=boot.num)
  F.boot.mat <- matrix(NA,nrow=boot.num,ncol=length(index))

  for(ii in 1:length(index)) {

    ## Get information from model each time the test is run (we modify
    ## these below hence the need to reset when multiple tests are
    ## requested)

    model.degree <- model$degree
    model.lambda <- model$lambda
    model.segments <- model$segments
    model.basis <- model$basis
    model.knots <- model$knots

    ## Placeholders for null degree/lambda/segments

    degree.restricted <- model.degree
    lambda.restricted <- model.lambda
    segments.restricted <- model.segments

    ## Determine whether variable(s) being considered is of type
    ## numeric or type factor

    xz.numeric <- FALSE
    if(is.numeric(model$xz[,index[ii]])) xz.numeric <- TRUE

    ## Test whether degree (or bandwidth) for the variable being
    ## tested is zero (or one). First, the predictors can be in any
    ## order (numeric/factor) so we need to be careful to set the
    ## degree/lambda for the predictor being considered

    if(xz.numeric) {

      ## Predictor being tested is numeric, so we need to determine
      ## the `degree' index (i.e. which numeric predictor is being
      ## tested, and what is its index among all predictors in the
      ## matrix of all predictors model$xz)

      degree.index <- 0
      for(jj in 1:index[ii]) if(is.numeric(model$xz[,jj])) degree.index <- degree.index + 1

      ## If the degree is zero (i.e. cross-validation has determined a
      ## variable is `irrelevant'), allow the predictor to be included
      ## using a degree 1/segments 1 fit (latter not necessary but
      ## zero cost to so doing) to avoid degeneracy of the statistic
      ## (if we did not the restricted and unrestricted models would
      ## coincide)

      model.degree[degree.index] <- ifelse(model$degree[degree.index]==0,1,model$degree[degree.index])
      model.segments[degree.index] <- ifelse(model$degree[degree.index]==0,1,model$segments[degree.index])

      ## Set the null degree to zero and segments to one (latter not
      ## necessary but zero cost to so doing)

      degree.restricted[degree.index] <- 0
      segments.restricted[degree.index] <- 1

    } else {

      ## Predictor being tested is categorical, so need to determine
      ## the `lambda' index

      lambda.index <- 0
      for(jj in 1:index[ii]) if(!is.numeric(model$xz[,jj])) lambda.index <- lambda.index + 1

      ## If the bandwidth lambda is one (i.e. cross-validation has
      ## determined a variable is `irrelevant'), allow the model to be
      ## included using the `frequency' fit (bandwidth of basically
      ## zero)

      model.lambda[lambda.index] <- ifelse(isTRUE(all.equal(model$lambda[lambda.index],1)),.Machine$double.eps,model$lambda[lambda.index])

      ## Set the null bandwidth to one

      lambda.restricted[lambda.index] <- 1

    }

    ## Now compute the unrestricted and restricted models with the
    ## caveats outlined above for the unrestricted model

    model.unrestricted <- crs(xz=model$xz,
                              y=model$y,
                              cv="none",
                              degree=model.degree,
                              lambda=model.lambda,
                              segments=model.segments,
                              basis=model.basis,
                              knots=model.knots)

    model.restricted <- crs(xz=model$xz,
                            y=model$y,
                            cv="none",
                            degree=degree.restricted,
                            lambda=lambda.restricted,
                            segments=segments.restricted,
                            basis=model.basis,
                            knots=model.knots)

    ## Compute the pseudo F-value under the null (here we could use
    ## the unrestricted residuals from the model fed to this function
    ## i.e.  residuals(model)) while we compute df for the model
    ## augmented (i.e. if degree==0 set degree==1) to avoid degeneracy
    ## when a variable is automatically removed.

    ## Note - if we test a categorical variable df1 will be zero (the
    ## categorical predictor does not enter as a variable, rather it
    ## enters as a weight and will not be recorded in model$k, the
    ## number of parameters). In this case we approximate the number
    ## of restrictions imposed using the value 1.

    uss.vec[ii] <- sum(residuals(model.unrestricted)^2)
    rss.vec[ii] <- sum(residuals(model.restricted)^2)

    df1.vec[ii] <- max(1,round(sum(model.unrestricted$hatvalues))-round(sum(model.restricted$hatvalues)))
    df2.vec[ii] <- model$nobs-round(sum(model.unrestricted$hatvalues))

    ## Compute the statistic and save each in F.vec (we allow multiple
    ## tests to be computed with one call to this function).

    F.df <- df2.vec[ii]/df1.vec[ii]
  
    F.pseudo <- F.df*(rss.vec[ii]-uss.vec[ii])/uss.vec[ii]
    F.vec[ii] <- F.pseudo

    ## If boot == TRUE conduct the bootstrap, otherwise only conduct
    ## the asymptotic P-value

    if(boot) {

      ## Bootstrap the P-values for the test statistic (we also record
      ## the asymptotic P-values)
      
      if(boot.type=="reorder") xz.boot <- model$xz
      
      for(bb in 1:boot.num) {
        
        ## Bootstrap sample under the null and recompute the
        ## `restricted' and `unrestricted' models for the null sample
        
        if(boot.type=="reorder") {
          
          xz.boot[,index[ii]] <- sample(model$xz[,index[ii]],replace=T)
          
          model.unrestricted.boot <- crs(xz=xz.boot,
                                         y=model$y,
                                         cv="none",
                                         degree=model.degree,
                                         lambda=model.lambda,
                                         segments=model.segments,
                                         basis=model.basis,
                                         knots=model.knots)
          
          model.restricted.boot <- crs(xz=xz.boot,
                                       y=model$y,
                                       cv="none",
                                       degree=degree.restricted,
                                       lambda=lambda.restricted,
                                       segments=segments.restricted,
                                       basis=model.basis,
                                       knots=model.knots)
          
        } else {
          
          y.boot <- fitted(model.restricted) + as.numeric(scale(sample(residuals(model),replace=T),center=TRUE,scale=FALSE))
          
          model.unrestricted.boot <- crs(xz=model$xz,
                                         y=y.boot,
                                         cv="none",
                                         degree=model.degree,
                                         lambda=model.lambda,
                                         segments=model.segments,
                                         basis=model.basis,
                                         knots=model.knots)
          
          model.restricted.boot <- crs(xz=model$xz,
                                       y=y.boot,
                                       cv="none",
                                       degree=degree.restricted,
                                       lambda=lambda.restricted,
                                       segments=segments.restricted,
                                       basis=model.basis,
                                       knots=model.knots)
          
        }
        
        
        ## Recompute the pseudo F-value under the null
        
        uss.boot <- sum(residuals(model.unrestricted.boot)^2)
        
        rss.boot <- sum(residuals(model.restricted.boot)^2)
        
        F.boot[bb] <- F.df*(rss.boot-uss.boot)/uss.boot
        
      }
      
      F.boot.mat[,ii] <- F.boot
      
      P.vec.boot[ii] <- mean(ifelse(F.boot > F.pseudo, 1, 0))

    }
      
    ## Compute the asymptotic P value
      
    P.vec.asy[ii] <- pf(F.pseudo,df1=df1.vec[ii],df2=df2.vec[ii],lower.tail=FALSE)
      
  }

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  ## Return a list containing test results
  
  return(sigtest.crs(index=index,
                 P=P.vec.boot,
                 P.asy=P.vec.asy,
                 F=F.vec,
                 F.boot=F.boot.mat,
                 df1=df1.vec,
                 df2=df2.vec,
                 rss=rss.vec,
                 uss=uss.vec,
                 boot.num=boot.num,
                 boot.type=boot.type,
                 xnames=names(model$xz)))

} 


sigtest.crs <- function(index,
                    P,
                    P.asy,
                    F,
                    F.boot,
                    df1,
                    df2,
                    rss,
                    uss,
                    boot.num,
                    boot.type,
                    xnames){

  tsig <- list(index=index,
               P=P,
               P.asy=P.asy,
               F=F,
               F.boot=F.boot,
               df1=df1,
               df2=df2,
               rss=rss,
               uss=uss,
               boot.num=boot.num,
               boot.type=switch(boot.type,
                 "residual" = "Residual",
                 "reorder" = "Reorder"),
               xnames=xnames)
                   
  tsig$reject <- rep('', length(F))
  tsig$rejectNum <- rep(NA, length(F))

  tsig$reject[a <- (P < 0.1)] <- '.'
  tsig$rejectNum[a] <- 10

  tsig$reject[a <- (P < 0.05)] <- '*'
  tsig$rejectNum[a] <- 5

  tsig$reject[a <- (P < 0.01)] <- '**'
  tsig$rejectNum[a] <- 1

  tsig$reject[a <- (P < 0.001)] <- '***'
  tsig$rejectNum[a] <- 0.1

  class(tsig) = "sigtest.crs"

  return(tsig)
  
}

print.sigtest.crs <- function(x, ...){
  cat("\nRegression Spline Significance Test",
      "\nTest Type: ", x$boot.type," (",x$boot.num,
      " replications)",
      "\nPredictors tested for significance:\n",
      paste(paste(x$xnames[x$index]," (",x$index,")", sep=""), collapse=", "),"\n\n",
      sep="")
      

  maxNameLen <- max(nc <- nchar(nm <- x$xnames[x$index]))

  cat("\nSignificance Test Summary\n")
  cat("P Value:", paste("\n", nm, ' ', blank(maxNameLen-nc), format.pval(x$P),
                        " ", formatC(x$reject,width=-4,format="s"),
                        "(F = ", formatC(x$F,digits=4,format="fg"),
                        ", df1 = ", x$df1,                        
                        ", df2 = ", x$df2,
                        ", rss = ", formatC(x$rss,digits=6,format="fg"),
                        ", uss = ", formatC(x$uss,digits=6,format="fg"),")",sep=''))
  cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
}

summary.sigtest.crs <- function(object, ...) {
  print(object)
}
