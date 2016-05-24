compute_initial_weights = function(y,x){
# The function returns initial weights.

    initial_weights = as.numeric(cor(y,x))^2
    initial_weights = initial_weights/(sum(initial_weights))
    return(initial_weights)
}


compute_scores = function(y,x,m,B,initial_weights=NULL){
# The function returns RSM scores.

    p = ncol(x)
    scores = numeric(p)
    ns = numeric(p)

    for(k in 1:B){
        submodel = sample(1:p,size=m,replace=FALSE,prob=initial_weights)
        lm1 = lm(y~x[,submodel])
        weights = as.numeric((summary(lm1)$coef[-1,3])^2)
        submodel = submodel[which(!is.na(coef(lm1)[-1]))] #UPDATE: deal with collinearity
        scores[submodel] =  scores[submodel] + weights
        ns[submodel] = ns[submodel] + 1
    }
    ns = ifelse(ns!=0,ns,1)
    scores = scores/ns

    return(scores)
}

compute_scores_parallel = function(y,x,m,B,initial_weights,nslaves){
#This function returns RSM scores computed using parallel function foreach %dopar%.
    p = ncol(x)
    requireNamespace("foreach",quietly = TRUE)
    taskCount=floor(B/nslaves)
    r1 = foreach(i = 1:(nslaves-1), .combine=cbind, .export=c("compute_scores_partial")) %dopar%{

        scores_ns = compute_scores_partial(y,x,m,taskCount,initial_weights)

        scores_ns
    }

    taskCount = B-taskCount*(nslaves-1)

    r2 = compute_scores_partial(y,x,m,taskCount,initial_weights)

    r3 = cbind(r1,r2)
    scores = apply(r3[1:p,],1,sum)
    ns = apply(r3[-(1:p),],1,sum)
    ns = ifelse(ns!=0,ns,1)
    scores = scores/ns

    return(scores)
}



compute_scores_partial = function(y,x,m,taskCount,initial_weights){
# The function returns RSM scores computed based on taskCount repetitions. Unlike the compute_scores this function returns concatenated vector of scores and ns.

    p = ncol(x)
    scores = numeric(p)
    ns = numeric(p)

    for(k in 1:taskCount){
        submodel = sample(1:p,size=m,replace=FALSE,prob=initial_weights)
        lm1 = lm(y~x[,submodel])
        weights = as.numeric((summary(lm1)$coef[-1,3])^2)
        submodel = submodel[which(!is.na(coef(lm1)[-1]))] #UPDATE: deal with collinearity
        scores[submodel] =  scores[submodel] + weights
        ns[submodel] = ns[submodel] + 1
    }

    return(c(scores,ns))
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
# The function checks if the given number is whole number.

    abs(x - round(x)) < tol
}

select_finalmodel_bic = function(y,x,order1,thrs,penalty){
# The function returns final model using Bayesian Information Criterion.

    n = length(y)
    lm0 = lm(y~1)
    beta0 = as.numeric(coef(lm0))
    rss0 = sum(lm0$residuals^2)
    bic0 = n*log(rss0/n)+1*penalty

    xo = x[,order1[1:thrs]]
    xo = cbind(1,xo)
    xo = as.matrix(xo)


#     qr1 = qr(xo)
#     Q = qr.Q(qr1)
#     R = qr.R(qr1)
#     R_inv = solve(R)
#     tQy = t(Q) %*% y

    qr1 = qr(xo)
    Q = qr.Q(qr1)
    R = qr.R(qr1)
    small_diag = which(abs(diag(R))<1e-7)
    if(length(small_diag)>0){
      #Deal with collinearities, coefficients corresponding to linearly dependent columns will be 0.
      del = qr1$pivot[small_diag]
      Q1 = Q[,-small_diag]
      R1= R[-small_diag,-small_diag]
      R1_inv = solve(R1)
      tQy1 = t(Q1) %*% y
      R_inv = matrix(0,ncol(xo),ncol(xo))
      tQy = numeric(ncol(xo))
      R_inv[-del,-del]=R1_inv
      tQy[-del] = tQy1
      warning("Multicollinearity in data. Only coeffients corresponding to linearly independent columns are estimated.")
    }else{
      R_inv = solve(R)
      tQy = t(Q) %*% y
    }


    betabj = vector("list",thrs)
    rssj =  numeric(thrs)
    bic = numeric(thrs)


    RSthrs = y-Q %*% tQy
    rssj[thrs] = t(RSthrs) %*% RSthrs
    for(j in thrs:2){
        rssj[j-1] = rssj[j]+ tQy[(j+1)]^2
        bic[j] = n*log(rssj[j]/n)+(j+1)*penalty
    }

    bic[1] =  n*log(rssj[1]/n)+2*penalty

    sel = which.min(bic)
    betab = as.numeric(R_inv[1:(sel+1),1:(sel+1)] %*% tQy[1:(sel+1)])
    model_sel = order1[1:sel]

    if(bic0<bic[sel]){
        betab = beta0
        model_sel = 0
    }

    Result = list(model=model_sel,informationCriterion=bic,coefficients=betab)
    return(Result)
}





select_finalmodel_qr = function(y,x,yval,xval,order1){
# The function returns final model corresponding to minimal prediction error on validation set.

    n = length(y)
    p = ncol(x)
    x = x[,order1]
    xval = xval[,order1]


    p1 = ifelse(p<=(n-1),p,n-1)

    x = x[,1:p1]
    xval = xval[,1:p1]

    x=cbind(1,x)
    xval=cbind(1,xval)

#     qr1 = qr(x)
#     Q = qr.Q(qr1)
#     R = qr.R(qr1)
#     R_inv = solve(R)
#     tQy = t(Q) %*% y

      qr1 = qr(x)
      Q = qr.Q(qr1)
      R = qr.R(qr1)
      small_diag = which(abs(diag(R))<1e-7)
      if(length(small_diag)>0){
        #Deal with collinearities, coefficients corresponding to linearly dependent columns will be 0.
        del = qr1$pivot[small_diag]
        Q1 = Q[,-small_diag]
        R1= R[-small_diag,-small_diag]
        R1_inv = solve(R1)
        tQy1 = t(Q1) %*% y
        R_inv = matrix(0,ncol(x),ncol(x))
        tQy = numeric(ncol(x))
        R_inv[-del,-del]=R1_inv
        tQy[-del] = tQy1
        warning("Multicollinearity in data. Only coefficients corresponding to linearly independent columns are estimated.")
      }else{
        R_inv = solve(R)
        tQy = t(Q) %*% y
      }

    betab = vector("list",p1)
    betab[[p1]] = R_inv %*% tQy

    pred_err = numeric(p1)
    xval_R_inv = xval %*% R_inv
    pred = xval_R_inv %*% tQy
    pred_err[p1] = mean( (yval-pred)^2 )

    for(j in p1:2){
       betab[[j-1]] = betab[[j]] - as.numeric(tQy[j+1] * R_inv[,(j+1)])
       pred = pred - tQy[j+1] * xval_R_inv[,(j+1)]
       pred_err[j-1] = mean( (yval - pred)^2 )
    }

    jsel=which.min(pred_err)
    model=order1[1:jsel]

    Result = list(model=model,predError=pred_err,coefficients=betab[[jsel]][1:(jsel+1)])
    return(Result)
}



new.regRSM <- function()
{

    regRSM=list(scores=NULL,model=NULL,time=list(user=0,system=0,elapsed=0),
    data_transfer=list(user=0,system=0,elapsed=0),
    coefficients=NULL, predError=NULL,input_data=list(x=NULL,y=NULL),
    control=list(useGIC=NULL,selval=NULL,screening=NULL,init_weights=FALSE,parallel=NULL,m=NULL,B=NULL))

    attr(regRSM,"class")="regRSM"
    return(regRSM)
}

#regRSM = function(x,y,yval,xval,m,B, parallel,nslaves,store_data,screening,init_weights,useGIC,thrs,penalty) UseMethod("regRSM")
regRSM = function(x,...) UseMethod("regRSM")

regRSM.default = function(x,y,yval=NULL,xval=NULL,m=NULL,B=NULL, parallel="NO",nslaves=c(4),
								store_data=FALSE,screening=NULL,init_weights=FALSE,useGIC=TRUE,thrs=NULL,penalty=NULL,...)
{
    data_x = x;
    x = as.matrix(x)
    y = as.numeric(y)
    n = length(y)
    p = ncol(x)
    scores = NULL


    startTime <- proc.time()

    #Check for missing values:
    complete_cases1 = complete.cases(x)
    which_row_na = which(!complete_cases1)
    if(length(which_row_na)!=0){
      stop("Missing values in rows: ",toString(which_row_na))
    }


    # Set default values of m and B
    if(is.null(m)){
        m = floor(min(n-1,p)/2)
    }else{
        if(m>(n-2)) stop("Parameter m cannot be larger than the number of observations minus two!")
        if(m<=0) stop("Parameter m must be a positive number!")
        if(!is.wholenumber(m)) stop("Parameter m must be a whole number!")
    }
    if(is.null(B)){
        B = 1000
    }else{
        if(B<=0) stop("Parameter B must be a positive number!")
        if(!is.wholenumber(B)) stop("Parameter B must be a whole number!")
    }

    #Check for screeneing
    if(!is.null(screening))
    {
    if((screening>=1)||(screening<=0)) stop("screening must be in (0,1)")

      iw =  compute_initial_weights(y,x)
      sel = which(iw>=quantile(iw,screening))
      if(m>length(sel)) stop('Parameter m cannot be larger than the number of attributes remaining after screening procedure!')
      x = x[,sel]
    }
    #Check for initial_weights
    if(init_weights){
        initial_weights = compute_initial_weights(y,x)
    }else{
        initial_weights = NULL
    }

    #RSM method esence
    if(parallel=="MPI"){
        if("Rmpi" %in% rownames(installed.packages()) == FALSE){
                stop("Please install package 'Rmpi' to use MPI parallel version of regRSM")
        }
        requireNamespace("Rmpi",quietly = TRUE)
        create_slaves(nslaves)
        d1=proc.time()
        send_data(y,x,m,initial_weights)
        d2=proc.time()
        p = ncol(x)
        n = nrow(x)
        scores=make_experimentStatic(n,p,B,m,initial_weights)
    }else{
        if(parallel=="POSIX"){
            if("doParallel" %in% rownames(installed.packages()) == FALSE){
                stop("Please install package 'doParallel' to use POSIX parallel version of regRSM")
            }
            requireNamespace("doParallel",quietly = TRUE)
            requireNamespace("parallel",quietly = TRUE)
            d1=proc.time()
            cl = makeCluster(nslaves)
            registerDoParallel(cl)
            #registerDoParallel(nslaves)
            d2=proc.time()
            scores = compute_scores_parallel(y,x,m,B,initial_weights,nslaves)
            stopCluster(cl)
            #stopImplicitCluster()
        }else{
             d1=d2=proc.time()
             scores = compute_scores(y,x,m,B,initial_weights)
        }
    }

   #Set score 0, when variable is not selected by screeneing
    if(!is.null(screening)){
        scores1 = numeric(ncol(data_x))
        scores1[sel] = scores
        scores = scores1
    }


    selval = ifelse(!is.null(yval) && !is.null(xval),TRUE,FALSE)
    if(selval) useGIC = FALSE

    if(useGIC){
        if(is.null(penalty)){
            penalty = log(length(y))
        }else{
            if(penalty<0) stop("Penalty must be positive!")
        }
        if(is.null(thrs)){
            thrs = ifelse(p<=floor(n/2),p,floor(n/2))
        }else{
            if(thrs>min(p,(n-2))) stop("Parameter thrs cannot be larger than min(p,(n-2))!")
            if(thrs<=1) stop("Parameter thrs must be greater than one!")
            if(!is.wholenumber(thrs)) stop("Parameter thrs must be a whole number!")
        }
        order1 = sort(scores,decreasing=TRUE,index.return=TRUE)$ix
        selected_model = select_finalmodel_bic(y,data_x,order1,thrs,penalty)
        model = selected_model$model
        coefficients =  as.numeric(selected_model$coefficients)
        predError = NULL
        informationCriterion = selected_model$informationCriterion
    }else{
        if(selval==TRUE){
            order1 = sort(scores,decreasing=TRUE,index.return=TRUE)$ix
            selected_model = select_finalmodel_qr(y,data_x,yval,xval,order1)
            model = selected_model$model
            coefficients =  as.numeric(selected_model$coefficients)
            predError = selected_model$predError
            informationCriterion = NULL
        }else{
            model = NULL
            coefficients =  NULL
            predError = NULL
            informationCriterion = NULL
        }
    }
    stopTime <- proc.time()

    regRSM = new.regRSM()
    regRSM$scores = scores
    regRSM$model = model
    regRSM$time = stopTime-startTime
    regRSM$coefficients = coefficients
    regRSM$predError = predError
    regRSM$informationCriterion = informationCriterion
    regRSM$data_transfer = d2-d1
    if(store_data) { regRSM$input_data$x=data_x; regRSM$input_data$y=y }

    regRSM$control$useGIC = useGIC
    regRSM$control$selval = selval
    regRSM$control$screening = screening
    regRSM$control$init_weights =  init_weights
    regRSM$control$parallel = parallel
    regRSM$control$m = m
    regRSM$control$B = B

    print(regRSM)
    return(regRSM)
}

regRSM.formula <- function(formula, data=NULL, ...)
{
  mf = model.frame(formula,data)
  x = model.matrix(attr(mf, "terms"), data=mf)
  #Remove column corresponding to the intercept
  x = as.matrix(x[,-1])
  y = as.numeric(model.response(mf))
  est = regRSM(x, y, ...)
  cl = match.call()
  est$call = cl
  est$formula = formula
  return(est)
}



predict.regRSM = function(object,xnew,...){
    if(is.null(object$coefficients)) stop("The final model is not selected!")
    ynew <- as.numeric(cbind(1,xnew[,object$model]) %*% coef(object))
    return(ynew)
}


plot.regRSM = function(x,...){

    object = x
    if(!is.null(object$informationCriterion)){
        p1 = length(object$informationCriterion)
        plot(1:p1,object$informationCriterion,type="b",xlab="Variables",ylab="Generalized Information Criterion",main="Generalized Information Criterion",...)
    }else{
        if(!is.null(object$predError)){
            p1 = length(object$predError)
            plot(1:p1,object$predError,type="b",xlab="Variables",ylab="Prediction Error",main="Prediction Error on Validation Set",...)
        }else{
            cat("The final model is not selected! \n")
        }
    }
}


ImpPlot = function(object) UseMethod("ImpPlot")

ImpPlot.default <- function(object)
    stop("No method implemented for this class of object")

ImpPlot.regRSM = function(object){
#The function produces dot plot showing variable importances of the variables.

  if (!inherits(object, "regRSM"))
        stop("This function only works for objects of class `regRSM'")

  dotchart(object$scores,main="Variable Importance measure",xlab="RSM scores",ylab="Variable number")
}



validate = function(object,yval,xval) UseMethod("validate")

validate.default <- function(object)
    stop("No method implemented for this class of object")

validate.regRSM = function(object,yval,xval){
#The function selects the new final model based on the final scores form the original fit.

    if(is.null(object$input_data$x) || is.null(object$input_data$y)) stop('Argument store_data must be TRUE!')

    order1 = sort(object$scores,decreasing=TRUE,index.return=TRUE)$ix

    selected_model = select_finalmodel_qr(object$input_data$y,object$input_data$x,yval,xval,order1)

    model = selected_model$model
    predError = selected_model$predError
    coefficients =  as.numeric(selected_model$coefficients)


    object$model=model
    object$coefficients=coefficients
    object$predError=predError

    return(object)
}


print.regRSM = function(x,...){
#The function prints information about the model.

    object = x

    cat("\n Model summary: \n\n")

    if(object$control$useGIC){
        cat("Selection method: Generalized Information Criterion \n")
    }else{
       if(object$control$selval==TRUE){
          cat("Selection method: Validation set \n")
       }else{
          cat("Selection method: Final model is not selected! \n")
       }
    }
    if(!is.null(object$control$screening)){
        cat("Screening: yes \n")
    }else{
       cat("Screening: no \n")
    }

    if(object$control$init_weights){
        cat("Initial weights: yes \n")
    }else{
       cat("Initial weights: no \n")
    }

    if(object$control$parallel=="MPI" | object$control$parallel=="POSIX"){
        cat("Version: parallel \n")
    }else{
        cat("Version: sequential \n")
    }

    cat("Subspace size:", object$control$m, "\n")
    cat("Number of simulations:", object$control$B, "\n")

    if(object$control$parallel=="MPI" ){
        cat("Remember to use: 'mpi.close.Rslaves()' before terminate the current R session. \n")
    }

}


summary.regRSM = function(object,...){
#The function prints information about the model.

    cat("\n Model summary: \n\n")

    if(object$control$useGIC){
        cat("Selection method: Generalized Information Criterion \n")
    }else{
       if(object$control$selval==TRUE){
          cat("Selection method: Validation set \n")
       }else{
          cat("Selection method: Final model is not selected! \n")
       }
    }
    if(!is.null(object$control$screening)){
        cat("Screening: yes \n")
    }else{
       cat("Screening: no \n")
    }

    if(object$control$init_weights){
        cat("Initial weights: yes \n")
    }else{
       cat("Initial weights: no \n")
    }

    if(object$control$parallel=="MPI" | object$control$parallel=="POSIX"){
        cat("Version: parallel \n")
    }else{
        cat("Version: sequential \n")
    }

    cat("Subspace size:", object$control$m, "\n")
    cat("Number of simulations:", object$control$B, "\n")

    if(object$control$parallel=="MPI" ){
        cat("Remember to use: 'mpi.close.Rslaves()' before terminate the current R session. \n")
    }

}

roc = function(object,truemodel,plotit,...) UseMethod("roc")

roc.default <- function(object)
    stop("No method implemented for this class of object")

roc.regRSM = function(object,truemodel,plotit=TRUE,...){
# The function produces ROC curve for ordering and computes AUC.

    ordering =  sort(object$scores,decreasing=TRUE,index.return=TRUE)$ix

    if(length(setdiff(truemodel,ordering))>=1){
        stop("Argument truemodel must be a subset of ordering!")
    }

    p = length(ordering)
    npos = length(truemodel)
    nneg = p - npos

    tpr = numeric(p+1)
    fpr = numeric(p+1)

    tpr[1] = 0
    fpr[1] = 0

    for(j in 2:(p+1)){
        if(ordering[j-1] %in% truemodel){
            tpr[j] = tpr[j-1] + 1
            fpr[j] = fpr[j-1]
        }else{
            tpr[j] = tpr[j-1]
            fpr[j] = fpr[j-1] + 1
        }

    }
    tpr = tpr/npos
    fpr = fpr/nneg

    counts = as.numeric(summary(as.factor(tpr[tpr!=0])))
    xs = (counts-1)*(1/nneg)
    values = unique(tpr[tpr!=0])
    auc = sum( xs * values )

    if(plotit){
        plot(fpr,tpr,type="b",main="ROC curve for RSM ordering")
        legend("bottomright",paste("AUC= ",round(auc,4)))
    }else{
        return(auc)
    }

}
