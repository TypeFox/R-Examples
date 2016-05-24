.onAttach <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("")
    packageStartupMessage("Package ", pkgname, " (",ver,") loaded.")
}

###########################################################
## show methods for the OutlierDM class
setGeneric("show")
setMethod("show", "OutlierDM", function(object){
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  cat("Outlier Detection for Multi-replicative High-throughput Data\n\n")
  cat(" Method: ")
  switch(object@method, 
        Zscore = cat("Z-score criterion"),
        iqr = cat("Interquartile range (IQR) criterion"),
        siqr = cat("Semiinterquartile range (SIQR) criterion"),
        dixon = cat("Dixon's test"),
        grubbs = cat("Grubbs test"),
        pair = cat("Pairwise OutlierD algoirthm"),
        diff = cat("Difference-based OutlierD algorithm"),
        proj = cat("Projection-based OutlierD algorithm")
    )
  cat(" (",object@method,")","\n")

  if( object@method %in% c("pair","diff","proj")){
    cat(" Reg: ")
    switch(object@quantreg, 
           constant = cat("constant quantile regression"),
           linear = cat("linear quantile regression"),
           nonlin = cat("non-linear quantile regression"),
           nonpar = cat("non-parametric quantile regression")
    )
    cat(" (",object@quantreg,")","\n")
  }

  if( object@method %in% c("siqr", "iqr", "pair","diff","proj")) {
    if (object@method == "siqr") cat(" k: ", object@k, "for 2k * SIQR \n")
    else cat(" k: ",object@k,"for k * IQR\n")

    if (object@method %in% c("pair", "diff", "proj")) {
        cat(" Upper Quantile: ", object@contrl.para$Upper,"\n")
        cat(" Lower Quantile: ",object@contrl.para$Lower, "\n")
    }
  } else if(object@method == "Zscore"){
      cat(" k: ", object@k, "for |Z| > k \n")
  }

  cat(" Number of Observations: ", nrow(object@res),"\n")
  cat(" Number of Outliers: ", object@n.outliers,"\n")
            
  if (object@method == "pair") cat(" Centering: ",object@contrl.para$centering,"\n")
  cat(" Transformation: ", object@contrl.para$trans ,"\n")
            
  if(object@method %in% c("Zscore", "grubbs")) {
    out = object@res[, -c(2:(ncol(object@raw.data)+1))]
    if (object@method == "grubbs") out[,-1] <- round(out[,-1],3)
  } else {
    out = object@res
  }

  out[is.na(out)] <- "."
  cat("\n Head of the Output Results \n")
  print(out[1:6,, drop = FALSE], digits = 3)
  cat("To see the full information for the result, use a command, 'output(your_object_name)'. \n")            
}
)

####################################################################
setGeneric("summary")
setMethod("summary", signature(object = "OutlierDM"), 
function(object) {
  options(warn = -1)
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  cat("Outlier Detection for Multi-replicative High-throughput Data\n\n")
  cat(" Method: ")

  switch(object@method, 
        Zscore = cat("Z-score criterion"),
        iqr = cat("Interquartile range (IQR) criterion"),
        siqr = cat("Semiinterquartile range (SIQR) criterion"),
        dixon = cat("Dixon's Q-test"),
        grubbs = cat("Grubbs test"),
        pair = cat("Pairwise OutlierD algoirthm"),
        diff = cat("Difference-based OutlierD algorithm"),
        proj = cat("Projection-based OutlierD algorithm")
    )
  cat(" (",object@method,")","\n")

  if( object@method %in% c("pair","diff","proj")){
    cat(" Regression: ")
    switch(object@quantreg, 
           constant = cat("constant quantile regression"),
           linear = cat("linear quantile regression"),
           nonlin = cat("non-linear quantile regression"),
           nonpar = cat("non-parametric quantile regression")
    )
    cat(" (",object@quantreg,")","\n")
  }

  if( object@method %in% c("siqr", "iqr", "pair","diff","proj")) {
    if (object@method == "siqr") cat(" k: ", object@k, "for 2k * SIQR \n")
    else cat(" k: ",object@k,"for k * IQR\n")

    if (object@method %in% c("pair", "diff", "proj")) {
        cat(" Upper Quantile: ", object@contrl.para$Upper,"\n")
        cat(" Lower Quantile: ",object@contrl.para$Lower, "\n")
    }
  } else if(object@method == "Zscore"){
      cat(" k: ", object@k, "for |Z| > k \n")
  }

  cat(" Number of Observations: ", nrow(object@res),"\n")
  cat(" Number of Outliers: ", object@n.outliers,"\n")
            
  if (object@method == "pair") cat(" Centering: ",object@contrl.para$centering,"\n")
  cat(" Transformation: ", object@contrl.para$trans ,"\n")
            
  cat("\n Head of the Input Data \n" )
  print(object@raw.data[1:6,], digits = 3)
  cat("To see the full information of the input dataset, use a command, 'input(your_object_name)'. \n")

  cat("\n Head of the Output Results \n")

  if(object@method%in% c("Zscore", "grubbs")) {
    out = object@res[, -c(2:(ncol(object@raw.data)+1))]
  } else {
    out = object@res
  }

  out[is.na(out)] <- "."
  print(out[1:6,], digits = 3)
  cat("To see the full information for the result, use a command, 'output(your_object_name)'. \n")            

  cat('\n Head of the Peptide Numbers detected to be an Outlier\n')
  out = rownames(object@res[object@res$Outlier,])
  print(head(out, 10))
  cat("To see the full information for the candidate outliers, use a command, 'outliers(your_object_name)'. \n")
}
)

####################################################################
setGeneric("input", function(object) {
  standardGeneric("input")
})
setMethod("input", signature(object = "OutlierDM"),
  function(object) {
    res = object@raw.data
    return(res)
  }
)

####################################################################
setGeneric("output", function(object) {
  standardGeneric("output")
})
setMethod("output", signature(object = "OutlierDM"),
    function(object) object@res
)

####################################################################
setGeneric("outliers", function(object) {
  standardGeneric("outliers")
})
setMethod("outliers", signature(object = "OutlierDM"),
    function(object) object@res[object@res$Outlier,]
)

####################################################################
## plot methods for the OutlierDM class
setGeneric("plot")
setMethod("plot", signature(x = "OutlierDM", y = "missing"), 
function(x, y = NA, pch1 = 20, pch2 = 8, cex = 0.5, legend.use = TRUE, main, ...) {
  options(warn = -1)
  res <- x@res

  # NOTICE: we need parallel computing to reduce computing time

  if(x@method == "Zscore"){
    sd.index = (2+ncol(x@raw.data)):(ncol(res)-3)
    sd.index.len = length(sd.index)
    if(missing(main)) main = "Outlier Detection by the Z-score criterion"

    matplot(res[,sd.index], type = "p", cex = cex, pch = pch1, main = main,
    ylim = c(min(-x@k -1, min(res[,sd.index])), max(x@k + 1, max(res[,sd.index]))), 
    xlab = "Sample Index", axes = FALSE, col = "black", ylab = expression(Z), ...)

    for(ind in 1:nrow(res)){
      if(res$Outlier[ind]) points(rep(ind, sd.index.len), res[ind,sd.index], col = "red", pch = pch2, cex = cex +.1)
    }
    abline(h = 0, col = "blue", lty = 1.5)          
    axis(2)
    axis(1)
    abline(h = x@k, col = "tomato", lty = 1, lwd = 1.5)
    abline(h = -x@k, col = "tomato", lty = 1, lwd = 1.5)
    # Add rectangle
    rect(0, -x@k, nrow(res), x@k, col = "#0000FF08", border = '#0000FF08')

  } else if(x@method == "grubbs"){

    if(missing(main)) main = "Outlier Detection by the Grubbs criterion"

    plot(res$pvalue, pch = ifelse(res$Outlier, pch2, pch1), cex = ifelse(res$Outlier, cex + 0.1, cex), main = main,
    ylim = c(0,1), axes = FALSE, xlab = "Sample Index", ylab = "p-value",
    col = ifelse(res$Outlier, "red","black"), ...)
    # Add rectangle
    rect(0, x@contrl.para$cri.pval, nrow(res), 1, col = "#0000FF08", border = '#0000FF08')

    abline(h = x@contrl.para$cri.pval, col = "tomato", lty = 1, lwd = 1.5)
    axis(2)
    axis(1)
  } else if(x@method %in% c("iqr","siqr")){

    iqr.index = 2:(ncol(res)-5)
    iqr.index.len = length(iqr.index)

    if(missing(main)){
      if(x@method == "iqr") main = "Outlier Detection by the IQR-based criterion"
      else main = "Outlier Detection by the SIQR-based criterion"
    } 
    matplot(res[,iqr.index], type = "p", cex = cex, pch = pch1, main = main,
    xlab = "Index", axes = FALSE, col = "black", ylab = expression(Y), ...)
    for(ind in 1:nrow(res)){
      if(res$Outlier[ind]) points(rep(ind, iqr.index.len), res[ind,iqr.index], col = "red", pch = pch2, cex = cex + 0.2)
    }

    axis(2)
    axis(1)
    points(res$LB, type = "l", col = "tomato", lty = 2, lwd = 1)
    points(res$UB, type = "l", col = "tomato", lty = 2, lwd = 1)

  } else if (x@method == "proj") {

    xlab = "A"
    ylab = "M"
    plot(res$A, res$M, 
    pch = ifelse(res$Outlier, pch2, pch1),
    cex = cex, xlab = xlab, ylab = ylab, 
    col = ifelse(res$Outlier, "red", "black"), ...)
    i <- sort.list(res$A)
    lines(res$A[i], res$UB[i], col = "red", lty = 1.5)
    lines(res$A[i], res$Q3[i], col = "blue", lty = 2)
    if (legend.use) legend("topright", c("Upper Bound", "Q3"), lty = 1:2, col= c("red", "blue"))

  } else if(x@method == "diff") {
    ## Need modify t
    ## now just work when n = 3
    xlab = "A"
    ylab = "M"
    plot(res$A, res$M1, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
    col = ifelse(res$Outlier, "red", "black"), ...)     
    points(res$A, res$M2, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
    col = ifelse(res$Outlier, "red", "black"))      
    points(res$A, res$M3, pch = pch, cex = cex, xlab = xlab, ylab = ylab, 
    col = ifelse(res$Outlier, "red", "black"))      
    i <- sort.list(res$A)
    lines(res$A[i], res$UB[i], col = "red", lty = 1)
    lines(res$A[i], res$Q3[i],, col = "blue", lty = 2)
    lines(res$A[i], res$LB[i],, col = "red", lty = 1)
    lines(res$A[i], res$Q1[i],, col = "blue", lty = 2)
    if (legend.use) legend("topright", c("Upper & Lower Bound", "Q3 and Q1"), lty = 1:2, col= c("red", "blue"))

  } else if(x@method == "pair"){

    xlab = "A"
    ylab = "M"

    plot(res$A, res$M, 
      pch = ifelse(res$Outlier, pch2, pch1),
      cex = cex, xlab = xlab, ylab = ylab, 
      col = ifelse(res$Outlier, "red", "black"), ...)
    kk <- sort.list(res$A)
    lines(res$A[kk], res$UB[kk], col = "red", lty = 1)
    lines(res$A[kk], res$Q3[kk], col = "blue", lty = 2)
    lines(res$A[kk], res$LB[kk], col = "red", lty = 1)
    lines(res$A[kk], res$Q1[kk], col = "blue", lty = 2)
    if (legend.use) legend("topright", c("Upper & Lower Bound", "Q3 and Q1"), lty = 1:2, col= c("red", "blue"))
  }
  invisible(NULL)
}
)

####################################################################
setGeneric("oneplot", function(object,i, ...) {
  standardGeneric("oneplot")
})
setMethod("oneplot", signature(object = "OutlierDM", i = "numeric"),
  function(object, i = 1, pick = 0){
  # outlier detection plot for a univariate case
  ####
  # input
  # i: sample index
  ####
  # output
  # dotplot
  if(ncol(object@raw.data) <= 6) {
    warning("WARNING: \n Criteria can be described when the sample size are larger than 6.")
      i.pep = unlist(object@res[i, 2:(1+ncol(object@raw.data))])
      stripchart(i.pep, 
        method="stack",
        col = "black",
        axes = FALSE,
        pch = 20, 
        offset=0.5, 
        xlab = expression(Y))     
      axis(1)
    if(pick){
      cat("Pick a sample!! \n")
      text(p <- locator(pick), names(i.pep)[which.min(p$x - i.pep)], adj=0)
    }
    return(i.pep)
  } else{
    if(object@method == "Zscore"){
      i.pep = unlist(object@res[i, (2+ncol(object@raw.data)):(ncol(object@res)-3)])
      i.pep.non = i.pep[abs(i.pep) < object@k]
      i.pep.out = i.pep[abs(i.pep) >= object@k]
      xlab1 = expression(Z)
    
    #} else if(object@type == "grubbs"){
      #i.pep = unlist(object@res[i, 2:(1+ncol(object@raw.data))])
      #i.pep.non = i.pep[!object@outlier[i,]]
      #i.pep.out = i.pep[object@outlier[i,]]
        #xlab1 = expression(Y) 
    } else { #else if(object@type %in% c("grubbs", "iqr", "siqr")){
      i.pep = unlist(object@res[i, 2:(1+ncol(object@raw.data))])
      i.pep.non = i.pep[!object@outlier[i,]]
      i.pep.out = i.pep[object@outlier[i,]]
        xlab1 = expression(Y) 
    } 
    stripchart(i.pep.non, 
        method="jitter",
        xlim = range(i.pep),
        col = "black",
        axes = FALSE,
        pch = 20, 
        offset=0.5, 
        xlab= xlab1)     

    stripchart(i.pep.out, 
                 method="jitter",
                 xlim = c(-object@k, object@k),
                 col = "red",
                 axes = FALSE,
                 pch = 8,
                 cex = 1.3,
                 offset=0.5,
                 add = TRUE)     
      #abline(v = -object@k, col = "red")
      #abline(v = -object@k, col = "red")
      
    axis(1)
    if(pick){
      cat("Pick a sample!! \n")
      text(p <- locator(pick), names(i.pep)[which.min(p$x - i.pep)], adj=0)
    }

    return(i.pep)
  }
} 
)

####################################################################
dist.two <- function(type, q1, q2, p){
    # distance between q=(q2-q1) vector and p
    # algorithm : (p'q / q'q)q
    # type
    ## Pred : new point, A_{j}
    ## Resp : A_{j}
    # additional assumption
    # A_j =
    ## |y*v|/norm,   if y*v/v'v >= 0
    ## -|y*v|/norm,  o.w.
    norm.two <- sum((q2 - q1)^2)
    #norm.two <- sum((q2 - c(0,0,0))^2)
    if(norm.two != 0){
        u <- sum((p-q1)*(q2-q1)) / norm.two
        new.pt <- q1 + u*(q2-q1)
#       if( u < 0) new.pt <- - new.pt
        A.tmp <- sqrt(sum((new.pt - q1)^2))
        if(type == "Resp") return(new.pt)
        else if(type == "Pred") return(c(M=sqrt(sum((new.pt - p)^2)),A=ifelse(u >= 0,A.tmp, -A.tmp),new.pt))
    }
    else{
        # R^{2} space
        return(sqrt(sum( q2 - p)^2))
    }
}

####################################################################
quant.const <- function(x, Lower, Upper, method){
    if(method == "pair"){
        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2

        quant.val <- (Lower + Upper) / 2
        Q2 <- quantile(M, probs= quant.val)
        M <- M-Q2
    }
    else if(method == "proj"){
        M <- x[,1]
        A <- x[,2]
    }

    Q3 = quantile(M, probs= Upper)
    Q3 <- rep(Q3, length(A))

    Q1 = quantile(M, probs= Lower)
    Q1 <- rep(Q1, length(A))

    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.linear <- function(x, Lower, Upper, method){
  if(method == "pair"){

    M = (x[,1]-x[,2])
    A = (x[,1]+x[,2])/2

    quant.val = (Lower + Upper) / 2
    Q2 = rq(M~A, tau = quant.val )
    Q2 = Q2$coef[1]+A*Q2$coef[2]
    M <- M-Q2

    Q3 = rq(abs(M)~A, tau = quant.val )
    Q3 <- Q3$coef[1]+A*Q3$coef[2]

    Q1 <- -1 * Q3       
  } else if(method == "proj"){

    M = x[,1]
    A = x[,2]

    Q3 <- rq(abs(M)~A, tau = Upper)
    Q3 <- Q3$coef[1]+A*Q3$coef[2]

    Q1 <- rq(abs(M)~A, tau = Lower)
    Q1 <- Q1$coef[1]+A*Q1$coef[2]
  }
  return(list(M = M, A = A, Q1 = Q1, Q3 = Q3))
}

####################################################################

quant.nonlin <- function(x, Lower, Upper, method, nonlin.method, nonlin.SS, nonlin.Frank){

    FrankModel <- function(x,delta,mu,sigma,df,tau){
        z <- qt(-log(1-(1-exp(-delta))/(1+exp(-delta*pt(x,df))*((1/tau)-1)))/delta,df)
        mu + sigma*z
    }

    if(method == "pair"){
        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2
        quant.val <- (Lower + Upper) / 2        

        tmp <- try(nlrq(M ~ SSlogis(A, Asym, mid, scal), tau = quant.val ), silent = TRUE)

        if(class(tmp) != "try-error") {
           Q2 <- nlrq(M ~ SSlogis(A, Asym, mid, scal), tau= quant.val)
           Q2 <- predict(Q2, newdata=list(x=A))
           M <- M-Q2
           Q3 <-  switch( nonlin.SS,
                Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = quant.val),
                    tau = quant.val, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
                    method = nonlin.method ), silent = TRUE),
                Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = quant.val, method = nonlin.method ), silent = TRUE),
                Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = quant.val, method = nonlin.method  ), silent = TRUE),
                Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = quant.val, method = nonlin.method  ), silent = TRUE),
                biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = quant.val, method = nonlin.method  ), silent = TRUE),
                AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = quant.val, method = nonlin.method  ), silent = TRUE)
            )
            if(class(Q3)=="try-error") {
               fit <- quant.linear(x, Lower, Upper, method)
               return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
            }
            Q3 <- predict(Q3, newdata=list(x=A))
            Q1 <- -Q3
        }
        else {
           fit <- quant.linear(x, Lower, Upper, method)
           return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
        }
    }
    else if(method == "proj" ){
        M <- x[,1]
        A <- x[,2]
        
        Q3 <-  switch( nonlin.SS,
                Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = Upper),
                    tau = Upper, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
                    method = nonlin.method ), silent = TRUE),
                Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = Upper, method = nonlin.method ), silent = TRUE),
                Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = Upper, method = nonlin.method  ), silent = TRUE),
                Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = Upper, method = nonlin.method  ), silent = TRUE),
                biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = Upper, method = nonlin.method  ), silent = TRUE),
                AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = Upper, method = nonlin.method  ), silent = TRUE)
                )
        if(class(Q3)=="try-error") {
           fit <- quant.linear(x, Lower, Upper, method)
           return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
        }

        Q3 <- predict(Q3, newdata=list(x=A))

        Q1 <-  switch( nonlin.SS,
                Frank = try( nlrq( abs(M) ~ FrankModel(A, delta, mu, sigma, df = nonlin.Frank[1], tau = Lower),
                    tau = Lower, start=list(delta=nonlin.Frank[2], mu = nonlin.Frank[3], sigma = nonlin.Frank[4]), 
                    method = nonlin.method ), silent = TRUE),
                Self = try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau = Lower, method = nonlin.method ), silent = TRUE),
                Asym = try(nlrq(abs(M) ~ SSasymp( A, Asym, R0, lrc), tau = Lower, method = nonlin.method  ), silent = TRUE),
                Orig = try(nlrq(abs(M) ~ SSasympOrig(A, Asym, lrc), tau = Lower, method = nonlin.method  ), silent = TRUE),
                biexp = try(nlrq(abs(M) ~ SSbiexp(A, A1, lrc1, A2, lrc2), tau = Lower, method = nonlin.method  ), silent = TRUE),
                AsymOff = try(nlrq(abs(M) ~ SSasympOff(A, Asym, lrc, c0), tau = Lower, method = nonlin.method  ), silent = TRUE)
                )
        if(class(Q1)=="try-error") {
           fit <- quant.linear(x, Lower, Upper, method)
           return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
        }

        Q1 <- predict(Q1, newdata=list(x=A))
    }   
    
    if(length(which(abs(M) < abs(Q1))) < length(M)/4 & method == "pair") {
       fit <- quant.linear(x, Lower, Upper, method)
       return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
    }
    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}


####################################################################

quant.nonpar <- function(x, Lower, Upper, method, lbda){

    if(method == "pair"){
        M <- (x[,2]-x[,1])
        A <- (x[,2]+x[,1])/2
    }
    else if(method == "proj"){
        M <- x[,1]
        A <- x[,2]
    }
    #First Quantile
    Q1 <- rqss(M~qss(A, lambda=lbda), tau = Lower )
    i <- match(A, Q1$qss[[1]]$xyz[-1,1])
    y  <- Q1$coef[-1]
    Q1 <- Q1$coef[1] + y[i]
    #Q1 <- -Q1

    #Third Quantile
    Q3 <- rqss(M~qss(A, lambda=lbda), tau = Upper )
    i <- match(A, Q3$qss[[1]]$xyz[-1,1])
    y  <- Q3$coef[-1]
    Q3 <- Q3$coef[1] + y[i]
    return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################

quant.const.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
    if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
    M <- x - A.med
    M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]   

    A <- apply(x, 1, mean, na.rm = TRUE)
    
    Q1 <- quantile(M.long, probs=Lower)
    Q3 <- quantile(M.long, probs=Upper)

    Q1 <- rep(Q1, length(A))
    Q3 <- rep(Q3, length(A))

    return(list(M = M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.linear.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
    if(is.data.frame(x)) x <- as.data.frame(x)
    A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
    M <- x - A.med
    M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]   

    A <- apply(x, 1, mean, na.rm = TRUE)
    A <- rep(A,  n.obs)
    
    quant.val <- (Lower + Upper) / 2
    Q3 <- rq(abs(M.long) ~ A, tau= quant.val)
    Q3 <- Q3$coef[1]+A*Q3$coef[2]

    Q1 <- -1 * Q3
    return(list(M = M, A = A[1:n.para], Q1 = Q1[1:n.para], Q3 = Q3[1:n.para]))
}


####################################################################

quant.nonlin.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para){
  if(is.data.frame(x)) x <- as.data.frame(x)
  A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
  M <- x - A.med
  M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]   

  A <- apply(x, 1, mean, na.rm = TRUE)
  A <- rep(A,  n.obs)

  quant.val <- (Lower + Upper) / 2
    
  tmp <- try(nlrq(abs(M.long) ~ SSlogis(A, Asym, mid, scal), tau= quant.val ), silent = TRUE)
  if(class(tmp) =="try-error") {
     fit <- quant.linear.D(x,dist.mthd, Lower, Upper, n.obs, n.para)
     return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
  }

  Q3 <- nlrq(abs(M.long) ~ SSlogis(A, Asym, mid, scal), tau= quant.val)
  Q3 <- predict(Q3, newdata=list(x=A))

  Q1 <- -1 * Q3
    
  return(list(M = M, A=A[1:n.para], Q1=Q1[1:n.para], Q3=Q3[1:n.para]))
}
####################################################################

quant.nonpar.D <- function(x, dist.mthd, Lower, Upper, n.obs, n.para, lbda){
  if(is.data.frame(x)) x <- as.data.frame(x)
  A.med <- apply(x,1,dist.mthd,na.rm=TRUE)
  M <- x - A.med
  M.long <- reshape(M, direction = "long", varying = list(1:n.obs))[,2]   

  A <- apply(x, 1, mean, na.rm = TRUE)
  A <- rep(A,  n.obs)
    
  #First Quantile
  Q1 <- rqss(M.long~qss(A, lambda=lbda), tau= Lower)
  i <- match(A, Q1$qss[[1]]$xyz[-1,1])
  y  <- Q1$coef[-1]
  Q1 <- Q1$coef[1] + y[i]

  #Third Quantile
  Q3 <- rqss(M.long~qss(A, lambda=lbda), tau= Upper )
  i <- match(A, Q3$qss[[1]]$xyz[-1,1])
  y  <- Q3$coef[-1]
  Q3 <- Q3$coef[1] + y[i]

  return(list(M = M, A=A[1:n.para], Q1=Q1[1:n.para], Q3=Q3[1:n.para]))
}

####################################################################
rgrubbs.test <- function(x, alpha = 0.05){
    #### input 
    # x: data
    # alpha: significance level alpha for a p-value
    #### output
    # x.res: Boolean vector for outliers
    x.res = matrix( NA, nrow = 2, ncol = length(x) )
    for(i in 1:(length(x)-2)) {

        n = length(x)
        xbar = mean(x)
        xsd = sd(x)

        Gmin = (min(x) - xbar) / xsd
        Gmax = (max(x) - xbar) / xsd
        pval.min = 1 - pgrubbs(Gmin, n, type = 10)
        pval.max = 1 - pgrubbs(Gmax, n, type = 10)

    if(min(pval.min, pval.max) > alpha) break
    
        if(pval.min <= pval.max){
            x.sel = which.min(x)
            x.res[1,x.sel] <- Gmin
            x.res[2,x.sel] <- pval.min
        } else{
            x.sel = which.max(x)
            x.res[1,x.sel] <- Gmax
            x.res[2,x.sel] <- pval.max
        }
        x <- x[-x.sel]
    }
    return(x.res)
}
####################################################################