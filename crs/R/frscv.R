## This function conducts factor regression spline cross-validation
## using exhaustive search. It takes as input a data.frame xz
## containing a mix of numeric and factor predictors and a vector y. A
## range of arguments can be provided, and one can do search on both
## the degree and knots ("degree-knots") or the degree holding the
## number of knots (segments+1) constant or the number of knots
## (segments+1) holding the degree constant. Three basis types are
## supported ("additive", "glp" or "tensor") and the argument "auto"
## will choose the basis type automatically.

frscv <- function(xz,
                  y,
                  degree.max=10,
                  segments.max=10,
                  degree.min=0,
                  segments.min=1,
                  complexity=c("degree-knots","degree","knots"),
                  knots=c("quantiles","uniform", "auto"),
                  basis=c("additive","tensor","glp","auto"),
                  cv.func=c("cv.ls","cv.gcv","cv.aic"),
                  degree=degree,
                  segments=segments,
                  tau=NULL,
                  weights=NULL,
                  singular.ok=FALSE) {

  complexity <- match.arg(complexity)
  knots <- match.arg(knots)
  basis <- match.arg(basis)
  cv.func <- match.arg(cv.func)

  t1 <- Sys.time()

  cv.objc <- function(input,
                      x,
                      z=NULL,
                      y,
                      restart,
                      num.restarts,
                      degree.max=degree.max,
                      segments.max=segments.max,
                      degree.min=degree.min,
                      segments.min=segments.min,
                      j=NULL,
                      nrow.KI.mat=NULL,
                      t2=NULL,
                      complexity=complexity,
                      knots=knots,
                      basis=basis,
                      cv.func=cv.func,
                      cv.df.min=1,
                      tau=tau,
                      weights=weights,
                      singular.ok=singular.ok) {

    if(missing(input) || missing(x) || missing(y)) stop(" you must provide input, x, y")

    ## Presumes x (continuous predictors) exist, but z
    ## (ordinal/nominal factors) can be optional

    n <- length(y)
    num.x <- NCOL(x)

    if(NROW(y) != NROW(x)) stop(" x and y have differing numbers of observations")

    ## K is a matrix, column 1 degree, column 2 segments, either or
    ## both can be determined via cv so need to take care to allow
    ## user to select knots (degree fixed), degree (knots fixed), or
    ## both degree and knots. The values used to evaluate the cv
    ## function are passed below.

    K <- round(cbind(input[1:num.x],input[(num.x+1):(2*num.x)]))

    if(!is.null(z)) {
      num.z <- NCOL(z)
      I <- round(input[(2*num.x+1):(2*num.x+num.z)])
    } else {
      num.z <- 0
      I <- NULL
    }

    cv <- cv.factor.spline.wrapper(x=x,
                                   y=y,
                                   z=z,
                                   K=K,
                                   I=I,
                                   knots=knots,
                                   basis=basis,
                                   cv.func=cv.func,
                                   cv.df.min=1,
                                   tau=tau,
                                   weights=weights,
                                   singular.ok=singular.ok)

    ## Some i/o unless options(crs.messages=FALSE)

    ## Degree is first column of K K[,1], segments second column K[,2]
    ## - could create a tmp vector for i/o, or could switch

    console <<- printClear(console)

    ## Format function...
    fw.format.2 <- function(input) sapply(input,sprintf,fmt="%#.2f")

    ## i/o depends on whether we are cross-validating degree, knots,
    ## or both

    if(complexity=="degree") {
        if(!is.null(j)) {
            if(j==1) {
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,", d[1]=",K[1,1],sep="")
            } else {
                dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
                tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                               fw.format.2(as.numeric((t2-t1),units="mins")),
                               "m",sep="")
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,tmp.0,", d[1]=",K[1,1],sep="")
            }
        } else {
            tmp.1 <- paste("d[1]=", K[1,1],sep="")
        }
        if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
    } else if(complexity=="knots") {
        if(!is.null(j)) {
            if(j==1) {
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,", s[1]=",K[1,2],sep="")
            } else {
                dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
                tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                               fw.format.2(as.numeric((t2-t1),units="mins")),
                               "m",sep="")
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,tmp.0,", s[1]=",K[1,2],sep="")
            }
        } else {
            tmp.1 <- paste("s[1]=", K[1,2],sep="")
        }
        if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")
    } else if(complexity=="degree-knots") {
        if(!is.null(j)) {
            if(j==1) {
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,", d[1]=",K[1,1],sep="")
            } else {
                dt <- (t2-t1)*(nrow.KI.mat-j+1)/j
                tmp.0 <- paste(", ",fw.format.2(as.numeric(dt,units="mins")),"/",
                               fw.format.2(as.numeric((t2-t1),units="mins")),
                               "m",sep="")
                tmp.1 <- paste("\r",j,"/",nrow.KI.mat,tmp.0,", d[1]=",K[1,1],sep="")
            }
        } else {
            tmp.1 <- paste("k[1]=", K[1,1],sep="")
        }
        if(num.x > 1) for(i in 2:num.x) tmp.1 <- paste(tmp.1, ", d[", i, "]=", K[i,1],sep="")
        for(i in 1:num.x) tmp.1 <- paste(tmp.1, ", s[", i, "]=", K[i,2],sep="")
    }

    ## For i/o for z variables...

    if(num.z > 0) for(i in 1:num.z) tmp.1 <- paste(tmp.1, ", I[", i, "]=", I[i],sep="")
    tmp.3 <- paste(", cv=", format(cv,digits=6), sep="")
    if(num.restarts > 0) {
        tmp.2 <- paste(", rs=", restart, "/", num.restarts,sep="")
        msg <- paste(tmp.1,tmp.2,tmp.3,sep="")
    } else {
        msg <- paste(tmp.1,tmp.3,sep="")
    }

    console <<- printPush(msg,console = console)

    return(cv)

  }

  console <- newLineConsole()
  console <- printPush("Working...",console = console)

  ## Take data frame x and parse into factors (z) and numeric (x)

  if(!is.data.frame(xz)) stop(" xz must be a data frame")

  xztmp <- splitFrame(xz)
  x <- xztmp$x
  z <- xztmp$z
  is.ordered.z <- xztmp$is.ordered.z
  if(is.null(z)) {
    include <- NULL
    num.z <- 0
  } else {
    num.z <- NCOL(z)
  }

  num.x <- ncol(x)
  n <- nrow(x)

  if(missing(x) || missing(y)) stop (" you must provide x and y")

  if(complexity=="degree") {
      if(missing(segments)||is.null(segments)) stop("segments missing for cross-validation of spline degree")
      if(length(segments)!=num.x) stop(" segments vector must be the same length as x")
  } else if(complexity=="knots") {
      if(missing(degree)||is.null(degree)) stop("degree missing for cross-validation of number of spline knots")
      if(length(degree)!=num.x) stop(" degree vector must be the same length as x")
  }

  ## For factor regression spline, if there is only one predictor
  ## (i.e. num.x + num.z = 1) disable auto, set to additive (tensor in
  ## this case, so don't waste time doing both).

  if((num.x+num.z==1) && (basis == "auto")) basis <- "additive"

  if(degree.min < 0 ) degree.min <- 0
  if(segments.min < 1 ) segments.min <- 1
  if(degree.max < degree.min) degree.max <- (degree.min + 1)
  if(segments.max < segments.min) segments.max <- (segments.min + 1)

  if(degree.max < 1 || segments.max < 1 ) stop(" degree.max or segments.max must be greater than or equal to 1")

  if(complexity == "degree-knots") {
    if(!is.null(z)) {
      KI.mat <- matrix.combn(K.vec1=degree.min:degree.max,K.vec2=segments.min:segments.max, num.x=num.x,num.z=num.z)
    } else {
      KI.mat <- matrix.combn(K.vec1=degree.min:degree.max, K.vec2=segments.min:segments.max,num.x=num.x)
    }
  ## We don't need to do the following again since we have specific values of segments.
  ## KI.mat[,(num.x+1):(2*num.x)] <- KI.mat[,(num.x+1):(2*num.x)]+1
  } else if(complexity == "degree") {
    if(!is.null(z)) {
      KI.mat <- matrix.combn(K.vec1=degree.min:degree.max,num.x=num.x,num.z=num.z)
      KI.mat <- cbind(KI.mat[,1:num.x],matrix(segments,nrow(KI.mat),length(segments),byrow=TRUE),KI.mat[,(num.x+1):(num.x+num.z)])
    } else {
      KI.mat <- matrix.combn(K.vec1=degree.min:degree.max,num.x=num.x)
      KI.mat <- cbind(KI.mat[,1:num.x],matrix(segments,nrow(KI.mat),length(segments),byrow=TRUE))
    }
  } else if (complexity == "knots"){
    if(!is.null(z)) {
      KI.mat <- matrix.combn(K.vec1=segments.min:segments.max,num.x=num.x,num.z=num.z)
      KI.mat <- cbind(matrix(degree,nrow(KI.mat),length(degree),byrow=TRUE),KI.mat[,1:num.x],KI.mat[,(num.x+1):(num.x+num.z)])
    } else {
      KI.mat <- matrix.combn(K.vec1=segments.min:segments.max,num.x=num.x)
      KI.mat <- cbind(matrix(degree,nrow(KI.mat),length(degree),byrow=TRUE),KI.mat[,1:num.x])
    }
  }

  ## For exhaustive search, we avoid redundant computation of cv
  ## function. Some rows will have all continuous predictors having
  ## degree 0 but different segments - computation is redundant here.
  ## Finally, if there are categorical predictors and we are using
  ## factor splines, some of these rows having no continuous
  ## predictors can differ in terms of their categorical predictors
  ## (inclusion), so determine which are unique in this case

  degree.zero.rows <- which(rowSums(KI.mat[,1:num.x,drop=FALSE])==0)

  if(length(degree.zero.rows) > 0) {

    degree.mat.zero.rows.unique <- numeric()

    if(num.z > 0) {

      KI.I.unique <- unique(KI.mat[degree.zero.rows,(2*num.x+1):(2*num.x+num.z),drop=FALSE])

      for(i in 1:nrow(KI.I.unique)) {
        degree.mat.zero.rows.unique[i] <- which(KI.mat[degree.zero.rows,(2*num.x+1):(2*num.x+num.z),drop=FALSE]==KI.I.unique[i,])[1]
      }

    } else {

      KI.I.unique <- unique(KI.mat[degree.zero.rows,1:num.x,drop=FALSE])

      for(i in 1:nrow(KI.I.unique)) {
        degree.mat.zero.rows.unique[i] <- which(KI.mat[degree.zero.rows,1:num.x,drop=FALSE]==KI.I.unique[i,])[1]
      }

    }

    degree.zero.rows <- degree.zero.rows[-degree.mat.zero.rows.unique]

    if(length(degree.zero.rows) > 0) KI.mat <- KI.mat[-degree.zero.rows,,drop=FALSE]

  }

  nrow.KI.mat <- NROW(KI.mat)
  basis.vec <- character(nrow.KI.mat)
  knots.vec <- character(nrow.KI.mat)

  ## Initialize

  cv.vec <- rep(.Machine$double.xmax,nrow.KI.mat)

  for(j in 1:nrow.KI.mat) {

    if(basis=="auto") {

      output <- cv.objc(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        degree.max=degree.max,
                        segments.max=segments.max,
                        degree.min=degree.min,
                        segments.min=segments.min,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="additive",
                        cv.func=cv.func,
                        cv.df.min=1,
                        tau=tau,
                        weights=weights,
                        singular.ok=singular.ok)

      if(output < cv.vec[j]) {
        basis.vec[j] <- "additive"
        knots.vec[j] <- attributes(output)$knots.opt
        attributes(output) <- NULL
        cv.vec[j] <- output
      }

      output <- cv.objc(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        degree.max=degree.max,
                        segments.max=segments.max,
                        degree.min=degree.min,
                        segments.min=segments.min,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="tensor",
                        cv.func=cv.func,
                        cv.df.min=1,
                        tau=tau,
                        weights=weights,
                        singular.ok=singular.ok)

      if(output < cv.vec[j]) {
        basis.vec[j] <- "tensor"
        knots.vec[j] <- attributes(output)$knots.opt
        attributes(output) <- NULL
        cv.vec[j] <- output
      }

      output <- cv.objc(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        degree.max=degree.max,
                        segments.max=segments.max,
                        degree.min=degree.min,
                        segments.min=segments.min,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis="glp",
                        cv.func=cv.func,
                        cv.df.min=1,
                        tau=tau,
                        weights=weights,
                        singular.ok=singular.ok)

      if(output < cv.vec[j]) {
        basis.vec[j] <- "glp"
        knots.vec[j] <- attributes(output)$knots.opt
        attributes(output) <- NULL
        cv.vec[j] <- output
      }

    } else {

      ## not auto, so use either "additive" or "tensor"

      output <- cv.objc(input=KI.mat[j,],
                        x=x,
                        y=y,
                        z=z,
                        degree.max=degree.max,
                        segments.max=segments.max,
                        degree.min=degree.min,
                        segments.min=segments.min,
                        restart=0,
                        num.restarts=0,
                        j=j,
                        nrow.KI.mat=nrow.KI.mat,
                        t2=Sys.time(),
                        complexity=complexity,
                        knots=knots,
                        basis=basis,
                        cv.func=cv.func,
                        cv.df.min=1,
                        tau=tau,
                        weights=weights,
                        singular.ok=singular.ok)

      if(output < cv.vec[j]) {
        basis.vec[j] <- basis
        knots.vec[j] <- attributes(output)$knots.opt
        attributes(output) <- NULL
        cv.vec[j] <- output
      }

    }

  }

  ## Sort on cv.vec

  ocv.vec <- order(cv.vec)

  cv.min <- cv.vec[ocv.vec][1]
  K.opt <- KI.mat[ocv.vec,,drop=FALSE][1,]
  basis.opt <- basis.vec[ocv.vec][1]
  knots.opt <- knots.vec[ocv.vec][1]
  degree <- K.opt[1:num.x]
  segments <- K.opt[(num.x+1):(2*num.x)]

  if(!is.null(z)) I.opt <- K.opt[(2*num.x+1):(2*num.x+num.z)]

  console <- printClear(console)
  console <- printPop(console)

  ## Set number of segments when degree==0 to 1 (or NA)

  segments[degree==0] <- 1

  if(any(degree==degree.max)) warning(paste(" optimal degree equals search maximum (", degree.max,"): rerun with larger degree.max",sep=""))
  if(any(segments==segments.max)) warning(paste(" optimal segment equals search maximum (", segments.max,"): rerun with larger segments.max",sep=""))

  if(is.null(z)) I.opt <- NULL

  crscv(K=K.opt,
        I=I.opt,
        basis=basis.opt,
        basis.vec=basis.vec,
        degree.max=degree.max,
        segments.max=segments.max,
        degree.min=degree.min,
        segments.min=segments.min,
        complexity=complexity,
        knots=knots.opt,
        degree=degree,
        segments=segments,
        K.mat=KI.mat,
        restarts=NULL,
        lambda=NULL,
        lambda.mat=NULL,
        cv.objc=cv.min,
        cv.objc.vec=as.matrix(cv.vec),
        num.x=num.x,
        cv.func=cv.func,
        tau=tau)

}
