#  *****************************************************************************
#   File : far.R
#         ************************************************************
#   Description :
#       Functional Autoregressive functions and methods
#   Version : 2.2
#   Date : 2007-10-01
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

#  *****************************************************************************
#   Title : far
#         ************************************************************
#   Description :
#       Modelization of Vectorized Functional Processes
#   Version : 2.1
#   Date : 2005-01-10
#  *****************************************************************************
far <- function(data, y, x, kn, center=TRUE, na.rm=TRUE, joined=FALSE)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata")) # test the type of data
        stop("data is not of class fdata")
    call <- match.call()

    # find dimensions
    n <- ncol(data[[1]])

    # find variables
    if (missing(y))
    {
        if (missing(x))
        {
            x <- NULL
            y <- names(data)
        } else {
            y <- x
            x <- NULL
        }
    } else {
        if (missing(x)) x <- NULL
    }

    # find dimensions and test
    variables <- c(y,x)
    nx <- length(x)
    ny <- length(y)
    r <- nx+ny

    if (joined) # if joined estimation
    {
        if (missing(kn)) kn <- r
        if (1 != length(kn)) {
            stop("Gives only one kn in joined estimation.")
        }
    } else {
        if (missing(kn)) kn <- rep(1,r)
        if (r != length(kn)) {
            stop("Gives a kn value for each variable. Dimension are different.")
        }
    }

    # adapt data to the model chosen
    data.adapt <- list()
    if (nx>0) n <- n-1
    for (i in 1:length(y))
        data.adapt[[y[i]]] <- (data[[y[i]]])[,1:n,drop=FALSE]
    if (nx>0) for (i in 1:length(x))
        data.adapt[[x[i]]] <- (data[[x[i]]])[,-1,drop=FALSE]
    class(data.adapt) <- "fdata"

    # Removing of non available data if required
    if (na.rm) {
        listobs <- c(apply(!is.na(data.adapt),2,all))
        listobs2 <- c(FALSE,listobs) * c(listobs,FALSE) == 1
    } else {
        listobs <- rep(TRUE,n)
        listobs2 <- c(FALSE,rep(TRUE,n-1),FALSE)
    }

    nbobs <- sum(listobs == TRUE)
    nbobs2 <- sum(listobs2 == TRUE)

    # centering
    if (center) {
        f1 <- function(x,listobs) matrix(apply(x[,listobs],1,mean),ncol=1)
        databar <- lapply(data.adapt,f1,listobs)
        class(databar)<-"fdata"
        data <- list()
        for (i in 1:length(variables))
            data[[variables[i]]] <-
                sweep(data.adapt[[variables[i]]],1,databar[[variables[i]]],"-")
        class(data) <- "fdata"
    } else {
        databar <- NULL
        data<-data.adapt
    }

    # Begining of the estimation
    # --------------------------

    if (joined)
    {
        eigenvector <- list()
        eigenvalues <- list()
        length(eigenvector) <- 1
        length(eigenvalues) <- 1

        # Calculation of the subspaces obtained from the covariance matrix
        nrowdata <- c(0,cumsum(unlist(lapply(data,nrow))))
        datacent <- matrix(0,nrow=nrowdata[r+1],ncol=nbobs)
        for (i in 1:r)
            datacent[(nrowdata[i]+1):nrowdata[i+1],] <-
                (data[[i]])[,listobs,drop=FALSE]
        sdbase <- eigen(datacent %*% t(datacent / nbobs))
        eigenvector[[1]] <- sdbase$vectors[, 1:kn,drop=FALSE]
        eigenvalues[[1]] <- as.matrix(sdbase$values/nrowdata[r+1])

        # Determination of the projection matrix
        datacent <- matrix(0,nrow=nrowdata[r+1],ncol=n)
        for (i in 1:r)
            datacent[(nrowdata[i]+1):nrowdata[i+1],] <- data[[i]]
        Proj <- t(eigenvector[[1]]) %*% datacent
    } else {
        eigenvector <- list()
        eigenvalues <- list()
        Projdata <- list()
        length(eigenvector) <- r
        length(eigenvalues) <- r
        length(Projdata) <- r

        # Calculation of the subspaces obtained from the covariance matrix
        for (i in 1:r)
        {
            datacent <- (data[[i]])[,listobs]
            sdbase <- eigen(datacent %*% t(datacent / nbobs))
            eigenvector[[i]] <- sdbase$vectors[, 1:kn[i],drop=FALSE]
            eigenvalues[[i]] <- as.matrix(sdbase$values/nrow(datacent))
            Projdata[[i]] <- t(eigenvector[[i]]) %*% data[[i]]
        }
          # Determination of the projection matrix
        Proj <- matrix(0,ncol=n,nrow=sum(kn))
        kkn <- c(0,kn)
        for (k in 1:r)
            Proj[sum(kkn[1:k])+(1:kkn[k+1]),] <- Projdata[[k]]
    }

    # Calculation of the correlation matrix rho
        Delta <- Proj[,listobs2[-(n+1)],drop=FALSE] %*%
                    t(Proj[,listobs2[-1],drop=FALSE])
        InvG <- invgen(Proj[,listobs,drop=FALSE] %*%
                    t(Proj[,listobs,drop=FALSE]))
        rho <- Delta %*% InvG * nbobs / nbobs2

    # result
    output <- list(
        call = call,
        data = data,
        databar = databar,
        y = y,
        x = x,
        v = eigenvector,
        values = eigenvalues,
        rho = rho,
        nbvar = r,
        kn = kn,
        joined = joined)
    class(output) <- "far"
    return(output)
}

#  *****************************************************************************
#   Title : print.far
#         ************************************************************
#   Description :
#       print method for the 'far' model
#   Version : 1.0
#   Date : 2001-03-27
#  *****************************************************************************
print.far<-function(x, ..., digits = max(3, getOption("digits") - 3),
            na.print = "", file="", append=TRUE)
{
    variables <- c(x$y,x$x)
    cat("Functional Autoregressive Model\n",file=file,append=append)
    cat("Call: ", deparse(x$call), "\n\n",file=file,append=append)

    if (x$joined)
    {
        cat("Joined variable\n",file=file,append=append)
        cat("Dimension of the subspace: ", format(x$kn, digits = digits),
        "\n",file=file,append=append)
        var.explained <- (x$values[[1]])^2
        cat("Explained Variance: ", format(sum(var.explained[1:x$kn[1]])/
                                           sum(var.explained)*100,
            digits = digits), "%\n",file=file,append=append)
        cat("Estimated first Eigen values of the Covariance: ",
            format((x$values[[1]])[1:x$kn], digits = digits),
            "\n\n",file=file,append=append)
    } else {
        # printed for each variable
        for (i in 1:length(x$kn))
        {
            cat("Variable: ", variables[i], "\n",file=file,append=append)
            cat("Dimension of the subspace: ", format(x$kn[[i]],
                 digits = digits), "\n",file=file,append=append)
            var.explained <- (x$values[[i]])^2
            cat("Explained Variance: ", format(sum(var.explained[1:x$kn[i]])/
                                               sum(var.explained)*100,
                digits = digits), "%\n",file=file,append=append)
            cat("Estimated first Eigen values of the Covariance: ",
                format((x$values[[i]])[1:x$kn[i]], digits = digits),
                "\n\n",file=file,append=append)
        }
    }

    cat("Estimated correlation Matrix in adequate subspace: \n",
        file=file,append=append)
    if (file=="")
        print(round(x$rho,3))
    else
        for (i in 1:nrow(x$rho))
            cat(format(x$rho[i,],digits=digits),"\n",file=file,append=append)
    cat("\n",file=file,append=append)
    invisible(x)
}

#  *****************************************************************************
#   Title : coef.far
#         ************************************************************
#   Description :
#       coef method for the 'far' model
#   Version : 1.1
#   Date : 2003-06-11
#  *****************************************************************************
coef.far<-function (object, ...)
{
    return(object$rho)
}

#  *****************************************************************************
#   Title : plot.far
#         ************************************************************
#   Description :
#       plot method for the 'far' model
#   Version : 2.0
#   Date : 2001-07-06
#  *****************************************************************************
plot.far <- function(x,...)
{
    xval <- rownames((x$data)[[1]])
    kn <- x$kn
    n <- length(x$kn)
    names <- names(x$data)
    if (x$joined)
    {
        matplot(x=1:nrow(x$v[[1]]),y=x$v[[1]],type='l',xlab="time",
            ylab="",main=paste(names,collapse=", "),...)
        range.plot <- (par()$usr[c(1,4)])
        legend(x=range.plot[1],y=range.plot[2],
                legend=paste("v",1:kn[1]),lty=1:kn[1],col=1:kn[1])

    } else {
        for (i in 1:n)
        {
            matplot(x=xval,y=x$v[[i]],type='l',xlab="time",
                ylab="",main=paste(names[i]),...)
            range.plot <- (par()$usr[c(1,4)])
            legend(x=range.plot[1],y=range.plot[2],
                    legend=paste("v",1:kn[i]),lty=1:kn[i],col=1:kn[i])
        }
    }
    invisible()
}

#  *****************************************************************************
#   Title : predict.far
#         ************************************************************
#   Description :
#       Computation of prediction for the class model "far"
#   Version : 2.0
#   Date : 2001-07-09
#  *****************************************************************************
predict.far<-function(object, ..., newdata = NULL, label, na.rm=TRUE,
                      positive=FALSE)
{
    if ( (!is.null(class(object)))
        && (class((object)) != "far")) # test the type of data
        stop("object is not of class far")
    if ( (!is.null(class(newdata)))
        && (class((newdata)) != "fdata")) # test the type of data
        stop("newdata is not of class fdata")

    x <- object$x
    y <- object$y
    nx <- length(x)
    ny <- length(y)
    r <- nx+ny

    n <- ncol(newdata[[object$y[1]]])
    if (nx>0) # if there is auxiliary variables
    {
        label <- (colnames(newdata[[y[1]]]))[-1]
        data <- list()
        if (is.null(object$databar))
        {
            for (i in 1:ny)
                data[[y[i]]] <- (newdata[[y[i]]])[,-n,drop=FALSE]
            for (i in 1:nx)
                data[[x[i]]] <- (newdata[[x[i]]])[,-1,drop=FALSE]
        } else {
            for (i in 1:ny)
                data[[y[i]]] <- sweep((newdata[[y[i]]])[,-n,drop=FALSE],1,
                                      object$databar[[y[i]]],"-")
            for (i in 1:nx)
                data[[x[i]]] <- sweep((newdata[[x[i]]])[,-1,drop=FALSE],1,
                                      object$databar[[x[i]]],"-")
        }
        class(data) <- "fdata"
        n <- (n-1)
    } else {
        if (missing(label))
            label <- c(colnames(newdata[[y[1]]])[-1],paste(n+1))
        else
            label <- c(colnames(newdata[[y[1]]])[-1],label)
        data <- list()
        if (is.null(object$databar))
        {
            for (i in 1:ny) {
                data[[y[i]]] <- newdata[[y[i]]]
            }
        } else {
            for (i in 1:ny) {
                data[[y[i]]] <- sweep((newdata[[y[i]]]),1,
                                      object$databar[[y[i]]],"-")
            }
        }
        class(data) <- "fdata"
    }

    kn <- object$kn
    if (na.rm) {
        listobs <- c(apply(!is.na(data),2,all))
    } else {
        listobs <- rep(TRUE,n)
    }
    nbobs <- sum(listobs==TRUE)

    if (object$joined)
    {
        nrowdata <- c(0,cumsum(unlist(lapply(data,nrow))))
        datacent <- matrix(0,ncol=nbobs,nrow=nrowdata[r+1])
        for (k in (1:r)) {
            datacent[(nrowdata[k]+1):nrowdata[k+1],] <-
                (data[[k]])[,listobs,drop=FALSE]
        }
        datacent <- t(object$v[[1]]) %*% datacent
        pred <- list()
        length(pred) <- ny
        pred2 <- object$v[[1]] %*% (object$rho %*% datacent)
        for (i in (1:ny)) {
            pred[[i]] <- pred2[(nrowdata[i]+1):nrowdata[i+1],,drop=FALSE]
        }
    } else {
        datacent <- matrix(0,ncol=nbobs,nrow=sum(kn))
        kkn <- c(0,kn)
        for (k in (1:r)) {
            datacent[sum(kkn[1:k])+(1:kkn[k+1]),] <-
                    t(object$v[[k]]) %*% ((data[[k]])[,listobs,drop=FALSE])
        }
        pred <- list()
        length(pred) <- ny
        for (i in (1:ny)) {
            pred[[i]] <- object$v[[i]] %*% (object$rho %*%
                         datacent)[sum(kkn[1:i])+(1:kkn[i+1]),,drop=FALSE]
        }
    }
    for (i in (1:ny))
    {
        if (!is.null(object$databar))
            pred[[i]] <- sweep((pred[[i]]),1,object$databar[[i]],"+")
        if (positive)
            pred[[i]] <- (pred[[i]]+abs(pred[[i]]))/2
        rownames(pred[[i]]) <- rownames(data[[i]])
        colnames(pred[[i]]) <- label[listobs]
    }

    names(pred) <- object$y
    class(pred) <- "fdata"
    return(pred)
}

#  *****************************************************************************
#   Title : far.cv
#         ************************************************************
#   Description :
#       Croos validation for the Model of Vectorized Functional Processes
#   Version : 1.2
#   Date : 2007-10-01
#  *****************************************************************************
far.cv <- function(data, y, x, kn, ncv, cvcrit, center=TRUE, na.rm=TRUE,
                   joined=FALSE)
{
    if (class(data) != "fdata") # test the type of data
        stop("data is not of class fdata")
    call <- match.call()

    # find dimensions
    n <- ncol(data[[1]])
    if (missing(ncv)) ncv <- round(n/5)
    n1 <- (n-ncv)

    # find variables
    if (missing(y))
    {
        if (missing(x))
        {
            x <- NULL
            y <- names(data)
        } else {
            y <- x
            x <- NULL
        }
    } else {
        if (missing(x)) x <- NULL
    }
    if (missing(cvcrit))
    {
        cvcrit <- y
    }

    # find dimensions and test
    variables <- c(y,x)
    nx <- length(x)
    ny <- length(y)
    ncrit <- length(cvcrit)
    r <- nx+ny
    dim1 <- unlist(lapply(data,nrow))
    dim1 <- dim1[variables]

    if (joined) # if joined estimation
    {
        if (missing(kn)) kn <- sum(dim1)
        if (1 != length(kn)) stop("Gives only one kn in joined estimation.")
    } else {
        if (missing(kn)) kn <- dim1
        if (r != length(kn))
            stop("Gives a kn value for each variable. Dimension are different.")
    }

    # adapt data to the model chosen
    data.apprent <- list()
    data.test <- list()
    if (nx>0)
    {
        n1 <- n1-1
    }
    for (i in 1:length(y))
    {
        data.apprent[[y[i]]] <- (data[[y[i]]])[,1:n1,drop=FALSE]
        data.test[[y[i]]] <- (data[[y[i]]])[,n1+(1:ncv),drop=FALSE]
    }
    if (nx>0) for (i in 1:length(x))
    {
        data.apprent[[x[i]]] <- (data[[x[i]]])[,1+(1:n1),drop=FALSE]
        data.test[[x[i]]] <- (data[[x[i]]])[,n1+1+(1:ncv),drop=FALSE]
    }
    class(data.apprent) <- "fdata"
    class(data.test) <- "fdata"

    # Removing non available data if required
    if (na.rm) {
        listobs <- c(apply(!is.na(data.apprent),2,all))
        listobs2 <- c(FALSE,listobs) * c(listobs,FALSE) == 1
        listobs.test <- c(apply(!is.na(data.test),2,all))
    } else {
        listobs <- rep(TRUE,n1)
        listobs2 <- c(FALSE,rep(TRUE,n1-1),FALSE)
        listobs.test <- rep(TRUE,ncv)
    }

    nbobs <- sum(listobs == TRUE)
    nbobs2 <- sum(listobs2 == TRUE)
    nbobs.test <- sum(listobs.test == TRUE)

    # centering
    if (center) {
        f0 <- function(x,listobs) matrix(apply(x[,listobs],1,mean),ncol=1)
        databar <- lapply(data.apprent,f0,listobs)
        class(databar)<-"fdata"
        data <- list()
        for (i in 1:length(variables))
            data[[variables[i]]] <-
                sweep(data.apprent[[variables[i]]],1,
                      databar[[variables[i]]],"-")
        class(data) <- "fdata"
        data2 <- list()
        for (i in 1:length(variables))
            data2[[variables[i]]] <-
                sweep(data.test[[variables[i]]],1,
                      databar[[variables[i]]],"-")
        class(data2) <- "fdata"
    } else {
        databar <- NULL
        data<-data.apprent
        data2<-data.test
    }

    # Begining the estimation
    # -----------------------

    if (joined)
    {
        eigenvector <- list()
        eigenvalues <- list()
        length(eigenvector) <- 1
        length(eigenvalues) <- 1

        # Calculation of the subspaces obtained from the covariance matrix
        nrowdata <- c(0,cumsum(dim1))
        datacent <- matrix(0,nrow=nrowdata[r+1],ncol=nbobs)
        for (i in 1:r)
            datacent[(nrowdata[i]+1):nrowdata[i+1],] <-
                (data[[i]])[,listobs,drop=FALSE]
        sdbase <- eigen(datacent %*% t(datacent / nbobs))
        eigenvector[[1]] <- sdbase$vectors[, 1:kn,drop=FALSE]
        eigenvalues[[1]] <- as.matrix(sdbase$values/nrowdata[r+1])

        # Determination of the projection matrix
        datacent <- matrix(0,nrow=nrowdata[r+1],ncol=n1)
        datacent2 <- matrix(0,nrow=nrowdata[r+1],ncol=nbobs.test)
        for (i in 1:r)
        {
            datacent[(nrowdata[i]+1):nrowdata[i+1],] <- data[[i]]
            datacent2[(nrowdata[i]+1):nrowdata[i+1],] <-
                data2[[i]][,listobs.test,drop=FALSE]
        }
    } else {
        eigenvector <- list()
        eigenvalues <- list()
        Projdata <- list()
        Projdata2 <- list()
        length(eigenvector) <- r
        length(eigenvalues) <- r
        length(Projdata) <- r
        length(Projdata2) <- r

        # Calculation of the subspaces obtained from the covariance matrix
        for (i in 1:r)
        {
            datacent <- (data[[i]])[,listobs]
            sdbase <- eigen(datacent %*% t(datacent / nbobs))
            eigenvector[[i]] <- sdbase$vectors[, 1:kn[i],drop=FALSE]
            eigenvalues[[i]] <- as.matrix(sdbase$values/nrow(datacent))
            Projdata[[i]] <- t(eigenvector[[i]]) %*% data[[i]]
            Projdata2[[i]] <- t(eigenvector[[i]]) %*%
                              data2[[i]][,listobs.test,drop=FALSE]
        }
    }

    # Begining of the cross validation
    # --------------------------------
    output <- matrix(0,ncol=length(kn)+6,nrow=prod(kn))
    f1<-function(x) mean(apply(abs(x),2,mean),na.rm=TRUE)
    f2<-function(x) mean(sqrt(apply(x^2,2,mean)),na.rm=TRUE)
    f3<-function(x) mean(apply(abs(x),2,max),na.rm=TRUE)
    f4<-function(x) mean(abs(x),na.rm=TRUE)
    f5<-function(x) sqrt(mean(x^2,na.rm=TRUE))
    f6<-function(x) max(abs(x),na.rm=TRUE)

    pk <- prod(kn)
    lk <- length(kn)

    if (joined)
    {
        for (k in 1:kn)
        {
            # projection
            Proj <- t(eigenvector[[1]][,1:k,drop=FALSE]) %*% datacent
            Proj2 <- t(eigenvector[[1]][,1:k,drop=FALSE]) %*% datacent2

            # Calculation of the correlation matrix rho
            Delta <- Proj[,listobs2[-(n1+1)],drop=FALSE] %*%
                     t(Proj[,listobs2[-1],drop=FALSE])
            InvG <- invgen(Proj[,listobs,drop=FALSE] %*%
                           t(Proj[,listobs,drop=FALSE]))
            rho <- Delta %*% InvG * nbobs / nbobs2

            # Prediction
            pred2 <- (eigenvector[[1]][,1:k,drop=FALSE]) %*% rho %*% Proj2
            pred <- list()
            pred.max <- list()
            for (i in 1:ny)
            {
                pred[[y[i]]] <- pred2[(nrowdata[i]+1):nrowdata[i+1],
                                      -nbobs.test,drop=FALSE]
                rownames(pred[[y[i]]]) <- rownames(data[[y[i]]])
                colnames(pred[[y[i]]]) <-
                   (colnames(data2[[y[i]]])[c(FALSE,listobs.test)])[-nbobs.test]
            }

            # Calculation of errors
            for (i in 1:ncrit)
            {
                pred.max[[cvcrit[i]]] <- apply((data2[[cvcrit[i]]])[,
                        colnames(pred[[cvcrit[i]]]),drop=FALSE],2,max)-
                        apply(pred[[cvcrit[i]]],2,max)
                pred[[cvcrit[i]]] <- ((data2[[cvcrit[i]]])[,
                        colnames(pred[[cvcrit[i]]]),drop=FALSE]-
                        pred[[cvcrit[i]]])
            }

            output[k,2] <- mean(unlist(lapply(pred,f1)))
            output[k,3] <- mean(unlist(lapply(pred,f2)))
            output[k,4] <- mean(unlist(lapply(pred,f3)))
            output[k,5] <- mean(unlist(lapply(pred.max,f4)))
            output[k,6] <- mean(unlist(lapply(pred.max,f5)))
            output[k,7] <- mean(unlist(lapply(pred.max,f6)))
            output[k,1] <- k
        }


    } else {
        kn2<-c(1,kn)
        for (i in 1:length(kn))
            output[,i]<-rep(rep(1:kn[i],rep(pk/prod(kn[1:i]),kn[i])),
                            prod(kn2[1:i]))

        for (j in 1:pk)
        {
            kn <- output[j,1:lk]

            # Determination of the projection matrix
            Proj <- matrix(0,ncol=n1,nrow=sum(kn))
            Proj2 <- matrix(0,ncol=nbobs.test,nrow=sum(kn))
            kkn <- c(0,kn)
            for (k in 1:r)
            {
                Proj[sum(kkn[1:k])+(1:kkn[k+1]),] <-
                    Projdata[[k]][1:kn[k],,drop=FALSE]
                Proj2[sum(kkn[1:k])+(1:kkn[k+1]),] <-
                    Projdata2[[k]][1:kn[k],,drop=FALSE]
            }

            # Calculation of the correlation matrix rho
            Delta <- Proj[,listobs2[-(n1+1)],drop=FALSE] %*%
                     t(Proj[,listobs2[-1],drop=FALSE])
            InvG <- invgen(Proj[,listobs,drop=FALSE] %*%
                           t(Proj[,listobs,drop=FALSE]))
            rho <- Delta %*% InvG * nbobs / nbobs2

            # Prediction
            pred <- list()
            pred.max <- list()
            for (i in 1:ny)
            {
                pred[[y[i]]] <- (eigenvector[[1]][,1:kn[i],drop=FALSE]) %*%
                    (rho %*% Proj2)[sum(kkn[1:i])+(1:kkn[i+1]),
                                    -nbobs.test,drop=FALSE]
                rownames(pred[[y[i]]]) <- rownames(data[[y[i]]])
                colnames(pred[[y[i]]]) <-
                   (colnames(data2[[y[i]]])[c(FALSE,listobs.test)])[-nbobs.test]
            }

            # Calculation of errors
            for (i in 1:ncrit)
            {
                pred.max[[cvcrit[i]]] <-
                    apply((data2[[cvcrit[i]]])[,colnames(pred[[cvcrit[i]]]),
                               drop=FALSE],2,max)-apply(pred[[cvcrit[i]]],2,max)
                pred[[cvcrit[i]]] <-
                    ((data2[[cvcrit[i]]])[,colnames(pred[[cvcrit[i]]]),
                               drop=FALSE]-pred[[cvcrit[i]]])
            }
            output[j,lk+1] <- mean(unlist(lapply(pred,f1)))
            output[j,lk+2] <- mean(unlist(lapply(pred,f2)))
            output[j,lk+3] <- mean(unlist(lapply(pred,f3)))
            output[j,lk+4] <- mean(unlist(lapply(pred.max,f4)))
            output[j,lk+5] <- mean(unlist(lapply(pred.max,f5)))
            output[j,lk+6] <- mean(unlist(lapply(pred.max,f6)))
        }
    }

    # result
    dimnames(output) <- list(NULL,c(paste("k",1:lk,sep=""),
                                    "L1","L2","Linf","L1max","L2max","Linfmax"))
    moutput1<-min(output[,length(kn)+1])
    moutput1<-output[output[,length(kn)+1]==moutput1,]
    moutput2<-min(output[,length(kn)+2])
    moutput2<-output[output[,length(kn)+2]==moutput2,]
    moutput3<-min(output[,length(kn)+3])
    moutput3<-output[output[,length(kn)+3]==moutput3,]
    moutput4<-min(output[,length(kn)+4])
    moutput4<-output[output[,length(kn)+4]==moutput4,]
    moutput5<-min(output[,length(kn)+5])
    moutput5<-output[output[,length(kn)+5]==moutput5,]
    moutput6<-min(output[,length(kn)+6])
    moutput6<-output[output[,length(kn)+6]==moutput6,]
    invisible(list("cv"=output,
                   "minL1"=moutput1,
                   "minL2"=moutput2,
                   "minLinf"=moutput3,
                   "minL1max"=moutput4,
                   "minL2max"=moutput5,
                   "minLinfmax"=moutput6))
}
