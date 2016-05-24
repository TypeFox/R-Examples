mean_centering<-function(Y,b){
################################################################################

if(missing(Y)){stop("Argument 'Y' missing, with no default\n")}

if(missing(b)){stop("Argument 'b' missing, with no default\n")}

if(class(Y)!='matrix'){stop("'Y' must be of class 'matrix'\n")}

if(class(b)!='factor'){stop("'b' must be of class 'factor'\n")}

if(any(is.na(b))){stop("NA values are not allowed in 'b'\n")}

if(length(b)!=nrow(Y)){stop("length(b) is different from nrow(Y)\n")}

if(any(apply(Y,2,mode)!='numeric')){stop('Array expression columns contain non-numeric values!\n')}

################################################################################

    Y1<-t(Y)
    dd <- model.matrix(~factor(b) - 1)
    
#####################################

    naY1 <- is.na(Y1)
    nsamples <- (!naY1)%*%dd
    Y1[naY1] <- 0
    Y1sum <- Y1%*%dd
    Y1bar <- Y1sum/nsamples
    mean_matrix<-Y1bar %*% t(dd)

    EAdj <- Y1 - mean_matrix
    return(t(EAdj))
}
