cat("###################################################################
########################### class matrix ##########################
################### Imputations des manquantes ####################
###################################################################\n")

#################
### global mean
###  - pour monotones et intermitentes

imput_mean_col <- function(col,force=FALSE){
    if(all(is.na(col))){
        if(is.numeric(force)){
            warning("[Imputation:crossMean] There is only NA on this column, FORCED imputation\n")
            return(rep(force,length(col)))
        }else{
            warning("[Imputation:crossMean] There is only NA on this column, impossible to impute\n")
            return(col)
        }
    }else{
        if(all(!is.na(col))){return(col)}else{}
    }
    col[is.na(col)] <- mean(col,na.rm=TRUE)
    return(col)
}

imput_crossMean <- function(longData,force=FALSE){
    if(identical(force,TRUE)){force <- mean(longData,na.rm=TRUE)}else{}
    return(apply(longData,2,imput_mean_col,force))
}


#################
### global median
###  - pour monotones et intermitentes

imput_median_col <- function(col,force=FALSE){
   if(all(is.na(col))){
        if(is.numeric(force)){
            warning("[Imputation:crossMedian] There is only NA on this column, FORCED imputation\n")
            return(rep(force,length(col)))
        }else{
            warning("[Imputation:crossMedian] There is only NA on this column, impossible to impute\n")
            return(col)
        }
    }else{
        if(all(!is.na(col))){return(col)}else{}
    }
    col[is.na(col)] <- median(col,na.rm=TRUE)
    return(col)
}

imput_crossMedian <- function(longData,force=FALSE){
    if(identical(force,TRUE)){force <- median(longData,na.rm=TRUE)}else{}
    return(apply(longData,2,imput_median_col,force))
}


#################
### hot deck
###  - pour monotones et intermitentes

imput_hotDeck_col <- function(col,force=FALSE){
    if(all(is.na(col))){
        if(is.numeric(force)){
            warning("[Imputation:crossHotDeck] There is only NA on this column, FORCED imputation\n")
            return(force[floor(runif(length(col),min=1,max=length(force)+1))])
        }else{
            warning("[Imputation:crossHotDeck] There is only NA on this column, impossible to impute\n")
            return(col)
        }
    }else{
        if(all(!is.na(col))){return(col)}else{}
    }
    missing <- is.na(col)
    nbPos <- sum(!missing)
    alea <- floor(runif(length(col)-nbPos,min=1,max=nbPos+1))
    col[missing] <- sapply(alea,function(x)col[cumsum(!missing)==x & !missing])
    return(col)
}

imput_crossHotDeck <- function(longData,force=FALSE){
    if(identical(force,TRUE)){force <- na.omit(as.numeric(longData))}else{}
    return(apply(longData,2,imput_hotDeck_col,force))
}


