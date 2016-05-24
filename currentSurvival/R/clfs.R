
clfs <- function(data, maxx = NULL, com.est = TRUE, conf.int = FALSE, conf.int.level = NULL, no.iter = NULL, points = NULL, fig = TRUE, strat = FALSE, pvals = FALSE, pval.test = NULL)
{

### check input parameters:

## check input parameter com.est:
if (!is.logical(com.est)) {
  stop(paste("","Invalid logical parameter 'com.est'! Its value must be set to","TRUE or FALSE. The default value is TRUE.",sep="\n"))
}

## check input parameter conf.int:
if (!is.logical(conf.int)) {
  stop(paste("","Invalid logical parameter 'conf.int'! Its value must be set to","TRUE or FALSE. The default value is FALSE.",sep="\n"))
}

## check input parameter conf.int.level:
if (conf.int) {
  if (is.null(conf.int.level)) {
    conf.int.level <- 0.95
  } else { 
    if (!is.numeric(conf.int.level) || conf.int.level<0.9 || conf.int.level>0.99) {
      stop(paste("","Invalid numerical parameter 'conf.int.level'! Its value must be","in the range 0.9-0.99. The default value is 0.95.",sep="\n"))
    }
  }
} else {
  if (!is.null(conf.int.level)) {
    if (!is.numeric(conf.int.level) || conf.int.level<0.9 || conf.int.level>0.99) {
      stop(paste("","Invalid numerical parameter 'conf.int.level'! Its value must be","in the range 0.9-0.99. The default value is 0.95. However, if","you want to calculate confidence intervals, you must also set","'conf.int' to TRUE. If not, do not specify 'conf.int.level'.",sep="\n"))
    } else {
      stop(paste("","Parameter 'conf.int' is missing or set to FALSE! If you want to","calculate confidence intervals, you must set 'conf.int' to TRUE.","If not, do not specify 'conf.int.level'.",sep="\n"))
    }
  }
}

## check input parameter no.iter:
if (conf.int) {
  if (is.null(no.iter)) {
    no.iter <- 100
  } else { 
    if (!is.numeric(no.iter) || no.iter<10 || no.iter>10000) {
      stop(paste("","Invalid numerical parameter 'no.iter'! Its value must be in the","range 10-10000. The default value is 100.",sep="\n"))
    }
  }
} else {
  if (!is.null(no.iter)) {
    if (!is.numeric(no.iter) || no.iter<10 || no.iter>10000) {
      stop(paste("","Invalid numerical parameter 'no.iter'! Its value must be in the","range 10-10000. The default value is 100. However, if you want","to calculate confidence intervals, you must also set 'conf.int'","to TRUE. If not, do not specify 'no.iter'.",sep="\n"))
    } else {
      stop(paste("","Parameter 'conf.int' is missing or set to FALSE! If you want to","calculate confidence intervals, you must set 'conf.int' to TRUE.","If not, do not specify 'no.iter'.",sep="\n"))
    }
  }
}

## check input parameter fig:
if (!is.logical(fig)) {
  stop(paste("","Invalid logical parameter 'fig'! Its value must be set to TRUE","or FALSE. The default value is TRUE.",sep="\n"))
}

## check input logical parameter strat:
if (!is.logical(strat)) {
  stop(paste("","Invalid logical parameter 'strat'! Its value must be set to TRUE","or FALSE. The default value is FALSE.",sep="\n"))
}

## check input parameter pvals:
if (!is.logical(pvals)) {
  if (!strat || !conf.int) {
    stop(paste("","Invalid logical parameter 'pvals'! Its value must be set to TRUE","or FALSE. The default value is FALSE. If you want to calculate","p-values for the stratified CLFS estimates, you must also set","'strat' and 'conf.int' to TRUE. If not, do not specify 'pvals'.",sep="\n"))
  } else {
    stop(paste("","Invalid logical parameter 'pvals'! Its value must be set to TRUE","or FALSE. The default value is FALSE.",sep="\n"))
  }
}
if (pvals && !strat && !conf.int) {
  stop(paste("","Parameters 'strat' and 'conf.int' are missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'strat' and 'conf.int' to TRUE.","If not, do not specify 'pvals'.",sep="\n"))
}
if (pvals && !strat) {
  stop(paste("","Parameter 'strat' is missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'strat' to TRUE. If not, do not","specify 'pvals'.",sep="\n"))
}
if (pvals && !conf.int) {
  stop(paste("","Parameter 'conf.int' is missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'conf.int' to TRUE because the","computation of p-values is based on the estimation of","confidence intervals. If not, do not specify 'pvals'.",sep="\n"))
}

## check input parameter pval.test:
if (!is.null(pval.test)) {
  check <- switch(pval.test,
                  naive = 1,
                  log = 1,
                  loglog = 1 )
  if (is.null(check)) {
    if (!pvals || !strat || !conf.int) {
      stop(paste("","Invalid string parameter 'pval.test'! Its value must be set","to 'naive', 'log' or 'loglog'. The default value is 'loglog'.","If you want to calculate p-values for the stratified CLFS","estimates, you must also set 'pvals', 'strat' and 'conf.int'","to TRUE. If not, do not specify 'pvals.test'.",sep="\n"))
    } else { 
      stop(paste("","Invalid string parameter 'pval.test'! Its value must be set","to 'naive', 'log' or 'loglog'. The default value is 'loglog'.",sep="\n"))
    }
  } else {
    if (!pvals && !strat && !conf.int) {
      stop(paste("","Parameters 'pvals', 'strat' and 'conf.int' are missing or set","to FALSE! If you want to calculate p-values for the stratified","CLFS estimates, you must set 'pvals', 'strat' and 'conf.int'","to TRUE. If not, do not specify 'pvals.test'.",sep="\n"))
    }
    if (!pvals && !strat) {
      stop(paste("","Parameters 'pvals' and 'strat' are missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'pvals' and 'strat' to TRUE. If not,","do not specify 'pvals.test'.",sep="\n"))
    }
    if (!pvals && !conf.int) {
      stop(paste("","Parameters 'pvals' and 'conf.int' are missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'pvals' and 'conf.int' to TRUE. If not,","do not specify 'pvals.test'.",sep="\n"))
    }
    if (!strat && !conf.int) {
      stop(paste("","Parameters 'strat' and 'conf.int' are missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'strat' and 'conf.int' to TRUE. If not,","do not specify 'pvals.test'.",sep="\n"))
    }
    if (!pvals) {
      stop(paste("","Parameter 'pvals' is missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'pvals' to TRUE. If not, do not","specify 'pvals.test'.",sep="\n"))
    }
    if (!strat) {
      stop(paste("","Parameter 'strat' is missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'strat' to TRUE. If not, do not","specify 'pvals.test'.",sep="\n"))
    }
    if (!conf.int) {
      stop(paste("","Parameter 'conf.int' is missing or set to FALSE!","If you want to calculate p-values for the stratified CLFS","estimates, you must set 'conf.int' to TRUE because the","computation of p-values is based on the estimation of","confidence intervals. If not, do not specify 'pvals.test'.",sep="\n"))
    }
  }
} else {
  if (pvals) {
    pval.test <- "loglog"
  }
}



### data pre-processing:

if (strat) {
  stratf <- data[,ncol(data)] # separate the stratification factor
  data <- data[,-ncol(data)] # remove the stratification factor from the data matrix

  # check whether the data contain only numeric values:
  isfac <- array(0,ncol(data)) # allocation of a vector for the identification of factor variables among the data columns 
  for (i in 1:ncol(data)) {
    isfac[i] <- is.factor(data[,i])
  }
  if (sum(isfac)>0) {
    stop(paste("","Invalid input data! Data contain string variables (apart from","the stratification factor in which string values are allowed).",sep="\n"))
  }

} else {

  # check whether the data contain only numeric values:
  isfac <- array(0,ncol(data)) # allocation of a vector for the identification of factor variables among the data columns 
  for (i in 1:ncol(data)) {
    isfac[i] <- is.factor(data[,i])
  }
  if (sum(isfac)>0) {
    stop(paste("","Invalid input data! Data contain string variables.","If the last column of the data matrix is a stratification factor,","in which string values are allowed, and you want to calculate the","stratified CLFS estimates, you must set 'strat' to TRUE. If not,","all columns in the data matrix must be numeric variables of type","Integer.",sep="\n"))
  }
}


### other data controls:

# allocate a vector for error messages regarding input data:
error.messages <- NULL

# check whether the data matrix does not contain a column with only NA values:
if (sum(colSums(is.na(data))==nrow(data))>0) {
  if (sum(colSums(is.na(data))==nrow(data))==1) {
    warning(paste(paste("Column no.",c(1:ncol(data))[colSums(is.na(data))==nrow(data)],"was excluded from the data matrix as it contains"),"only NA values.",sep="\n"))
  } else {
    warning(paste(paste("Columns no.",paste(c(1:ncol(data))[colSums(is.na(data))==nrow(data)],collapse=", ")," were excluded from the data matrix"),"as they contain only NA values.",sep="\n"))
  }
  data <- data[,(colSums(is.na(data))<nrow(data))]
}

# check whether the times to events and the follow-up time are integers:
if (sum(!(abs(data[,-ncol(data)]-round(data[,-ncol(data)]))<.Machine$double.eps),na.rm=TRUE)>0) {
  error.messages <- c(error.messages,"All times to events as well as follow-up times must be in days","(i.e. integer values)!")
}

# check whether the times to events and the follow-up time contain only values higher than 0:
if (sum(data[,-ncol(data)]<=0,na.rm=TRUE)>0) {
  error.messages <- c(error.messages,"All times to events as well as follow-up times must be higher than 0!")
}


# allocate another vector for error messages regarding input data:
error.messages2 <- NULL

# check whether the event times are in ascending order:
E <- data[,1:(ncol(data)-2)]
for (i in 1:nrow(E)) {
  E[i,is.na(E[i,])] <- data[i,ncol(data)-1]
}
Ediff <- t(diff(t(E))) # compute a matrix of differences between the subsequent pairs of data columns
Ediff[is.na(data[,2:(ncol(data)-2)])] <- NA
if (sum(Ediff<=0,na.rm=TRUE)>0) {
  error.messages2 <- c(error.messages2,"The event times must be ascending and not equal for each patient!") # check whether all differences are >= 0
}
rm(list=c("E","Ediff"))


# remove patients who did not achieve the first disease remission:
if (strat) {
  stratf <- stratf[!is.na(data[,1])]
}
data <- data[!is.na(data[,1]),]


# check if there are no NA values in the censoring indicator: 
if(sum(is.na(data[,ncol(data)]))>0) {
  warning(paste("","NA values were detected in the censoring indicator.",paste("Number of patients excluded from the analysis due to NAs:",sum(is.na(data[,ncol(data)]))),sep="\n"))
  if (strat) {
    stratf <- stratf[!is.na(data[,ncol(data)])]
  }
  data <- data[!is.na(data[,ncol(data)]),]
}

# check if there are no NA values in the follow-up time: 
if(sum(is.na(data[,ncol(data)-1]))>0) {
  warning(paste("","NA values were detected in the follow-up time.",paste("Number of patients excluded from the analysis due to NAs:",sum(is.na(data[,ncol(data)-1]))),sep="\n"))
  if (strat) {
    stratf <- stratf[!is.na(data[,ncol(data)-1])]
  }
  data <- data[!is.na(data[,ncol(data)-1]),]
}

# check whether the censoring indicator has only values 0 and 1:
if ((sum(data[,ncol(data)]==0) + sum(data[,ncol(data)]==1)) != nrow(data)) {
  error.messages2 <- c(error.messages2,"Invalid censoring indicator; it must be 0 (censored) or 1 (dead)","for each patient!") 
}

# check whether the follow-up time is higher than the preceding times to events:
fup <- array(0,nrow(data))
for (i in 1:nrow(data)) {
  if (sum(!is.na(data[i,1:(ncol(data)-2)]))>0) {
    fup[i] <- (max(data[i,1:(ncol(data)-2)],na.rm=TRUE)>data[i,ncol(data)-1]) # the follow-up time must be higher or equal to individual event times
  }
}
if (sum(fup)>0) {
  error.messages2 <- c(error.messages2,"Invalid follow-up time; it must be higher than all event times","for each patient!") 
}

# check the number of levels for stratification:
if (strat) {
  if (length(levels(as.factor(stratf)))<2) {
    error.messages2 <- c(error.messages2,"Stratification factor must have at least 2 levels!")
  } else {
    if (length(levels(as.factor(stratf)))>8) {
      error.messages2 <- c(error.messages2,"Stratification factor must have no more than 8 levels!")
    }
  }
}

# check if there are no NA values in the stratification factor: 
if (strat) {
  if(sum(is.na(stratf))>0) {
    warning(paste("","NA values were detected in the stratification factor.",paste("Number of patients excluded from the analysis due to NAs:",sum(is.na(stratf))),sep="\n"))
    data <- data[!is.na(stratf),]
    stratf <- stratf[!is.na(stratf)]
  }
}

# create final error message:
if (strat) {
  if (!is.null(error.messages)) {
    if (!is.null(error.messages2)) {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages,collapse="\n"),paste(error.messages2,collapse="\n"),"As the parameter 'strat' is set to TRUE, check whether the last","three data columns contain follow-up times, censoring indicators,","and stratification levels, respectively.",sep="\n"))
    } else {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages,collapse="\n"),sep="\n"))
    }
  } else {
    if (!is.null(error.messages2)) {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages2,collapse="\n"),"As the parameter 'strat' is set to TRUE, check whether the last","three data columns contain follow-up times, censoring indicators,","and stratification levels, respectively.",sep="\n"))
    }
  }
} else {
  if (!is.null(error.messages)) {
    if (!is.null(error.messages2)) {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages,collapse="\n"),paste(error.messages2,collapse="\n"),"As the parameter 'strat' is missing or set to FALSE, check","whether the last two data columns contain follow-up times and","censoring indicators, respectively.",sep="\n"))
    } else {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages,collapse="\n"),sep="\n"))
    }
  } else {
    if (!is.null(error.messages2)) {
      stop(paste("","Following problem(s) were found in the data matrix:",paste(error.messages2,collapse="\n"),"As the parameter 'strat' is missing or set to FALSE, check","whether the last two data columns contain follow-up times and","censoring indicators, respectively.",sep="\n"))
    }
  }
}



### check input parameters maxx and points:

# check input parameter maxx and convert it to days:
LastContact.clfs <- data[,ncol(data)-1] - data[,1]
if (is.null(maxx)) {
  maxx <- max(LastContact.clfs,na.rm=TRUE)
} else { 
  if (!is.numeric(maxx) || maxx<1 || maxx>(max(LastContact.clfs,na.rm=TRUE)/365)) {
    stop(paste("","Invalid numerical parameter 'maxx'! It must be in range from","1 year to the maximum follow-up time except the time from therapy","initiation to achievement of the first disease remission. The","default value is the maximum follow-up time except the time from","therapy initiation to achievement of the first disease remission.",sep="\n"))
  } else {
    maxx <- floor(maxx*365)
  }
}

# check input parameter points:
if (is.null(points)) {
  points <- seq(12,floor(maxx/(365/12)),12)
} else { 
  if (!is.vector(points,mode="numeric") || sum(points<0)>0 || sum(points>floor(maxx/(365/12)))>0) {
    stop(paste("","Invalid numerical vector 'points'! Its values must be in range","from 0 months to the maximum follow-up time except the time from","therapy initiation to achievement of the first disease remission.","The default is a vector of 0, 12, 24, ..., floor(maxx/(365/12))","months.",sep="\n"))
  } else {
    points <- sort(points)
    if (sum(points==0)>0) {
      points <- points[-1] # if 0 is included in the time points, it is removed because it will be added later automatically
    }
  }
}


### call functions clfs.strat or clfs.nostrat:

if (strat) {
  clfs <- clfs.strat(data, stratf, maxx, com.est, conf.int, conf.int.level, no.iter, points, fig, pvals, pval.test)
} else {
  clfs <- clfs.nostrat(data, maxx, com.est, conf.int, conf.int.level, no.iter, points, fig)
}

}
