

#'Filling the last cell
#'@param x a data vector.
fill_last <- function(x) {
  na_idx <- which(is.na(x))
  unique_elements <- unique(x[-na_idx])
  set_diff <- unique_elements[length(unique_elements)]
  x[na_idx] <- set_diff
  x
}

#'Data pre-processing for analysis with stochastic process model methodology.
#'@param x A path to the table with follow-up oservations (longitudinal study). Formats: csv, sas7bdat
#'@param y A path to the table with vital statistics (mortality). File formats: csv, sas7bdat
#'@param col.id A name of column containing subject ID. 
#'This ID should be the same in both longdat and vitstat tables.
#'If not provided, the first column in the x and y will be used by default.
#'@param col.status A name of the column containing status variable (0/1 which indicate alive/dead). 
#'If not provided - then the column #2 from the vital statistics dataset will be used.
#'@param col.age A name of age column (also called 't1'). 
#'If not provided then the 3rd column from the longitudinal dataset (x) will be used.
#'@param col.age.event A name of 'event' column.
#'The event column indicates a time when the even occured (e.g. system failure).
#'If not provided then the 3rd column from the vital statistics dataset will be used.
#'@param covariates A list of covariates. 
#'If covariates not provided, then all columns from longitudinal table having index > 3 will be used as covariates. 
#'@param interval A number of breaks between observations for discrete model. Default = 1 unit of time.
#'@param verbose A verbosing output indicator. Default=FALSE.
#'@return A list of two elements: first element contains a preprocessed data for continuous model, with arbitrary intervals between observations  and 
#'second element contains a prepocessed data table for a discrete model (with constant intervals between observations).
#'@examples \dontrun{ 
#'library(stpm) 
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#'					   y=system.file("data","vitstat.csv",package="stpm"))
#'head(data[[1]])
#'head(data[[2]])
#'}
prepare_data <- function(x, y, 
                         col.id=NULL, 
                         col.status=NULL,
                         col.age=NULL, 
                         col.age.event=NULL, 
                         covariates=NULL, 
                         interval=1, 
                         verbose=FALSE) {
  
  
  if(file_ext(x) == "csv") {
    longdat <- read.csv(x)
  } else if(file_ext(x) == "sas7bdat") {
    longdat <- read.sas7bdat(x)
  } else {
  	stop(paste(x, ":", "unknown file format, it must be csv or sas7bdat."))
  }
  
  if(file_ext(y) == "csv") {
    vitstat <- read.csv(y)
  } else if(file_ext(y) == "sas7bdat") {
    vitstat <- read.sas7bdat(y)
  } else {
  	stop(paste(y, ":", "unknown file format, it must be csv or sas7bdat."))
  }
  
  # Parsing input parameters in order to check for errors:
  if( !is.null(col.status) ) {
    if( !(col.status %in% colnames(vitstat)) ) {
      stop(paste("Status column",col.status, "not found in vitstat table. Aborting."))
    }
    col.status.ind <- grep(paste("\\b", col.status, "\\b", sep=""), colnames(vitstat))
  } else if(is.null(col.status)) {
    col.status.ind <- 2
  }
  
  if( !is.null(col.id) ) { 
    if( !(col.id %in% colnames(vitstat)) || !(col.id %in% colnames(longdat)) ) {
      stop(paste("ID column",col.id, "not found in vitstat and/or longdat tables. Aborting."))
    }
    col.id.ind <- grep(paste("\\b", col.id, "\\b", sep=""), colnames(vitstat))
  } else if(is.null(col.id)) {
    col.id.ind <- 1
  }
  
  if( !is.null(col.age) ) {
    if( !(col.age %in% colnames(longdat)) ) {
      stop(paste("Age column",col.age, "not found in longdat table. Aborting."))
    }
    col.age.ind <- grep(paste("\\b", col.age, "\\b", sep=""), colnames(longdat))
  } else if(is.null(col.age)) {
    col.age.ind <- 3
  } 
  
  if( !is.null(col.age.event) ) { 
    if( !(col.age.event %in% colnames(vitstat)) ) {
      stop(paste("Event column",col.age.event, "not found in vitstat table. Aborting."))
    }
    col.age.event.ind <- grep(paste("\\b", col.age.event, "\\b", sep=""), colnames(vitstat))
  } else if(is.null(col.age.event)) {
    col.age.event.ind <- 3
  }
  
  if(!is.null(covariates)) {
    col.covar.ind <- c()
    for(c in covariates) {
      if( !(c %in% colnames(longdat)) ) {
        stop(paste("Covariate",c, "not found. Aborting."))
      }
      col.covar.ind <- c(col.covar.ind, grep(paste("\\b", c, "\\b", sep=""), colnames(longdat)))
    }
  } else if(is.null(covariates)) {
    col.covar.ind <- 4:dim(longdat)[2]
  }
  
  if((interval == 0) || (interval > 1)) {
    interval <- 1
  }
  
  #-----------Done parsing imput parameters---------------------#
  # First time of data pre-processing:
  longdat <- longdat[which(!is.na(longdat[ , col.age.ind])),]
  
  # Prepare data for continuous optimisation:
  data_cont <- prepare_data_cont(longdat, vitstat, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose)
  
  # Prepare data for fast discrete optimization:
  data_discr <- prepare_data_discr(longdat, vitstat, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose)
  
  list(model.continuous=data_cont, model.discrete=data_discr)
}

#'Prepares continuouts-time dataset.
#'@param longdat a longitudinal study dataset.
#'@param vitstat vital (mortality) statistics.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
prepare_data_cont <- function(longdat, 
                              vitstat, 
                              col.status.ind, 
                              col.id.ind, 
                              col.age.ind, 
                              col.age.event.ind, 
                              col.covar.ind, 
                              verbose) {
  
  # Split records by ID:
  prep.dat <- matrix(ncol=(4+2*length(col.covar.ind)),nrow=0)
  splitted <- split(longdat, longdat[ , col.id.ind])
  vitstat.splitted <- split(vitstat, vitstat[ , col.id.ind])
  
  for(iii in 1:length(splitted)) {
    nrows <- length(splitted[[iii]][ , col.id.ind])
    id <- splitted[[iii]][ , col.id.ind]
    case <- rep(0, nrows)
    case[nrows] <- vitstat.splitted[[iii]][, col.status.ind]
    t1 <- splitted[[iii]][ , col.age.ind]
    t2 <- c(splitted[[iii]][ , col.age.ind][-1], vitstat.splitted[[iii]][ , col.age.event.ind])
    
    tmp.frame <- cbind(id, case, t1, t2)
    # Adding covariates:
    for(ind in col.covar.ind) {
      tmp.frame <- cbind(tmp.frame, 
                         splitted[[iii]][, ind], 
                         c(splitted[[iii]][, ind][-1], NA))
      
    }
    prep.dat <- rbind(prep.dat, tmp.frame)
  }
  
  
  #prep.dat <- prep.dat[rowSums( matrix(is.na(prep.dat[,5:dim(prep.dat)[2]]), ncol=2*length(covariates),byrow=T)) !=2*length(covariates),]
  prep.dat <- prep.dat[which(is.na(prep.dat[,4])==FALSE),]
  
  if(verbose) {
    head(prep.dat)  
  }
  
  ans_final <- prep.dat
  if(length(which(is.na(prep.dat[,5:dim(prep.dat)[2]]) == TRUE)) > 0) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tmp_ans <- mice(prep.dat[,5:dim(prep.dat)[2]], printFlag=ifelse(verbose, TRUE, FALSE),m = 2, maxit=2)
    ans1 <- complete(tmp_ans)
    #ans_final <- cbind(prep.dat[,1:3], ans1)
    ans_final <- cbind(prep.dat[,1:4], ans1)
  }
  
  if(verbose)
    cat("Making final table...\n")
  
  # Database should be in appropriate format:
  for(i in 1:(dim(ans_final)[1])) {
    if(ans_final[i,2] > 1) {
      ans_final[i,2] <- 1
    }
  }
  
  # Finalizing:
  ans_final <- ans_final[which(ans_final[,3] != ans_final[,4]),] # t1 must be different from t3
  ans_final <- ans_final[which(ans_final[,3] < ans_final[,4]),] # t1 must be less than t3
  # t1 must be equal t3 on previous step, if status = 0 and id is the same
  for(i in 2:dim(ans_final)[1]) {
    if((ans_final[i,3] != ans_final[(i-1),4]) & (ans_final[i,2] == 0) & (ans_final[i,1] == ans_final[(i-1),1])) {
      ans_final[i,3] <- ans_final[(i-1),4]
    }
  }
  
  colnames(ans_final) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(col.covar.ind), function(n) {c(names(longdat)[col.covar.ind[n]], 
                                                                                                  paste(names(longdat)[col.covar.ind[n]],".next",sep=""))} )) )
  ans_final
  
}

#'Prepares discrete-time dataset.
#'@param longdat a longitudinal study dataset.
#'@param vitstat vital (mortality) statistics.
#'@param interval interval between observations.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
prepare_data_discr <- function(longdat, vitstat, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose) {
  
  # Interpolation
  dt <- interval
  tt <- matrix(nrow=0, ncol=4)
  par <- matrix(nrow=0, ncol=length(col.covar.ind))
  
  # Split records by ID:
  splitted <- split(longdat, longdat[, col.id.ind])
  vitstat.splitted <- split(vitstat, vitstat[, col.id.ind])
  #iii <- 1
  # For each particular person's record:
  for(iii in 1:length(splitted)) {
    if(!is.na(vitstat.splitted[[iii]][ , col.age.event.ind]) & !is.na(vitstat.splitted[[iii]][ , col.status.ind]) ) {
      if(verbose) {
        print(iii)
      }
      id <- splitted[[iii]][ , col.id.ind][1]
      nrows <- (tail(splitted[[iii]][ , col.age.ind], n=1) - splitted[[iii]][ , col.age.ind][1])/dt + 1
      # Perform approximation:
      t1.approx <- matrix(ncol=4, nrow=nrows)
      t1.approx[,1] <- id
      t1.approx[,2] <- 0
      t1.approx[nrows,2] <- vitstat.splitted[[iii]][ , col.status.ind][1] #Last value
      t1.approx[,3] <- seq(splitted[[iii]][ , col.age.ind][1], splitted[[iii]][ , col.age.ind][length(splitted[[iii]][ , col.age.ind])], by=dt)
      if(nrows > 1) {
        t1.approx[,4] <- c(t1.approx[,3][2:nrows], vitstat.splitted[[iii]][ , col.age.event.ind][1])
      } else {
        t1.approx[,4] <- vitstat.splitted[[iii]][ , col.age.event.ind][1]
      }
      
      tt <- rbind(tt,t1.approx)
      par1.approx <- matrix(ncol=length(col.covar.ind), nrow=nrows, NA)
      
      j <- 1
      for(ind in col.covar.ind) {
        #ind <- col.covar.ind[j]
        if ( (length(splitted[[iii]][, ind]) > 1) & (length(which(!is.na(splitted[[iii]][, ind]))) > 0) ) {
          if(length(which(!is.na(splitted[[iii]][, ind]))) == 1) {
            splitted[[iii]][, ind] <- fill_last(splitted[[iii]][, ind])
          }
          # Fill NAs by linear approximation with approx():
          nn <- length(splitted[[iii]][, ind])
          splitted[[iii]][, ind] <- approx(splitted[[iii]][, ind],n=nn)$y
          par1.approx[,j] <-  approx(splitted[[iii]][, ind], n=nrows)$y
        }
        
        j <- j + 1
        
      }
      par <- rbind(par,par1.approx)
    }
  }
  
  ans=cbind(tt,par)
  colnames(ans) <- c("id", "case", "t1", "t2", names(longdat)[col.covar.ind])
  
  ans <- ans[rowSums( matrix(is.na(ans[,5:dim(ans)[2]]), ncol=length(col.covar.ind),byrow=T)) !=length(col.covar.ind),]
  
  ans_final <- ans
  if(length(which(is.na(ans[,5:dim(ans)[2]]) == TRUE)) > 0) {
    if(verbose)
      cat("Filing missing values with multiple imputations:\n")
    
    tmp_ans <- mice(ans[,5:dim(ans)[2]], printFlag=ifelse(verbose, TRUE, FALSE), m = 2, maxit = 2)
    ans1 <- complete(tmp_ans)
    ans_final <- cbind(ans[,1:4], ans1)
  }
  
  if(verbose)
    cat("Making final table...\n")
  ndim <- length(col.covar.ind)
  averages = matrix(nrow=1,ncol=length(col.covar.ind))
  
  dat <- ans_final[,1] #pid
  dat <- cbind(dat, ans_final[,2]) #sta (outcome)
  dat <- cbind(dat, ans_final[,3]) #tt1 (t1)
  dat <- cbind(dat, ans_final[,4]) #tt3 (t2)
  
  
  j <- 0
  i <- 0
  for(i in 0:(length(col.covar.ind)-1)) {
    dat <- cbind(dat, ans_final[,(5+i)]) 
    dat[2:dim(dat)[1],(5+j)] <- dat[1:(dim(dat)[1]-1),(5+j)]
    dat <- cbind(dat, ans_final[,(5+i)]) 
    averages[1,(i+1)] = dat[1,(5+j)]
    j <- j + 2
  }
  
  # Database should be in appropriate format:
  pid=dat[1,1]
  for(i in 1:(dim(dat)[1]-1)) {
    if(dat[i,1] != pid) {
      for(ii in seq(0,(ndim-1),2)) {
        dat[(i+1),(5+ii)] = dat[i,(6+ii)]
      }
      pid = dat[i,1]
    }
    if(dat[i,2] > 1) {
      dat[i,2] <- 1
    }
  }
  
  colnames(dat) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(col.covar.ind), function(n) {c(names(longdat)[col.covar.ind[n]], 
                                                                                                  paste(names(longdat)[col.covar.ind[n]],".next",sep="")
                                                                                                  )} 
                                                             )
                                                      ) 
                     )
  rownames(dat) <- 1:dim(dat)[1]
  dat
}

