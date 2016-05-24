## syntax:
##  "col1, col2, col3 ~ row1, row2, row3"
## NOTE: variable names CANNOT contain commas,  tildes, or whitespace
##     
"get.mat" <- function(fstring, data){


  
  if (length(grep(",(\\s)*~", fstring)) > 0){
    stop("formula string has ',' adjacent to '~'")
  }
  
  if (length(grep("~(\\s)*,", fstring)) > 0){
    stop("formula string has '~' adjacent to ','")
  }


  
  data <- as.data.frame(data)
  
  ## check to see if any variables in data have commas in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep(",", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains ','") 
  }

  ## check to see if any variables in data have tildes in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep("~", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains '~'") 
  }

  ## check to see if any variables in data have whitespace in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep(" ", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains ' '") 
  }

  
  ## separate into row and column vars (split on ~)
  fstring.split <- strsplit(fstring, "~")
  colvars <- fstring.split[[1]][1]
  rowvars <- fstring.split[[1]][2]

  ## separate into individual variables by splitting on ',' and
  ## getting rid of whitespace
  colvars <- strsplit(colvars, ",")
  rowvars <- strsplit(rowvars, ",")
  ncv <- length(colvars[[1]])
  nrv <- length(rowvars[[1]])
  colvars.vec <- NULL
  rowvars.vec <- NULL
  for (i in 1:ncv){
    colvars[[1]][i] <- gsub(" ", "", colvars[[1]][i])
    colvars.vec <- c(colvars.vec, colvars[[1]][i])
  }
  for (i in 1:nrv){
    rowvars[[1]][i] <- gsub(" ", "", rowvars[[1]][i])
    rowvars.vec <- c(rowvars.vec, rowvars[[1]][i])
  }
  
  
  colmat <- NULL
  colmat.names <- NULL
  rowmat <- NULL
  rowmat.names <- NULL
  for (i in 1:ncv){
    if (length(grep(paste("^", colvars.vec[i], sep=""), colnames(data))) > 1){
      stop(paste("dataframe contains multiple variables starting with '",
                 colvars.vec[i], "'", sep="")) 
    }
    holder <- eval(parse(text=paste("data$", colvars.vec[i], sep="")))
    colmat <- cbind(colmat, holder)
    holder <- colnames(data)[grep(paste("^", colvars.vec[i], sep=""),
                                  colnames(data))]
    colmat.names <- c(colmat.names, holder)
  }
  for (i in 1:nrv){
    if (length(grep(paste("^", rowvars.vec[i], sep=""), colnames(data))) > 1){
      stop(paste("dataframe contains multiple variables starting with '",
                 rowvars.vec[i], "'", sep="")) 
    }
    holder <- eval(parse(text=paste("data$", rowvars.vec[i], sep="")))
    rowmat <- cbind(rowmat, holder)
    holder <- colnames(data)[grep(paste("^", rowvars.vec[i], sep=""),
                                  colnames(data))]
    rowmat.names <- c(rowmat.names, holder)
  }

  colnames(colmat) <- colmat.names
  colnames(rowmat) <- rowmat.names
  

  return(list(colmat=colmat, rowmat=rowmat))
  
}

"get.Xmatrix" <- function(covstring,data,NN.cellnames,
                          rownames.local,colnames.local,standardize)
{
  # Make into a list if not already:
  if (!is.list(covstring)){
    covstring <- list(covstring)
  }

  for (j in 1:length(covstring)){

    # Error check each element of the list:

    if (length(grep(",(\\s)*~", covstring[[j]])) > 0){
      stop("formula string has ',' adjacent to '~'")
    }
  
    if (length(grep("~(\\s)*,", covstring[[j]])) > 0){
      stop("formula string has '~' adjacent to ','")
    }

    # New checks: 
    if (length(grep(",(\\s)*=", covstring[[j]])) > 0){
      stop("formula string has ',' adjacent to '='")
    }
  
    if (length(grep("=(\\s)*,", covstring[[j]])) > 0){
      stop("formula string has '=' adjacent to ','")
    }

    if ( (length(grep("=",covstring[[j]]))>0) && (length(grep("~",covstring[[j]])))){
      stop("formula string has both '~' and '='")
    }

  }
  
  data <- as.data.frame(data)

  ## check to see if any variables in data have commas in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep(",", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains ','") 
  }

  ## check to see if any variables in data have tildes in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep("~", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains '~'") 
  }

  ## check to see if any variables in data have whitespace in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep(" ", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains ' '") 
  }

  # New check:

  ## check to see if any variables in data have equals in their names
  ## abort if true
  data.names <- colnames(data)
  if (length( grep("=", data.names) ) > 0){
    stop("dataframe contains a variable whose name contains '='") 
  }

  # For each element of the list:
  response.vars.vec   <- list(NULL)
  covariate.vars.vec  <- list(NULL)
  equals.or.tilde.vec <- list(NULL)

  for (j in 1:length(covstring)){

    # Empty holders to start with:
    response.vars.vec[[j]]   <- NULL
    covariate.vars.vec[[j]]  <- NULL
    equals.or.tilde.vec[[j]] <- NULL

    ## Separate into response and covariates (split on either ~ or =)
    ## Determine if this has a tilde or an equals as the separator

    if (length(grep("~",covstring[[j]]))){

      ## Expression is a ~ formula:
      covstring.split  <- strsplit(covstring[[j]], "~")
      equals.or.tilde.vec[[j]] <- "tilde"

    } else {

      ## Expression is an = formula:
      covstring.split  <- strsplit(covstring[[j]], "=")
      equals.or.tilde.vec[[j]] <- "equals"

    }

    response.vars  <- covstring.split[[1]][1]
    covariate.vars <- covstring.split[[1]][2]

    ## separate into individual variables by splitting on ',' and
    ## getting rid of whitespace
    response.vars  <- strsplit(response.vars, ",")
    covariate.vars <- strsplit(covariate.vars, ",")
    ncv <- length(response.vars[[1]])
    nrv <- length(covariate.vars[[1]])
    tmp.response   <- NULL
    tmp.covariates <- NULL

    for (i in 1:ncv){
      response.vars[[1]][i]  <- gsub(" ", "", response.vars[[1]][i])
      tmp.response <- c(tmp.response, response.vars[[1]][i])
    }
    for (i in 1:nrv){
      covariate.vars[[1]][i]  <- gsub(" ", "", covariate.vars[[1]][i])
      tmp.covariates <- c(tmp.covariates, covariate.vars[[1]][i])
    }

    response.vars.vec[[j]]  <- tmp.response
    covariate.vars.vec[[j]] <- tmp.covariates
  
  }

  ## Check that the response names match those supplied in the NN formula:

##cat("NN.cellnames:\n")
##print(NN.cellnames)
##cat("\n")

  for (j in 1:length(covstring)){

##cat("Checking the following vector of responses:\n")
##print(response.vars.vec[[j]])
##cat("Each element must match one of the elements in NN.cellnames.\n")

    for (k in 1:length(response.vars.vec[[j]])){    
      if (length(grep(response.vars.vec[[j]][k],NN.cellnames))==0){
        cat(paste("Error: response variable '",response.vars.vec[[j]][k],"' not found in:\n",sep=""))
        print(NN.cellnames)
        stop("invalid response variable in covstring")
      }
    }
  }

  ## Intercepts are always unique and untouched:
  retlist <- list(NULL)
  for (k in 1:nrow(data)){
    retlist[[k]] <- diag(1,length(NN.cellnames))
  }

  ## Take care of the names:
  names(retlist) <- rownames(data)


  ## Intercepts are always included first:
  retlist.colnames <- paste("mu", rep(rownames.local, each=(length(colnames.local)-1)),
                            rep(colnames.local[-length(colnames.local)],
                            length(rownames.local)), sep=".")

  nprecincts <- nrow(data)

  ## Take care of the equals equations first:
 
  for (j in 1:length(covstring)){
    if (equals.or.tilde.vec[[j]]=="equals"){

############################
##### Equals formulae ######
############################

      ## Go covariate-by-covariate:
      tmp.covariates <- covariate.vars.vec[[j]]
      tmp.response   <- response.vars.vec[[j]]

##cat(paste("Formula ",j," is an equals formula:\n",sep=""))
##cat(paste("Covariates in formula ",j,":\n",sep=""))
##print(tmp.covariates)
##cat(paste("Responses in formula ",j,":\n",sep=""))
##print(tmp.response)

      for (i in 1:length(tmp.covariates)){

        if (length(grep(paste("^", tmp.covariates[i], sep=""), colnames(data))) > 1){
          ## Duplicate match:
          stop(paste("dataframe contains multiple variables starting with '",
                     tmp.covariates[i], "'", sep="")) 
        }

        ## Extract the values across precincts:
        all.precincts <- eval(parse(text=paste("data$", tmp.covariates[i], sep="")))

        ## Find out which rows have this affects:
        affected.rows <- rep(FALSE,length(NN.cellnames))

        ## handle the variable names
        tmp.varname <- tmp.covariates[i]

        ## remember to omit the last column: 
        rownamevec <- rep(rownames.local,times=length(colnames.local)-1)
        colnamevec <- rep(colnames.local[1:(length(colnames.local)-1)],each=length(rownames.local))

        for (k in 1:length(tmp.response)){
          for (l in 1:length(NN.cellnames)){

            ## NOTE: This is bad, NEED TO PREVENT PARTIAL MATCHING HERE (ALLOWED CURRENTLY)
            try.k <- grep(tmp.response[k],NN.cellnames[l])
            if (length(try.k)>0){

              ## Cell k has equal covariate:
              affected.rows[l] <- TRUE

              ## Append the cell name:
              tmp.varname <- paste(tmp.varname,colnamevec[l],rownamevec[l],sep=".")

            }
 
          }
        }

##cat("Affected rows:\n")
##print(affected.rows)

        ## Should have at least one affected row o/w something is amiss:
        if (sum(affected.rows)==0)
          stop("Error: Invalid response(s) that eluded error check(s)")

        ## Create a matrix containing the columns of precinct design matrices
        ## as its rows:
        tmp.matrix <- matrix(0,nrow=nprecincts,ncol=length(NN.cellnames))

        ## In each of the affected columns we put in the covariates:
        for (k in 1:length(NN.cellnames)){
          if (affected.rows[k]==TRUE){
            tmp.matrix[,k] <- all.precincts
          }
        }

        ## Now transfer the temporary matrix into the full design matrices:
        for (k in 1:nprecincts){
          retlist[[k]] <- cbind(retlist[[k]],matrix(tmp.matrix[k,],ncol=1))
        }

        ## handle the variable name (only one column here):
        retlist.colnames <- c(retlist.colnames,tmp.varname)


      } ## end: for (i in 1:length(tmp.covariates)){...}

    } else {  ## end: if (equals.or.tilde.vec[[j]]=="equals"){...}

############################
##### Tilde formulae #######
############################

      ## Go covariate-by-covariate:
      tmp.covariates <- covariate.vars.vec[[j]]
      tmp.response   <- response.vars.vec[[j]]

##cat(paste("Formula ",j," is an equals formula:\n",sep=""))
##cat(paste("Covariates in formula ",j,":\n",sep=""))
##print(tmp.covariates)
##cat(paste("Responses in formula ",j,":\n",sep=""))
##print(tmp.response)

      for (i in 1:length(tmp.covariates)){

        if (length(grep(paste("^", tmp.covariates[i], sep=""), colnames(data))) > 1){
          ## Duplicate match:
          stop(paste("dataframe contains multiple variables starting with '",
                     tmp.covariates[i], "'", sep="")) 
        }

        ## Extract the values across precincts:
        all.precincts <- eval(parse(text=paste("data$", tmp.covariates[i], sep="")))

        ## Find out which rows have this affects:
        affected.rows <- rep(FALSE,length(NN.cellnames))

        ## handle the variable names
        tmp.varnames <- NULL

        ## remember to omit the last column: 
        rownamevec <- rep(rownames.local,times=length(colnames.local)-1)
        colnamevec <- rep(colnames.local[1:(length(colnames.local)-1)],each=length(rownames.local))

        for (k in 1:length(tmp.response)){
          for (l in 1:length(NN.cellnames)){

            ## NOTE: This is bad, NEED TO PREVENT PARTIAL MATCHING HERE (ALLOWED CURRENTLY)
            try.k <- grep(tmp.response[k],NN.cellnames[l])
            if (length(try.k)>0){

              ## Cell k has equal covariate:
              affected.rows[l] <- TRUE

              ## Append the cell name:
              tmp.varnames <- c(tmp.varnames,paste(tmp.covariates[i],colnamevec[l],rownamevec[l],sep="."))

            }
 
          }
        }

##cat("Affected rows:\n")
##print(affected.rows)

        ## Should have at least one affected row o/w something is amiss:
        if (sum(affected.rows)==0)
          stop("Error: Invalid response(s) that eluded error check(s)")

        ## handle the variable names:
        retlist.colnames <- c(retlist.colnames,tmp.varnames)

        ## In each of the affected columns we put in the covariates:
        for (k in 1:length(NN.cellnames)){
          if (affected.rows[k]==TRUE){

            ## Create a new matrix containing the new column:
            tmp.matrix <- matrix(0,nrow=nprecincts,ncol=length(NN.cellnames))
            tmp.matrix[,k] <- all.precincts
        
            ## Now transfer the temporary matrix into the full design matrices:
            for (l in 1:nprecincts){
              retlist[[l]] <- cbind(retlist[[l]],matrix(tmp.matrix[l,],ncol=1))
            }

          } ## end: if (affected.rows[k]==TRUE){...}
        } ## end: for (k in 1:length(NN.cellnames){...}

      } ## end: for (i in 1:length(tmp.covariates)){...}
    
    } ## end: ## if (equals.or.tilde.vec[[j]]=="equals"){...} else {...}

  } ## end: for (j in 1:length(covstring)){...}

  ## Put column names in place:
  for (i in 1:nprecincts){
    colnames(retlist[[i]]) <- retlist.colnames
  }

  ### STANDARDIZE THE COVARIATES BY DEFAULT ###


  return(retlist)

}


