seqblock <- function(object = NULL, id.vars, id.vals, exact.vars = NULL, exact.vals = NULL, exact.restr = NULL, exact.alg = "single", covar.vars = NULL, covar.vals = NULL, covar.restr = NULL, covars.ord = NULL, n.tr = 2, tr.names = NULL, assg.prob = NULL, seed = NULL, seed.dist, assg.prob.stat = NULL, trim = NULL, assg.prob.method = NULL, assg.prob.kfac = NULL, distance = NULL, file.name = NULL, query = FALSE, verbose = TRUE, ...){
  
  if(is.null(object)){
    if(query==TRUE){
      init <- readline("Is this the first unit that is being assigned in this experiment?  [y/n]  ")     
      if(!(substr(init,1,1) %in% c("y", "n"))){    
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
          init <- "y"
        }else(stop("The seqblock function requires the ''object'' argument to be set to a valid file name for subsequent unit assignment."))
      }
      if(substr(init,1,1) == "n"){
        object <- readline("Enter the name of the input file without quotation marks.  [E.g., sbout1.RData]  ")
      }
    }
  }
  
  if(is.null(object)){  ## If assigning first unit:
    if(query == TRUE){
      ## ID
      lid <- readline("How many identification variables are there?  ")
      id.vars <- id.vals <- NULL
      for(ii in 1:lid){ # For each identification variables,
        ## tmp stores the name of ID variable:
        tmp <- readline(cat("Enter the name of ID variable ", ii, " without quotation marks.", sep=""))
        ## Add the ii-th name of ID variable to a vector id.vars:
        id.vars <- append(id.vars, tmp)
        ## tmp2 stores the value of each ID variable:
        tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
        ## Add the ii-th value of ID variable to a vector id.vals
        id.vals <- append(id.vals, tmp2)
      } ## End loop over ID variables
      
      nnn <- readline("Should the ID values be numeric? [y/n]  ")
      if(!(substr(nnn,1,1) %in% c("y", "n"))){ 
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline("Should the ID values be numeric? [y/n]  ")
        }
      }
      if(substr(nnn,1,1) != "n"){    # If the respondent answers that the ID values should be numeric,
        id.vals <- as.numeric(id.vals) # change id.vals to numerics.
      }
      if(sum(is.na(id.vals)) > 0){ # If there is at least one missing value in id.vals, 
        tmp <- readline(cat("Warning: ID value(s)", which(is.na(id.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
        if(!(substr(tmp,1,1) %in% c("y", "n"))){
          ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
          if(substr(ddd, 1, 1) != "n"){
          }else{
            tmp <- readline(cat("Warning: ID value(s)", which(is.na(id.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))  			   
          }
        }
        if(substr(tmp,1,1) == "n"){
          stop()
        }
      }
      
      ## Specifying EXACT blocking variables:
      lex <- readline("`Exact blocking variables' are those on which you require covariate values to be identical.  How many exact blocking variables are there?  ")
      exact.vars <- exact.vals <- exact.restr <- NULL
      if(lex > 0){ ## If there are exact blocking variables
        for(ii in 1:lex){ ## Then, for each exact blocking variable
          ## tmp stores the name of exact blocking variable:
          tmp <- readline(cat("Enter the name of exact block variable ", ii, " without quotation marks.", sep="")) 
          exact.vars <- append(exact.vars, tmp) ## then add to the vector exact.vars
          tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
          tmp3 <- readline(cat("Should '", tmp, "' be restricted to certain values? [n/y]  ", sep=""))
          if(!(substr(tmp3,1,1) %in% c("y", "n"))){        ## if the answer is neither yes nor no, 
            ddd <- readline("The default is 'no'.  Continue? [y/n]  ") 
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp3 <- readline(cat("Should '", tmp, "' be restricted to certain values? [n/y]  "))
            }
          }  ## Close 'if'
          
          if(substr(tmp3, 1, 1) == "y"){ ## if the ii-th exact variable should be restricted to certain values, tmp4 stores the values that variable can take that the respondent inputs.
            tmp4 <- readline(cat("How many values can '", tmp, "' take?  ", sep=""))
            ## if the value is less than one, execution is halted.
            if(tmp4 < 1){
              stop("Restricted variables must have at least one valid value.  ")	
            }
            exact.restr <- as.list(exact.restr)
            exact.restr[[tmp]] <- NA
            for(rrr in 1:tmp4){ ## for each value that the variable tmp can take,
              tmp5 <- readline(cat("Enter possible value for '", tmp, "' number ", rrr, sep=""))
              exact.restr[[tmp]] <- append(exact.restr[[tmp]], tmp5) ## the object named tmp in the list exact.restr will be a vector of the values that the variable can take. 					
            } ## End loop over values 'tmp' can take
          } ## Close 'if'
          ## if the ii-th exact variable should NOT be restricted to certain values, exact.vals simply contains the values of exact blocking variables. 
          exact.vals <- append(exact.vals, tmp2)			
        } ## Close loop over exact blocking variables 'ii in 1:lex'
        nnn <- readline("Should the exact block values be numeric? [y/n]  ")
        if(!(substr(nnn, 1, 1) %in% c("y", "n"))){
          ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
          if(substr(ddd, 1, 1) != "n"){
          }else{
            nnn <- readline("Should the exact block values be numeric? [y/n]  ")						
          }
        }
        ## If the answer is yes (the exact block values should be numeric), 
        if(substr(nnn,1,1) != "n"){
          ## change values in exact.vals into numerics:
          exact.vals <- as.numeric(exact.vals)
          ## If there is a restriction on the values of exact variable, 
          if(length(exact.restr) > 0){
            for(rrr in 1:length(exact.restr)){
              ## Each element of those values are changed to numerics.
              exact.restr[[rrr]] <- as.numeric(exact.restr[[rrr]])
            }
          }
        }
        
        if(sum(is.na(exact.vals)) > 0){
          tmp <- readline(cat("Warning: Exact blocking value(s)", which(is.na(exact.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
          if(!(substr(tmp, 1, 1) %in% c("y", "n"))){
            ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp <- readline(cat("Warning: Exact blocking value(s)", which(is.na(exact.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
            }
          }
          if(substr(tmp,1,1) == "n"){
            stop()
          }
        } ## End confirm NA values for exact blocking vars
        ## CHECK that each exact.val is in exact.restr here.
      } ## End 'if(lex > 0)'
      
      ## COVARIATES
      lcv <- readline("How many other blocking variables are there?  ")
      covar.vars <- covar.vals <- covar.restr <- NULL
      if(lcv > 0){ ## If there is at least one non-exact blocking variable,
        for(ii in 1:lcv){  ## For each of those variables,
          tmp <- readline(cat("Enter the name of blocking variable ", ii, " without quotation marks.", sep="")) 
          ## Save its name to a vector 'covar.vars'
          covar.vars <- append(covar.vars, tmp) 
          tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
          tmp3 <- readline(cat("Should '", tmp, "' be restricted to certain values? [n/y]  ", sep=""))
          if(!(substr(tmp3, 1, 1) %in% c("y", "n"))){
            ddd <- readline("The default is 'no'.  Continue? [y/n]  ")
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp3 <- readline(cat("Should '", tmp, "' be restricted to certain values? [n/y]  ", sep=""))
            }
          } ## close if
          ## If the value of the blocking variable should be restricted to certain values,
          if(substr(tmp3,1,1) == "y"){
            tmp4 <- readline(cat("How many values can '", tmp, "' take?  ", sep=""))
            if(tmp4 < 1){
              stop("Restricted variables must have at least one valid value.  ")	
            }
            covar.restr <- as.list(covar.restr)
            covar.restr[[tmp]] <- NA
            for(rrr in 1:tmp4){
              tmp5 <- readline(cat("Enter possible value for '", tmp, "' number ", rrr, sep=""))
              # Make a list covar.restr in which each object represents restricted values of blocking variable that should be restricted to certain values. 
              covar.restr[[tmp]] <- append(covar.restr[[tmp]], tmp5) 						
            }
          } 
          ## If the values should NOT be restricted to certain values, covar.vals contains the values of the blocking variable:
          covar.vals <- append(covar.vals, tmp2)
        } ## End loop over blocking vars 'ii in 1:lcv' 
        
        nnn <- readline("Should the blocking values be numeric? [y/n]  ")
        
        if(!(substr(nnn,1,1) %in% c("y", "n"))){
          ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
          if(substr(ddd, 1, 1) != "n"){
          }else{
            nnn <- readline("Should the blocking values be numeric? [y/n]  ")						
          }
        }
        
        if(substr(nnn,1,1) != "n"){
          covar.vals <- as.numeric(covar.vals)
        }else{
          stop("The blocking variables are stored as character strings.\nCurrent implementation does not support distance calculations for character strings.")
        }  ## (Change to warning() when categorical distances incorporated)
        
        if(sum(is.na(covar.vals)) > 0){
          tmp <- readline(cat("Warning: Blocking value(s)", which(is.na(covar.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
          
          if(!(substr(tmp,1,1) %in% c("y", "n"))){
            ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp <- readline(cat("Warning: Blocking value(s)", which(is.na(covar.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
            }
          }						
          if(substr(tmp,1,1) == "n"){
            stop()
          }
        }
        
        ## Covariates' order
        nnn <- readline("Would you like to specify the order in which blocking variables will be used while (number units considered) < (number covariates + 2) ?  (If not, covariates will be added as possible in the order in which they were entered.)  [n/y]  ")
        if(!(substr(nnn,1,1) %in% c("y", "n"))){
          ddd <- readline("The default is `no'.  Continue? [y/n]  ")
          if(substr(ddd, 1, 1) != "n"){
          }else{
            nnn <- readline(cat("Would you like to specify the order in which blocking variables will be used while (number units considered) < (number covariates +2) ? [n/y]  "))
          }
        }							
        if(substr(nnn,1,1) == "y"){
          for(ii in 1:lcv){
            tmp <- readline(cat("Enter the name of blocking variable prioritized number ", ii, ".", sep=""))
            #					if(!(tmp %in% names(x)[(lid+lex+1):(lid+lex+lcv)])){
            if(!(tmp %in% covar.vars)){
              tmp <- readline(cat("`", tmp, "' is not a blocking variable.  Please try again.  Enter the name of blocking variable prioritized number ", ii, ".", sep=""))
            }
            if(tmp %in% covars.ord){
              tmp <- readline(cat("`", tmp, "' has already been prioritized.  Please try again.  Enter the name of blocking variable prioritized number ", ii, ".", sep=""))
            }
            covars.ord <- append(covars.ord, tmp)
          }
        }
      }  ## End 'if(lcv > 0)'
      
      ## Number of treatment conditions:
      n.tr <- readline("How many experimental/treatment conditions are there?  ")
      n.tr <- as.numeric(n.tr)
      
      ## Names of treatment conditions:
      nnn <- readline(cat("Would you like to specify the names of the experimental conditions? [y/n]  [If not, `Treatment 1', ..., `Treatment ", n.tr, "', will be used]  ", sep=""))
      if(substr(nnn,1,1) != "n"){
        tr.names <- NULL
        for(ii in 1:n.tr){
          tmp <- readline(cat("Enter the name of experimental condition ", ii, " without quotation marks.", sep="")) 
          tr.names <- append(tr.names, tmp)
        }					
      }
      
      ## Initial condition assignment probabilities:
      nnn <- readline("Would you like to specify the initial assignment probabilities?  [y/n]  [If not, all treatment assignment probabilities for the first unit will be equal.]  ")
      
      if(!(substr(nnn,1,1) %in% c("y", "n"))){
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline("Would you like to specify the initial assignment probabilities?  [y/n]  [If not, all treatment assignment probabilities for the first unit will be equal.]  ")
        }
      }
      
      if(substr(nnn,1,1) != "n"){
        assg.prob <- NULL
        for(ii in 1:n.tr){
          tmp <- readline(cat("Enter the assignment probability of experimental condition ", ii, ".", sep="")) 
          assg.prob <- append(assg.prob, tmp)
        }
        assg.prob <- as.numeric(assg.prob)
      }		
      
      ## Set random seed:
      nnn <- readline("Would you like to specify a random seed for the initial experimental condition assignment?  [y/n]  ")
      
      if(!(substr(nnn,1,1) %in% c("y", "n"))){
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline("Would you like to specify a random seed for the initial experimental condition assignment?  [y/n]  ")
        }
      }
      
      if(substr(nnn,1,1) != "n"){
        seed <- as.numeric(readline("Enter the random seed.  "))
      }
      
      ## ASSIGNMENT PROBABILITY STATISTIC:
      nnn <- readline("Would you like to specify the assignment probability summary statistic?  [n/y]  [If not, `mean' will be used.]  ")
      if(substr(nnn,1,1) == "y"){
        assg.prob.stat <- readline("Enter the summary statistic name.  [mean, median, trimmean]  ")
        if(assg.prob.stat == "trimmean"){
          ooo <- readline("Specify the proportion of blocks to be dropped from each end of the distribution before the mean is calculated.  [Defaults to 0.1]  ")
          if(substr(ooo,1,1) == ""){
            trim <- 0.1
          }else{
            trim <- ooo
          }
        }
      }else{
        assg.prob.stat <- "mean"
      }
      
      ## ASSIGNMENT PROBABILITY METHOD:
      nnn <- readline("Would you like to specify the assignment probability algorithm?  [n/y]  [If not, `ktimes' with k=2 will be used.]  ")
      if(substr(nnn,1,1) == "y"){
        assg.prob.method <- readline("Enter the algorithm name.  [ktimes, fixed, prop, prop2, wprop]  ")
        if(assg.prob.method == "ktimes"){
          assg.prob.kfac <- as.numeric(readline("Enter the value of k, the factor by which the most likely experimental condition will be multiplied, relative to the other conditions.  "))
          }      
        }else{
        assg.prob.method <- "ktimes"
      }
      
      if(assg.prob.method == "fixed"){
        assgpr <- readline("Enter 1st probability.  ")
        for(pc in 2:n.tr){
          assg.prob <- append(assg.prob, readline(cat("Enter probablity number ", pc, ".  ", sep="")))
        }
        assg.prob <- as.numeric(assg.prob)
      }
      
      ## MULTIVARIATE DISTANCE CALCULATION:
      nnn <- readline("Would you like to specify how the multivariate distance used for blocking is calculated?  [n/y]  [If not, `mahalanobis' will be used.]  ")
      if(substr(nnn,1,1) == "y"){
        distance <- readline("Enter the calculation type.  [mahalanobis, mcd, mve, euclidean]  ")
      }else{
        distance <- "mahalanobis"
      }
      
      ## OUTPUT FILE NAME:
      nnn <- readline("Would you like to specify the output data file name?  [y/n]  [If not, `sbout1.RData' will be used.]  ")
      if(!(substr(nnn,1,1) %in% c("y", "n"))){
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline("Would you like to specify the output data file name?  [y/n]  [If not, `sbout1.RData' will be used.]  ")
        }
      }		
      
      if(substr(nnn,1,1) != "n"){
        file.name <- readline("Enter the file name.  Note: seqblock2k() expects a file of type `.RData'.  ")
      }						
    } ## END (query == T)
    
    ## create variables
    lid <- length(id.vars)
    lex <- length(exact.vars)
    lcv <- length(covar.vars)
    
    ## create vector of treatment names
    if(is.null(tr.names)){
      tr.names <- paste("Treatment", 1:n.tr)
    }
    
    ## create kfac, if null
    if(!(is.null(assg.prob.kfac))){
      kfac <- assg.prob.kfac
    } else {kfac <- assg.prob.kfac <- 2}
       
    ## create assg.prob.stat, if null
    if(!(is.null(assg.prob.stat))){
      assg.prob.stat <- assg.prob.stat
    } else {assg.prob.stat <- "mean"}
    
    ## create trim, if null
    if(!(is.null(trim))){
      trim <- trim
    } else {trim <- 0.1}
    
    ## create assg.prob.method, if null
    if(!(is.null(assg.prob.method))){
      assg.prob.method <- assg.prob.method
    } else {assg.prob.method <- "ktimes"}
    
    ## create distance, if null
    if(is.null(distance)){
      distance <- "mahalanobis"
    }
    
    ## create vector of (equal) initial treatment probabilities
    if(is.null(assg.prob)){
      assg.prob <- rep(1/n.tr, n.tr)
    }
    
    ## check that assg.prob sum to 1  
    if(sum(assg.prob)!=1){
      stop("Initial assignment probabilities do not sum to 1.\nRespecify and try again.")
    } 
    
    ## check that if exact named, then supplied
    if((!is.null(exact.vars)) && is.null(exact.vals)){
      stop("Exact (grouped) blocking requested, but no values given.\nSpecify exact block values and try again.")
    } 
    
    ## check that exact vars and vals are same length
    if(lex != length(exact.vals)){
      stop(cat("Number of exact blocking variables (", lex, ") does not equal number of values (", length(exact.vals), ").\nSpecify exact block values and try again.\n", sep=""))	
    } 
    
    ## check that if covars named, then supplied
    if((!is.null(covar.vars)) && is.null(covar.vals)){
      stop("Covariate blocking requested, but no values given.\nSpecify covariate values and try again.")
    } 
    
    ## check that covar vars and vals are same length
    if(length(covar.vars) != length(covar.vals)){
      stop("Number of covariate blocking variables does not equal number of values.\nSpecify covariate block values and try again.")
    }
    
    ## Add NA to each existing restriction set, 
    ##   then check if any exact vals violate exact.restr
    if((lex > 0) && !(is.null(exact.restr))){
      for(i in names(exact.restr)){
        if(!(NA %in% exact.restr[[i]])){
          exact.restr[[i]] <- append(exact.restr[[i]], NA)
        }
        wexrstr <- which(exact.vars == i)
        if(!(exact.vals[wexrstr] %in% exact.restr[[exact.vars[wexrstr]]])){
          stop("Exact blocking value number ", wexrstr, " (", exact.vals[wexrstr], ") is not in the set of values to which exact blocking variable `", i, "' is restricted.  Respecify either the value, the restriction set, or both as needed." )
        }
      }
    }   
    
    ## check if any covar vals violate covar.restr
    if((lcv > 0) && (!is.null(covar.restr))){
      for(i in names(covar.restr)){
        if(!(NA %in% covar.restr[[i]])){
          covar.restr[[i]] <- append(covar.restr[[i]], NA)
        }
        wcvrstr <- which(covar.vars == i)
        if(!(covar.vals[wcvrstr] %in% covar.restr[[covar.vars[wcvrstr]]])){
          stop("Covariate blocking value number ", wcvrstr, " (", covar.vals[wcvrstr], ") is not in the set of values to which covariate blocking variable `", i, "' is restricted.  Respecify either the value, the restriction set, or both as needed." )
        }
      }	
    }
    
    ## set user-assigned seed
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    ## Draw initial treatment assignment
    init.tr <- sample(tr.names, 1, replace = FALSE, prob = assg.prob)
    
    ## Create data frame of unit 1 data
    ## create id placeholders
    id.tmp <- rep(NA, lid)
    names(id.tmp) <- id.vars
    ## create group placeholders
    ex.tmp <- rep(NA, lex)
    names(ex.tmp) <- exact.vars
    
    ## create data frame
    x <- unlist(c(id.tmp, ex.tmp, unname(covar.vals), NA))
    x <- data.frame(t(x))
    colnames(x) <- c(id.vars, exact.vars, covar.vars, "Tr")
    ## add id vals.  NOTE: will be all char or all num
    x[1, 1:lid] <- id.vals
    ## add group vals.  NOTE: will be all char or all num
    if(!is.null(exact.vars)){
      x[1, (lid+1):(lid+lex)] <- exact.vals
    }
    ## add initial treatment assignment
    x[1,ncol(x)] <- init.tr
    
    if(is.null(file.name)){
      file.name <- "sbout1.RData"
    }
    
    if(substr(file.name, start=nchar(file.name)-5, stop=nchar(file.name)) !=".RData"){
      warning("Destination file name does not end in .RData.")
    } 
    
    ## Create storage list 'bdata'
    bdata <- list()
    bdata$x <- x
    bdata$nid <- id.vars
    bdata$nex <- exact.vars
    bdata$ncv <- covar.vars
    bdata$rex <- exact.restr  ## EXACT restricted values list
    if(!is.null(bdata$rex)){
      names(bdata$rex) <- exact.vars
    }
    bdata$rcv <- covar.restr  ## BLOCK restricted values list
    if(!is.null(bdata$rcv)){
      names(bdata$rcv) <- covar.vars
    }
    bdata$ocv <- covars.ord   ## block covariates order
    bdata$trn <- tr.names
    bdata$apstat <- assg.prob.stat 
    bdata$mtrim <- trim
    bdata$apmeth <- assg.prob.method
    bdata$kfac <- kfac ## multiple for method 'ktimes'
    bdata$assgpr <- assg.prob
    bdata$distance <- distance
    bdata$datetime <- date()
    bdata$orig <- x
    
    ## Write bdata to file
    save(bdata, file = file.name)
    if(verbose == TRUE){
      cat("Unit 1 data stored as file ", file.name, ".\nThe current working directory is ", getwd(), "\n", sep="")
      cat("Unit ", nrow(x), " assigned to ", init.tr, ".\n", sep="")
      cat("The new data as entered:\n")
      print(x)
    }	
  } ## End 'if(is.null(object))'.  I.e., end assigning first unit.
  
  ## Next case: if object is not NULL (I.e., begin assignment for subsequent unit after the first one)
  if(!is.null(object)){ 
    
    load(object)  ## loads, doesn't export to wkspace
    
    ## Rename object elements	
    x <- bdata$x  
    nid <- bdata$nid  ## names of id.vars
    nex <- bdata$nex  ## names of exact.vars
    ncv <- bdata$ncv  ## names of covar.vars
    rex <- bdata$rex  ## EXACT restricted values list
    rcv <- bdata$rcv  ## BLOCK restricted values list
    ocv <- bdata$ocv  ## block covariates order
    trn <- bdata$trn  ## treatment condition names
    apstat <- bdata$apstat ## assignment prob statistic
    mtrim <- bdata$mtrim ## trim fraction for apmeth=trimmean
    apmeth <- bdata$apmeth ## assignment prob method
    kfac <- bdata$kfac     ## assg multiple for method 'ktimes'
    assgpr <- bdata$assgpr ## assg probs for fixed prob
    dist <- bdata$distance ## distance calculation method
    orig <- bdata$orig  ## original unjittered data
    n.tr <- length(trn) ## number of treatment conditions

    lid <- length(nid) #number of id variables
    lex <- length(nex) #number of exact variables
    lcv <- length(ncv) #number of covariates  
    
    if(!(is.null(exact.restr))){
      rex <- exact.restr	
    }
    if(!(is.null(covar.restr))){
      rcv <- covar.restr	
    }
    if(!(is.null(covars.ord))){
      ocv <- covars.ord	
    }
    if(!(is.null(assg.prob.stat))){
      if(all.equal(assg.prob.stat,apstat) != TRUE){
      warning(paste("'assg.prob.stat' of previous assignments was ", apstat, ". 'assg.prob.stat' for this assignment is ", assg.prob.stat, ".", sep = ""))
      }
      apstat <- assg.prob.stat
    }
    if(!(is.null(trim))){
      if(all.equal(trim,mtrim) != TRUE){
      warning(paste("'trim' of previous assignments was ", mtrim, ". 'trim' for this assignment is ", trim, ".", sep = ""))
      }
      mtrim <- trim	
    }
    if(!(is.null(assg.prob.method))){
      if(all.equal(assg.prob.method, apmeth) != TRUE){
      warning(paste("'assg.prob.method' of previous assignments was ", apmeth, ". 'assg.prob.method' for this assignment is ", assg.prob.method, ".", sep = ""))
      }
      apmeth <- assg.prob.method
    }	
    if(!(is.null((assg.prob.kfac)))){
      if(all.equal(assg.prob.kfac,kfac) != TRUE){
      warning(paste("'kfac' of previous assignments was ", kfac, ". 'kfac' for this assignment is ", assg.prob.kfac, ".", sep = ""))
      }
      kfac <- assg.prob.kfac
    }		
    if(!(is.null(assg.prob))){
      assgpr <- assg.prob
    }
    if(!(is.null((distance)))){
      if(all.equal(distance,dist) != TRUE){
        warning(paste("Distance option of previous assignments was '", dist, "'. Distance option for this assignment is '", distance, "'.", sep = ""))
      }
      dist <- distance
    }  
      
    if(query == TRUE){
      ## ID
      id.vals <- NULL
      for(ii in 1:lid){
        tmp <- names(x)[ii]
        tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
        id.vals <- append(id.vals, tmp2)
      }
      class(id.vals) <- class(x[,lid])  ## lid could be any of 1:lid
      
      if(sum(is.na(id.vals)) > 0){ # If there is at least one missing value in id.vals, 
        tmp <- readline(cat("Warning: ID value(s)", which(is.na(id.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))		
        if(!(substr(tmp, 1, 1) %in% c("y", "n"))){
          ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
          if(substr(ddd, 1, 1) != "n"){
          }else{
            tmp <- readline(cat("Warning: ID value(s)", which(is.na(id.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
          }
        }							
        if(substr(tmp, 1, 1) == "n"){
          stop()
        }
      }
      
      ## Check for duplicated IDs if single id.var
      if((lid == 1) && (id.vals %in% x[, nid])){
        tmp <- readline(cat("Warning: ID value ", id.vals, " already used in the data.  Allow duplication and proceed? [y/n]", sep = ""))
        if(substr(tmp, 1, 1) == "n"){
          stop()	
        }		
      }
      
      ## EXACT
      exact.vals <- NULL
      if(lex > 0){
        for(ii in 1:lex){
          tmp <- nex[ii]
          tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
          exact.vals <- append(exact.vals, tmp2)
        }
        class(exact.vals) <- class(x[,(lid+lex)]) ## could be any of (lid+1):(lid+lex) If the exact blocking variable was set to be numeric in seqblock1(), new exact values should be also numeric. Otherwise, the values of the exact blocking variable will be coerced to NAs.		
        if(sum(is.na(exact.vals)) > 0){
          tmp <- readline(cat("Warning: Exact blocking value(s)", which(is.na(exact.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
          
          if(!(substr(tmp,1,1) %in% c("y", "n"))){
            ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp <- readline(cat("Warning: Exact blocking value(s)", which(is.na(exact.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
            }
          }							
          if(substr(tmp, 1, 1) == "n"){
            stop()
          }
        }
      }
      
      ## EXACT ALGORITHM
      ##   "single": if lex > 1, creates unique categories (4 for 2x2, e.g.)
      exact.alg <- "single"
      
      ## COVARIATES
      covar.vals <- NULL
      if(lcv > 0){
        for(ii in 1:lcv){
          tmp <- ncv[ii]
          tmp2 <- readline(cat("Enter the value of '", tmp, "'.  ", sep="")) 
          covar.vals <- append(covar.vals, tmp2)
        }
        class(covar.vals) <- class(x[,lid+lex+lcv]) ## could be any of (lid+lex+1):(lid+lex+lcv)
        
        if(sum(is.na(covar.vals)) > 0){
          tmp <- readline(cat("Warning: Blocking value(s)", which(is.na(covar.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
          if(!(substr(tmp,1,1) %in% c("y", "n"))){
            ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
            if(substr(ddd, 1, 1) != "n"){
            }else{
              tmp <- readline(cat("Warning: Blocking value(s)", which(is.na(covar.vals)), "was/were coerced to and stored as `NA'.  Proceed? [y/n] ", sep=" "))
            }
          }				
          if(substr(tmp,1,1) == "n"){
            stop()
          }
        }
      }
      
      ## SET RANDOM SEED:
      nnn <- readline("Would you like to specify a random seed for this experimental condition assignment?  [y/n]  ")
      
      if(!(substr(nnn,1,1) %in% c("y", "n"))){
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline(cat("Would you like to specify a random seed for this experimental condition assignment?  [y/n]  "))
        }
      }						
      if(substr(nnn,1,1) != "n"){
        seed <- as.numeric(readline("Enter the random seed.  "))
      }
      
      ## OUTPUT FILE NAME:
      nnn <- readline("Would you like to specify the output data file name?  [y/n]  [If not, `sbout2k.RData' will be used.]  ")
      if(!(substr(nnn,1,1) %in% c("y", "n"))){
        ddd <- readline("The default is `yes'.  Continue? [y/n]  ")
        if(substr(ddd, 1, 1) != "n"){
        }else{
          nnn <- readline(cat("Would you like to specify the output data file name?  [y/n]  [If not, `sbout2k.RData' will be used.]  "))
        }
      }						
      if(substr(nnn,1,1) != "n"){
        file.name <- readline("Enter the file name.  Note: seqblock2k() expects a file of type `.RData'.  ")
      }else{
        file.name <- NULL
      }
    } ## Close query==T 
    
    ## Warn if duplicated ID values if single ID var
    if((lid == 1) && (id.vals %in% x[, nid])){
      warning("Warning: ID value ", id.vals, " already used in the data.")
    }
    
    ## check # exact vars vs. # exact vals
    if(length(nex) != length(exact.vals)){
      stop(cat("Number of exact blocking variables (", lex, ") does not equal number of values (", length(exact.vals), ").\nSpecify exact block values and try again.\n", sep=""))
    } 
    
    ## check if any exact val violates exact.restr
    if((lex > 0) && !(is.null(rex))){
      for(i in names(rex)){
        if(!(NA %in% rex[[i]])){
          rex[[i]] <- append(rex[[i]], NA)
        }
        wexrstr <- which(nex == i)
        if(!(exact.vals[wexrstr] %in% rex[[nex[wexrstr]]])){
          tmp <- readline(cat("Exact blocking value number ", wexrstr, " (", exact.vals[wexrstr], ") is not in the set of values to which exact blocking variable `", i, "' is restricted.  Would you like to add ", exact.vals[wexrstr], " to the set?  [y/n]", sep=""))
          if(substr(tmp,1,1) != "n"){
            rex[[i]] <- append(rex[[i]], exact.vals[wexrstr])
          }else{
            stop("Exact blocking value number ", wexrstr, " (", exact.vals[wexrstr], ") is not in the set of values to which exact blocking variable `", nex[wexrstr], "' is restricted.  Respecify either the value, the restriction set, or both as needed." )
          }
        }
      }
    }
    
    ## check if any block val violates covar.restr  
    if((lcv > 0) && (!is.null(rcv))){		for(i in names(rcv)){
      if(!(NA %in% rcv[[i]])){
        rcv[[i]] <- append(rcv[[i]], NA)
      }
      wcvrstr <- which(ncv == i)
      if(!(covar.vals[wcvrstr] %in% rcv[[ncv[wcvrstr]]])){
        tmp <- readline(cat("Covariate blocking value number ", wcvrstr, " (", covar.vals[wcvrstr], ") is not in the set of values to which covariate blocking variable `", i, "' is restricted.  Would you like to add ", covar.vals[wcvrstr], " to the set?  [y/n]", sep=""))
        if(substr(tmp,1,1) != "n"){
          rcv[[i]] <- append(rcv[[i]], covar.vals[wcvrstr])
        }else{
          stop("Covariate blocking value number ", wcvrstr, " (", covar.vals[wcvrstr], ") is not in the set of values to which covariate blocking variable `", ncv[wcvrstr], "' is restricted.  Respecify either the value, the restriction set, or both as needed.")
        }
      }
    }
    }	
    
    if(!(length(ocv) == length(unique(ocv)))){
      stop("Prioritization of covariates for (number units considered) < (number covariates + 2) includes duplicate variable names.  Respecify 'covars.ord'.")	
    }	
    
    if(lex > 0){
      if(exact.alg == "single"){
        x.ex <- x[apply(t(t(x[,(lid+1):(lid+lex)])==exact.vals), 1, sum) == lex,]
      }
    }  
    
    ## When no exact matches, use all previous data:
    if(is.null(exact.vals) || (nrow(x.ex)==0)){
      x.ex <- x
    } 
    
    prev <- nrow(x.ex)
    
    ## Create data frame of unit 1 data
    ## create id placeholders
    id.tmp <- rep(NA, lid)
    names(id.tmp) <- nid
    ## create exact block placeholders
    ex.tmp <- rep(NA, lex)
    names(ex.tmp) <- nex
    
    qqq <- rbind(x.ex[,1:(ncol(x.ex)-1)], c(id.tmp, ex.tmp, covar.vals))
    qqq[nrow(qqq), 1:lid] <- id.vals
    if(!is.null(exact.vals)){
      qqq[nrow(qqq), (lid+1):(lid+lex)] <- exact.vals
    }
    
    orig <- rbind(orig, c(qqq[nrow(qqq), ], Tr=NA))
    
    ## When non-exact blocking covariates are used:
    ## (Default: covar order = order in dataframe)
    if(lcv > 0){
      if(is.null(ocv)){  
        ocv <- ncv  
      }
      
      ## Select covars: if prev = 1, pick 1.  
      if(prev == 1){
        cov.tmp <- ocv[1]  
      }
      ## if prev>=2, pick prev-1 covars.
      ##   to define MD, need num covars, lcv >= prev. 
      ##   for distinct MD, need num covars, lcv >= prev+2
      if(prev >= 2){
        if(prev-2 >= lcv){
          cov.tmp <- ocv  ##changed to ocv from covars.ord
        }else{
          cov.tmp <- ocv[1:(prev-1)]
        } 
      }
      qqq.c <- data.frame(qqq[ , cov.tmp]) 
      names(qqq.c) <- cov.tmp
      
      ## If there is a variable with no variation, remove it from the vcov matrix for this assignment:
      wh.cut <- NULL
      for(jj in 1:(ncol(qqq.c))){
        if(isTRUE(var(qqq.c[, jj]) == 0)){
          wh.cut <- append(wh.cut, jj)
        }
      }
      if(length(wh.cut) > 0){
        n.qqq.c <- names(qqq.c)[-wh.cut]
        qqq.c <- data.frame(qqq.c[, -wh.cut])
        names(qqq.c) <- n.qqq.c
      }
      
      ## If qqq.c has zero columns, then, go recreate qqq.c.
      ## (This occurs if second unit is identical to first on first nonexact blocking variable, e.g.)
      if(ncol(qqq.c) == 0){
        qqq.c <- data.frame(qqq[ , cov.tmp]) 
        names(qqq.c) <- cov.tmp
      }
        
      ## If new unit identical to old unit, jitter new unit's covar vals
      for(jj in 1:(nrow(qqq.c)-1)){
        if(isTRUE(sum((qqq.c[jj,] - qqq.c[nrow(qqq.c),])==0) == length(qqq.c[1,]))){ 
          jit.val <- array()
          for(kk in 1:ncol(qqq.c)){
            jit.val[kk] <- jitter(qqq.c[nrow(qqq.c), kk])
          }
          qqq[nrow(qqq), cov.tmp] <- qqq.c[nrow(qqq.c), ] <- jit.val 
        }
      }  ## To implement: jitter for identical rows w/ identical NAs.
      
      ## replace NA's with variable means
      for(jj in 1:ncol(qqq.c)){
        if(sum(is.na(qqq.c[, jj])) > 0 ){
          qqq.c[which(is.na(qqq.c[, jj])), jj] <- mean(qqq.c[, jj], na.rm = TRUE)
        }
      }
      
      ## Different distance measures
      if(is.character(dist)){
        ## Since mve and mcd require conditions on IQR, check, warn, use nonresistant:
        if(dist %in% c("mve", "mcd")){
          quan <- floor((nrow(qqq.c)+ncol(qqq.c)+1)/2)
          if(quan < (ncol(qqq.c)+1)){
            warning(paste("'Quantile' must be at least ", ncol(qqq.c)+1, " when using option '", dist, "'. Blocking will proceed using nonresistant Mahalanobis distance scaling matrix.", sep = ""))
            dist <- "mahalanobis"
          } else if(quan > (nrow(qqq.c)-1)){
            warning(paste("'Quantile' must be at most ", nrow(qqq.c)-1, " when using option '", dist, "'. Blocking will proceed using nonresistant Mahalanobis distance scaling matrix.", sep = ""))
            dist <- "mahalanobis"
          }
        }
        if(dist %in% c("mve", "mcd")){
          iqr.idx <- 0
          while(iqr.idx < ncol(qqq.c)){
            iqr.idx <- iqr.idx + 1
            iqr.tmp <- unname(quantile(qqq.c[, iqr.idx], c(.25, .75)))
            if(isTRUE(all.equal(iqr.tmp[1], iqr.tmp[2]))){
              warning(paste("Variable ", colnames(qqq.c)[iqr.idx], " has IQR 0; blocking will proceed using nonresistant Mahalanobis distance scaling matrix.", sep = ""))
              dist <- "mahalanobis"
              iqr.idx <- ncol(qqq.c)
            }
          }
        }
        if(dist == "mahalanobis"){
          vc.all <- var(qqq.c)
        } else if(dist == "mcd"){
          vc.all <- cov.rob(qqq.c, method="mcd", seed = seed.dist, ...)$cov
        } else if(dist == "mve"){
          vc.all <- cov.rob(qqq.c, method="mve", seed = seed.dist, ...)$cov
        } else if(dist == "euclidean"){
          vc.all <- diag(ncol(qqq.c))
        } else{
          dist <- "mahalanobis"
          vc.all <- var(qqq.c)
        }
      }  ## End 'if(is.character(dist))'
      
      mahmat <- mahal(qqq.c, vc.all)	
      
      tr.dist <- list()
      lrow <- mahmat[nrow(mahmat), 1:(ncol(mahmat)-1)]
      ##	for(ii in levels(data1[, ncol(data1)])){}  ## good for factor TR
      for(ii in unique(x.ex[, ncol(x.ex)])){  ## good for character TR
        tr.dist[[ii]] <- lrow[which(x.ex[, ncol(x.ex)]==ii)]
      }      
      
      ## Calculate summary for comparison
      if(apstat == "mean"){
        ms <- lapply(tr.dist, mean)  ## only prev assg'd tr's
      }
      if(apstat == "median"){
        ms <- lapply(tr.dist, median)  
      }	
      if(apstat == "trimmean"){
        ms <- lapply(tr.dist, mean, trim = mtrim)
      }
      
      ## Sort previously assigned treatments, largest first
      tr.sort <- names(sort(unlist(ms), decreasing = TRUE))  
      ## Add treatments not yet assigned, in random order
      if(length(tr.sort) != n.tr){
        trn.unassg <- sample(trn[!(trn %in% names(ms))])
        nunassg <- length(trn.unassg)
        tr.sort <- c(trn.unassg, tr.sort)
      }
      
      if((length(ms) != length(tr.sort)) && (apmeth %in% c("prop", "prop2", "wprop"))){
        ms <- append(ms, rep(NA, nunassg), after=0)
        names(ms)[1:nunassg] <- trn.unassg
        ## Give unassigned tr's double distance of max condition dist
        ## (Quasi-minimization)
        ms[1:nunassg] <- 2*max(unlist(ms), na.rm=T)
      }
      
      if((length(tr.dist) != length(tr.sort)) && (apmeth == "wprop")){
        tr.dist <- append(tr.dist, rep(NA, nunassg), after=0)
        names(tr.dist)[1:nunassg] <- trn.unassg
        ## Give unassigned tr's double distance of max condition dist, 
        ##   with weight of one unit.
        ## (Quasi-minimization)
        tr.dist[1:nunassg] <- max(unlist(ms), na.rm=T) 
      }
    } ## End 'if(lcv > 0)'
    
    if(lcv == 0){ ## if there are ONLY exact covariates used:
      tr.dist <- list()
      tr.sort <- character()      
      tr.counts <- table(x.ex$Tr)
      trn.unassg <- trn[!(trn %in% c(x.ex$Tr))]
      if(length(trn.unassg) > 0){
        tr.counts <- c(tr.counts, rep(0, length(trn.unassg)))
        names(tr.counts)[(n.tr-length(trn.unassg)+1):n.tr] <- trn.unassg
      }
      p.lcv0 <- (1-tr.counts/sum(tr.counts))/sum(1-tr.counts/sum(tr.counts))
    } ## End 'if(lcv == 0)'
    
    ## Set assignment probability method
    if(apmeth == "ktimes"){
      p <- 	c(kfac/(kfac+n.tr-1), rep(1/(kfac+n.tr-1), n.tr-1))
    }
    
    if(apmeth=="fixed"){
      if(length(assgpr) != n.tr){
        stop("assg.prob.method is `fixed', but assg.prob (vector of assignment probabilities) is length ", length(assgpr), ", while n.tr (number of treatment conditions) is ", n.tr, ".  Respecify and try again.")
      }		
      if(sum(assgpr) != 1){
        stop("assg.prob.method is `fixed', but assg.prob (vector of assignment probabilities) sums to ", sum(assgpr), " instead of 1.  Respecify and try again.")
      }
      p <- assgpr
    }		
    if(apmeth=="prop"){
      mssum <- sapply(ms, sum)
      p <- c(unname((mssum/sum(mssum))[tr.sort]))       
    }			 
    if(apmeth =="prop2"){  
      sumsq <- sapply(ms, sum)^2
      p <- c(unname((sumsq/sum(sumsq))[tr.sort]))
    }
    ## Biases toward tr's w/fewer units
    if(apmeth=="wprop"){
      nnnj <- lapply(tr.dist, length)
      nnn <- sum(unlist(nnnj))
      div.func <- function(x){return(nnn/x)}
      nnnw <- lapply(nnnj, div.func)
      wdist <- unlist(ms)*unlist(nnnw)
      p <- c(unname((wdist/sum(wdist))[tr.sort]))
    }
    
    ## set user-assigned seed
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    ## Assign new TR, given at least 1 "treated as continuous" covars
    if(lcv > 0){
      tr.new <- sample(tr.sort, 1, prob=p)
    }
    ## Assign new TR using correct order for probs and Tr condition names, EXACT covars only
    if(lcv == 0){
      tr.new <- sample(names(tr.counts), 1, prob=p.lcv0)
    }
    ## Assign new TR using same order for fixed probs and Tr condition names
    if(apmeth == "fixed"){
      tr.new <- sample(trn, 1, prob=p)
    }
    
    ## Append to x
    x <- rbind(x, c(qqq[nrow(qqq),], Tr=tr.new))
    orig[nrow(orig),"Tr"] <- tr.new
    tr.counts <- table(x$Tr)
  
    if(is.null(file.name)){
      file.name <- "sbout2k.RData"
    }
    
    if(substr(file.name, start=nchar(file.name)-5, stop=nchar(file.name)) !=".RData"){
      warning("Destination file name does not end in .RData.")
    }
    
    prev.dates <- bdata$datetime
    
    ## Create storage list 
    bdata <- list()
    bdata$x <- x
    bdata$nid <- nid
    bdata$nex <- nex
    bdata$ncv <- ncv
    bdata$rex <- rex  ## EXACT restricted values list
    if(!is.null(bdata$rex)){
      names(bdata$rex) <- nex    
    }
    bdata$rcv <- rcv  ## BLOCK restricted values list
    if(!is.null(bdata$rcv)){
      names(bdata$rcv) <- covar.vars
    }  
    bdata$ocv <- ocv  ## block covariates order
    bdata$trn <- trn
    bdata$apstat <- apstat ## assignment prob statistic
    bdata$mtrim <- mtrim ## trim fraction for apmeth=trimmean
    bdata$apmeth <- apmeth ## assignment prob method
    bdata$kfac <- kfac     ## assg probs for method=ktimes
    bdata$assgpr <- assgpr ## assg probs for fixed prob
    bdata$distance <- dist ## distance calculation method
    bdata$trd <- tr.dist
    bdata$tr.sort <- tr.sort
    if (lcv > 0){
      bdata$p <- p      
    }else{
      bdata$p <- p.lcv0
    }
    bdata$trcount <- tr.counts
    bdata$datetime <- append(prev.dates, date())
    bdata$orig <- orig
    
    save(bdata, file = file.name)
    if(verbose == TRUE){
      cat("Unit 1:", nrow(x), " data stored as file ", file.name, ".\nThe current working directory is ", getwd(), "\n", sep="")
      cat("Unit ", nrow(x), " assigned to ", tr.new, ".\n", sep="")
      cat("The new data as entered:\n")
      print(orig[nrow(orig), ])
    }	
  } ## End 'if(!is.null(object))'
  
  return(bdata)
}