#
# Template for a data mining process (both with or w/o CV)
# (dataset: CPU, method: Random Forest or SVM)
#
# Use "browser()" if you want to look at the variables inside
#
# Author: Wolfgang Konen, FHK, Oct'2009 - Apr'2014
#
# Example usage:
#       result <- main_cpu();
#
main_cpu <- function(opts=NULL,dset=NULL,tset=NULL) {           

    if (is.null(opts)) source("cpu_00.apd", local=TRUE);
    if (opts$rgain.type=="rgain") opts$rgain.type="rmae";
    opts <- tdmOptsDefaultsSet(opts);  # fill in all opts params which are not yet set (see tdmOptsDefaults.r)
    
    gdObj <- tdmGraAndLogInitialize(opts);     # init graphics and log file
        
    #===============================================
    # PART 1: READ DATA
    #===============================================
    if (is.null(dset)) {
      cat1(opts,opts$filename,": Read data ...\n")
      dset <- tdmReadData(opts);
    }

    # which variable is response variable:
    response.variables <- c("ERP")  
    ID.variable <- "ID"
    
    #===============================================
    # PART 2a: DATA PREPARATION
    #===============================================
    # special for cpu-dataset: force columns to be numeric (otherwise R thinks 
    # that they are factors and rowSums below does not work)
    for (n in 2:7) {    
        dset[,n] <- as.numeric(dset[,n])  
    }
    # PREPROC: diminuish the skewness of the target variable  
    #dset[,response.variables] <- log(dset[,response.variables]+1) 
    
    # set input variables (everything what is not response.variables and not "ID"):
    input.variables <- setdiff(names(dset), c(response.variables,ID.variable))
    
    #===============================================
    # PART 2b: DATA RECORD SELECTION
    #===============================================
    # disregard records which contain extreme values in response.variable (outliers)
    dset <- dset[dset[,response.variables]<opts$OCUT,] 
    opts$lim = c(min(dset[,response.variables],na.rm=T),
                 max(dset[,response.variables],na.rm=T))
    
    cat1(opts,opts$filename,":", length(dset[,1]), "records used.\n")

    #===============================================
    # PART 3 - 6
    #===============================================
    result <- tdmRegressLoop(dset,response.variables,input.variables,opts);

    # print summary output and attach certain columns (here: y,sd.y,dset) to list result:
    result <- tdmRegressSummary(result,opts,dset);
    
    tdmGraAndLogFinalize(opts,gdObj);      # close graphics and log file
    
    result;
    
}                                                   

readCmdCpu <- function(filename,opts) {
  read.csv2(file=paste(opts$dir.data, filename, sep=""), dec=".", na.string="-1",nrow=opts$READ.NROW);
}

cpu.postproc <- function(x,d,opts) { 
        x[x<0] <- 0; 		# clip all negative predictions to value 0	
        x;
}

#result = main_cpu()                            
