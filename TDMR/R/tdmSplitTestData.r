######################################################################################
# tdmSplitTestData
#
#'   Read and split the task data and return a \code{\link{TDMdata}} object.
#'
#'   Read the task data using \code{\link{tdmReadData}} and split them into a test part and 
#'   a training/validation-part.
#'
#'   @param opts a list from which we need here the elements
#'     \itemize{
#'       \item \code{READ.INI}:  [T] =T: do read and split, =F: return NULL   
#'       \item \code{READ.*}:  other settings for \code{\link{tdmReadData}}  
#'       \item \code{filename}:  needed for \code{\link{tdmReadData}}  
#'       \item \code{filetest}:  needed for \code{\link{tdmReadData}}  
#'       \item \code{TST.testFrac}:  [0.1] set this fraction of the daa aside for testing 
#'       \item \code{TST.COL}:   string with name for the partitioning column, if tdm$umode is not "SP_T".
#'                      (If tdm$umode=="SP_T", then TST.COL="tdmSplit" is used.)
#'     }
#'   @param tdm a list from  which we need here the elements 
#'     \itemize{
#'       \item \code{mainFile}:  if not NULL, set working dir to \code{dir(mainFile)} before executing  \code{\link{tdmReadData}} 
#'       \item \code{umode}:  [ "RSUB" | "CV" | "TST" | "SP_T" ], how to divide in training/validation data for tuning
#'                      and test data for the unbiased runs  
#'       \item \code{SPLIT.SEED}:  if not NULL, set random number generator to this value --> deterministic test set split
#'       \item \code{stratified}: [NULL] if not NULL, do the split by stratified sampling (at least one record of each class level
#'                      found in \code{dset[,tdm$stratified]} shall appear in the train-vali-set, recommended for classification)
#'     }
#'   @param nExp  [0] experiment counter, used to select a reproducible different seed, if tdm$SPLIT.SEED!=NULL
#'   @param dset  [NULL] if non-NULL, reading of dset is skipped and the given data frame dset is used.
#'
#'   @return \code{dataObj}, either NULL (if opts$READ.INI==FALSE) or an object of class \code{\link{TDMdata}} containing
#'      \item{dset}{ a data frame with the complete data set}
#'      \item{TST.COL}{ string, the name of the column in \code{dset} which has a 1 for 
#'                      records belonging to the test set and a 0 for train/vali records. If tdm$umode=="SP_T", then 
#'                      TST.COL="tdmSplit", else TST.COL=opts$TST.COL. }
#'      \item{filename}{ \code{opts$filename}, from where the data were read}
#'    Use the accessor functions  \code{\link{dsetTrnVa.TDMdata}} and \code{\link{dsetTest.TDMdata}} to extract the train/vali and 
#'    the test data, resp., from \code{dataObj}.
#'   
#'    Known caller: \code{\link{tdmBigLoop}}
#'
#' @seealso   \code{\link{dsetTrnVa.TDMdata}}, \code{\link{dsetTest.TDMdata}}, \code{\link{tdmReadData}}, \code{\link{tdmBigLoop}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), 2012 - 2013
#' @aliases TDMdata 
#' @export
######################################################################################
tdmSplitTestData <- function(opts,tdm,nExp=0,dset=NULL) {
	if (opts$READ.INI) {
  	if (exists(".Random.seed")) SAVESEED<-.Random.seed	   #save the Random Number Generator RNG status
		if (is.null(dset)) {
    		oldwd = getwd();                                              #
        if (!is.null(tdm$mainFile)) setwd(dirname(tdm$mainFile));     # save & change working dir

        dset <- tdmReadData(opts);
    
    		setwd(oldwd);       
    }                                          # restore working dir
  	
		if (tdm$umode=="SP_T") {              
        if (is.null(tdm$SPLIT.SEED)) {       
          set.seed(tdmRandomSeed());
        } else {
          newseed=tdm$SPLIT.SEED+nExp;
          set.seed(newseed);  # if you want reproducably the same training/test sets,
        }                     # but different for each experiment nExp

    		L = nrow(dset);
    		if (is.null(L)) stop("No data");    		if (L==0) stop("Empty data frame");
    		
        if (!is.null(tdm$stratified)) {
          # the division is done by ***stratified*** random sampling (recommended for classification):
      		cat1(opts,opts$filename,": Stratified random test-trainVali-index w.r.t. variable",tdm$stratified
              ,"and with tdm$TST.testFrac = ",tdm$TST.testFrac*100,"%\n");
      		if (!any(names(dset)==tdm$stratified)) stop("The value of tdm$stratified does not match any column name in dset!");
          rv <- dset[,tdm$stratified];
          #lrv = length(response.variables);
          #if (lrv>1) warning(sprintf("Stratified sampling is only done w.r.t. 1st response variable. It is not guaranteed to work for all %d response variables",lrv))  
          urv <- unique(rv);    # the class levels in rv
          # calculate tfr, the number of training set records for each level of variable tdm$stratified
          tfr <- sapply(urv,function(x) { round((1-tdm$TST.testFrac)*length(which(rv==x))) });
          # (this seems complicated, but the simpler command: tfr <- round((1-opts$TST.valiFrac)*table(rv));
          # does not work, because 'table' orders the levels alphabetically but 'unique'  requires them in the 
          # order they appear in column rv.) 
          #
          tfr[tfr<1] <- 1;      # ensure that there is at least one record for each class level in the trainVali-set
          cvi <- rep(1,L);
          for (i in 1:length(urv))  cvi[ sample(which(rv==urv[i]), tfr[i]) ] <- 0;
        } else {  # i.e. tdm$stratified is NULL
          # simple random sampling (recommended for regression):
          p <- sample(L)                  # random permutation of indices 1:L  
      		# calculate tfr, the record where the test set starts (TST.testFrac)
	        tfr <- (1-tdm$TST.testFrac)*L;
      		cat1(opts,opts$filename,": Setting data aside for testing with tdm$TST.testFrac = ",tdm$TST.testFrac*100,"%\n");
          cvi <- rep(0,L);
	        cvi[p[(tfr+1):L]] <- 1;          # test set index ( TST.testFrac  percent of the data)
        }
        
        cat1(opts,"*** from tdmSplitTestData: ***  \n");
        wI <- which(cvi==1);
        cat1(opts,"dset contains",L,"records; we put ",length(wI),"records aside into test set\n")
        if (length(wI)<300) print1(opts,wI);

    		TST.COL = "tdmSplit";
        if (any(names(dset)==TST.COL)) stop(sprintf("Name clash in dset, which has already a column \"%s\". Please consider renaming it.",TST.COL));            
    		dset[,TST.COL] = cvi;			       # add a new column "tdmSplit" with 1 for test data, 0 else
    		
		} else if (tdm$umode=="TST") {
    		TST.COL = opts$TST.COL;          # take the partition delivered by tdmReadData     
        if (!any(names(dset)==TST.COL)) {
          stop(sprintf("Data frame dset does not contain a column opts$TST.COL=\"%s\". \n%s",
                        opts$TST.COL,"This might be due to a missing opts$READ.TST==T when using tdm$umode==\"TST\".")); 
        }
		} else {  # i.e. if tdm$umode=="RSUB" or =="CV"
		    cvi = rep(0,nrow(dset));         # all data in dset used for training and validation
    		TST.COL = opts$TST.COL;          
        if (any(names(dset)==TST.COL)) stop(sprintf("Name clash in dset, which has already a column \"%s\". Please consider renaming it.",TST.COL));            
        dset[,TST.COL] = cvi;			       # add a new column with only 0's   --> dsetTest(dataObj) will return NULL
		}
		dataObj <- list(
			  dset=dset
			, TST.COL=TST.COL
			, filename=opts$filename
			, tdm = tdm
			);
		class(dataObj) <- "TDMdata";
   	if (exists("SAVESEED")) assign(".Random.seed", SAVESEED, envir=globalenv()); 		#load the saved RNG status

  	checkData(dataObj);
	} else {   # i.e. if opts$READ.INI==FALSE
		dataObj <- NULL;
	}	# if (opts$READ.INI)	
	
	dataObj;
}

# checkData:
#     Helper fct for tdmSplitTestData:
#     Perform some plausibility checks on the data read
#
checkData <- function(dataObj) {
  checkFactor <- function(set,txt,fact) {
      for (i in which(fact==TRUE)) {
        if (nlevels(set[,i])>32) {
          strng = sprintf("Column %s of %s has %d levels. Consider to use tdmPreGroupLevels() or as.numeric().",names(dset)[i], txt,nlevels(set[,i]));
          cat("NOTE:",strng);
          warning(strng); 
        }
      }
  }
  dset <- dsetTrnVa(dataObj);
  tset <- dsetTest(dataObj);
  
  dfactor <- sapply(dset,is.factor);
  tfactor <- sapply(tset,is.factor);
  if (any(dfactor!=tfactor)) {
    w = which(dfactor!=tfactor);
    strng = paste("dataObj has columns with different mode in train-vali- and test part:",sprintf("%s,",names(dset)[w]));
    cat("NOTE:",strng);
    warning(strng);
  }
  checkFactor(dset,"train-validation set",dfactor);
  checkFactor(tset,"test set",tfactor);
}

dsetTrnVa <- function(x)  UseMethod("dsetTrnVa");
dsetTrnVa.default <- function(x)  stop("Method dsetTrnVa only allowed for objects of class TDMdata. Consider opts$READ.INI=TRUE");
######################################################################################
# dsetTrnVa.TDMdata
#
#'   Return the train-validation part of a \code{\link{TDMdata}}  object containing the task data.
#'
#'   @method dsetTrnVa TDMdata
#'   @param x  return value from a prior call to \code{\link{tdmSplitTestData}}, an object of class \code{\link{TDMdata}}.
#'   @return \code{dset}, a data frame with all train-validation records
#'
#' @seealso   \code{\link{dsetTest.TDMdata}} \code{\link{tdmSplitTestData}}
#' @author Wolfgang Konen, FHK, Feb'2012
#' @export
#### @keywords internal
######################################################################################
dsetTrnVa.TDMdata <- function(x)  {
  if (!any(names(x$dset)==x$TST.COL)) {
    stop(sprintf("Data frame dset does not contain a column x$TST.COL=\"%s\". \n%s",
                 x$TST.COL,"Are you sure to have opts$READ.TST==T when using tdm$umode==\"TST\"?")); 
  }
  # return the training-validation  part of x$dset, excluding the column x$TST.COL (="tdmSplit" usually)
  ind=x$dset[,x$TST.COL]==0;
  if (x$tdm$umode == "SP_T") {
    x$dset[ind, -which(names(x$dset)==x$TST.COL) ];	
  } else {
    x$dset[ind,];
  }
  
}

dsetTest <- function(x)  UseMethod("dsetTest");
dsetTest.default <- function(x)  stop("Method dsetTest only allowed for objects of class TDMdata");
######################################################################################
# dsetTest.TDMdata
#
#'   Return the test part of a \code{\link{TDMdata}} object containing the task data.
#'
#'   @method dsetTest TDMdata
#'   @param x  return value from a prior call to \code{\link{tdmSplitTestData}}, an object of class \code{\link{TDMdata}}.
#'   @return \code{tset}, a data frame with all test records. If there are 0 test records, return NULL.
#'
#' @seealso   \code{\link{unbiasedRun}} \code{\link{dsetTrnVa.TDMdata}} 
#' @author Wolfgang Konen, FHK, Feb'2012
#' @export
#### @keywords internal
######################################################################################
dsetTest.TDMdata <- function(x)  {
  ind = which(x$dset[,x$TST.COL]==1);
  if (length(ind)==0) {
    NULL;
  } else {
    # return the test part of data frame dset, excluding the column x$TST.COL (="tdmSplit" usually)
    if (x$tdm$umode == "SP_T") {
      x$dset[ind, -which(names(x$dset)==x$TST.COL) ];	
    } else {
      x$dset[ind,];
    }
  }
}

######################################################################################
# print.TDMdata
#
#'   Print an overview for a \code{\link{TDMdata}} object.
#'
#'   Print number of rows and number of columns of the data frame \code{dset} contained in the \code{\link{TDMdata}} object.
#'
#'   @method print TDMdata
#'   @param x  return value from a prior call to \code{\link{tdmSplitTestData}}, an object of class \code{\link{TDMdata}}.
#'   @param ... currently not used
#'
#' @seealso   \code{\link{tdmSplitTestData}}
#' @author Wolfgang Konen, FHK
#' @export
######################################################################################
print.TDMdata <- function(x,...) {
  cat(sprintf("TDMdata object with %d records and %d variables:\n",nrow(x$dset),length(x$dset)));
  print(names(x$dset));
}