######################################################################################
# tdmRegressLoop
#
#'    Core regression double loop of TDMR returning a \code{\link{TDMregressor}} object. 
#'
#'    tdmRegressLoop contains a double loop (opts$NRUN and CV-folds)
#'    and calls \code{\link{tdmRegress}}. It is called  by all R-functions main_*. \cr
#'    It returns an object of class \code{\link{TDMregressor}}.
#'
#'   @param dset    the data frame for which cvi is needed
#'   @param response.variables   name of column which carries the target variable - or -
#'                   vector of names specifying multiple target columns
#'                   (these columns are not used during prediction, only for training and
#'                   for evaluating the predicted result)
#'   @param input.variables     vector with names of input columns
#'   @param opts    a list from which we need here the following entries
#'     \describe{
#'       \item{\code{NRUN}}{ number of runs (outer loop)}
#'       \item{\code{TST.SEED}}{ =NULL: leave the random number seed as it is. =any value: set the random number seed
#'               to this value to get reproducible random numbers and thus reproducible training-test-set-selection.
#'               (only relevant in case TST.kind=="cv" or "rand") (see also MOD.SEED in \code{\link{tdmClassify}})}
#'       \item{\code{TST.kind}}{ how to create cvi, handed over to \code{\link{tdmModCreateCVindex}}. If TST.kind="col", then cvi is taken from dset[,opts$TST.col].}
#'       \item{\code{GD.RESTART}}{ [TRUE] =TRUE/FALSE: do/don't restart graphic devices}
#'       \item{\code{GRAPHDEV}}{ ["non"| other ]}
#'     }
#'
#'   @return \code{result},  an object of class \code{\link{TDMregressor}}, this is a list with results, containing
#     \describe{
#'       \item{opts}{ the res$opts from \code{\link{tdmRegress}}}
#'       \item{lastRes}{ last run, last fold: result from \code{\link{tdmRegress}}}
#'       \item{R_train}{ RMAE / RMSE on training set (vector of length NRUN), depending on opts$rgain.type=="rmae" or "rmse"}
#'       \item{S_train}{ RMSE on training set (vector of length NRUN)}
#'       \item{T_train}{ Theil's U for RMAE on training set (vector of length NRUN)}
#'       \item{*_test}{ --- similar, with test set instead of training set ---  }
#'       \item{Err}{ a data frame with as many rows as opts$NRUN and columns = (rmae.trn, rmse.trn
#'                   made.trn, rmae.theil.trn, ntrn, rmae.tst, rmse.tst, made.tst, rmae.theil.tst,
#'                   ntst) }
#'       \item{predictions}{ last run: data frame with dimensions [nrow(dset),length(response.variable)]. In case of CV, all 
#'              validation set predictions (for each record in dset), in other cases mixed validation / train set predictions.  }
#    }
#'
#' @seealso   \code{\link{tdmRegress}}, \code{\link{tdmClassifyLoop}}, \code{\link{tdmClassify}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), FHK, Sep'2010 - Jun'2012
#' @aliases TDMregressor 
#' @example demo/demo00-1regress.r
#' @export
######################################################################################
tdmRegressLoop <- function(dset,response.variables,input.variables,opts) {
  	if (exists(".Random.seed")) SAVESEED<-.Random.seed	   #save the Random Number Generator RNG status
    if (class(opts)[1] != "tdmOpts")  stop("Class of object opts is not tdmOpts");
    if (is.null(opts$PRE.PCA.numericV)) opts$PRE.PCA.numericV <- input.variables;
    if (opts$NRUN<=0) stop(sprintf("opts$NRUN has to be positive, but it is %d",opts$NRUN));
    if (opts$PRE.PCA!="none" & opts$PRE.SFA!="none") stop("It is not allowed to activate opts$PRE.PCA and opts$PRE.SFA simultaneously.")
    rgain.types = c("rmae","rmse","made");
    if (! (opts$rgain.type %in% rgain.types)) {
      warning(sprintf("opts$rgain.type='%s' has not one of the values {%s} allowed for regression. TDMR will set it to 'rmae'.",
                      opts$rgain.type,paste(rgain.types,collapse=",")));
      opts$rgain.type="rmae";
    }
    if (!all(response.variables %in% names(dset)))
      stop(sprintf("Not all response.variables are in names(dset)!\n  %s\n  response.variables=%s",
                   "Note that response.variables have to be strings (names of columns in dset), not the columns themselves.",
                   paste(response.variables,collapse=","))
           )
    
    R_train <- R_vali <- T_train <- T_test <- NULL       # reserve names (dynamic extension of
    S_train <- S_test <- NULL                            # these vectors at and of for-i-loop)
    Err <- NULL
    for (i in 1:opts$NRUN) {
        if (opts$NRUN>1) {
          if (opts$GD.RESTART) tdmGraphicCloseDev(opts);
          tdmGraphicInit(opts);
        }
        opts$i = i;

        #=============================================================================
        # PART 3: CREATE NFOLD CROSSVALIDATION INDEX  OR DIVIDE IN TRAINING / TEST SET
        #=============================================================================
        if (is.null(opts$TST.SEED)) {
          # NEW: when called via SPOT, the RNG might be at (different but) fixed seed in each call.
          #      But we want different seeds (for test set selection) to see the variability
          set.seed(tdmRandomSeed());
        } else if (opts$TST.SEED=="algSeed") {    # use the seed from SPOT
          newseed=opts$ALG.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed); 
        } else {
          newseed=opts$TST.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed);  # if you want reproducably the same training/test sets,
        }                     # but different for each run i and each repeat (opts$rep)

        cvi <- tdmModCreateCVindex(dset,response.variables,opts);
        nfold = max(cvi,1);     # Why ',1'? - In the special case where all records belong to train sample (max(cvi)=0), we want 
                                # to have exactly one fold  (one pass through k-loop) and NOT k in { 1, 0 }

        #=============================================
        # PART 4: MODELING, EVALUATION, VISUALIZATION
        #=============================================
        alltrn = alltst = NULL;
        predictions <- as.data.frame(0*dset[,response.variables]);
        names(predictions) <- response.variables;
        for (k in 1:nfold) {
            opts$k=k;
            opts$the.nfold = nfold;
            d_test  <- dset[cvi==k, ];             # validation set
            d_train <- dset[cvi!=k & cvi>=0, ];    # training set (disregard records with cvi<0)
            ntst=nrow(d_test);
            ntrn=nrow(d_train);

            d_preproc <- NULL;                # no preprocessing -> no mem-consumption for d_preproc  (if condition on next line is false)
            if (opts$PRE.PCA!="none" | opts$PRE.SFA!="none") {
              if (opts$PRE.allNonVali) {
                d_preproc <- dset[cvi!=k, ]   # all non-validation-data are used for preproc
              } else {
                d_preproc <- d_train 
              }
            } 
            if (opts$PRE.PCA!="none") {
              # a) do PCA on the numeric variables of d_train, if opts$PRE.PCA!="none"
              # b) add monomials of degree 2 for the first opts$PRE.PCA.npc numeric variables
              # c) apply this PCA and monomials to d_test and d_dis in the same way
              other.variables <- setdiff(input.variables,opts$PRE.PCA.numericV);
              pca <- tdmPrePCA.train(d_train,opts);                 # see tdmPreprocUtils.r
              d_train <- tdmPrePCA.apply(d_train,pca$pcaList,opts,d_train)$dset;
              d_test <- tdmPrePCA.apply(d_test,pca$pcaList,opts,d_train)$dset;
              d_dis <- tdmPrePCA.apply(d_dis,pca$pcaList,opts,d_train)$dset;

              input.variables <- union(pca$numeric.variables,other.variables);
              #opts$PRE.SFA.numericV <- pca$numeric.variables;        # -- now obsolete, no SFA after PCA allowed --
              if (length(setdiff(input.variables,names(d_train)))>0) 
                  stop("Some elements of input.variables are not columns of d_train");
            }

            res <- tdmRegress(d_train,d_test,d_preproc,response.variables,input.variables,opts)
            #res <- regress_lm(d_train,d_test,response.variables,input.variables,opts)

            avgRMAE = as.data.frame(t(colSums(res$allRMAE)))/nrow(res$allRMAE)
            avgRMSE = as.data.frame(t(colSums(res$allRMSE)))/nrow(res$allRMSE)
            cat1(opts,sprintf("k=%d  RMAE_train=%7.5f  RMAE_test=%7.5f\n",k,avgRMAE$train,avgRMAE$test))
            alltrn = rbind(alltrn,as.data.frame(list(rmae.trn=avgRMAE$rmae.train * ntrn
                                                    ,rmse.trn=avgRMSE$rmse.train * ntrn
                                                    ,made.trn=avgRMAE$made.tr * ntrn
                                                    ,rmae.theil.trn=avgRMAE$theil.train * ntrn
                                                    ,ntrn=ntrn
                                                    )));
            alltst = rbind(alltst,as.data.frame(list(rmae.tst=avgRMAE$rmae.test * ntst
                                                    ,rmse.tst=avgRMSE$rmse.test * ntst
                                                    ,made.tst=avgRMAE$made.te * ntst
                                                    ,rmae.theil.tst=avgRMAE$theil.test * ntst
                                                    ,ntst=ntst
                                                    )));
            predictions[cvi==k,response.variables] <- res$d_test[,paste("pred_",response.variables,sep="")];
            if (!(opts$TST.kind=="cv")) {
              predictions[cvi!=k & cvi>=0,response.variables] <- res$d_train[,paste("pred_",response.variables,sep="")];
            } 
        } # for (k)
        
        # each measure is a weighted average over folds, weighted with #cases in each fold:
        Err = rbind(Err,c(colSums(alltrn)/sum(alltrn$ntrn),colSums(alltst)/sum(alltst$ntst)));
        Err[i,"ntrn"]=ifelse(opts$TST.kind=="cv",nrow(dset),nrow(d_train));   # the number of distinct cases 
        Err[i,"ntst"]=ifelse(opts$TST.kind=="cv",nrow(dset),nrow(d_test));    # used for this line of measures
        
        col.trn = paste(opts$rgain.type,".trn",sep="");
        col.tst = paste(opts$rgain.type,".tst",sep="");
        cat1(opts,"\n",ifelse(opts$TST.kind=="cv","CV",""),opts$rgain.string,"on training set   ",Err[i,col.trn],"\n")
        cat1(opts,"",  ifelse(opts$TST.kind=="cv","CV",""),opts$rgain.string,"on     test set   ",Err[i,col.tst],"\n\n")

        R_train[i] = Err[i,col.trn];
        #if (opts$rgain.type=="rmse") R_train[i] = Err["rmse.trn"];
        #if (opts$rgain.type=="made") R_train[i] = Err["made.trn"];
        S_train[i] = Err[i,"rmse.trn"];
        T_train[i] = Err[i,"rmae.theil.trn"];
        R_vali[i] = Err[i,col.tst];
        #if (opts$rgain.type=="rmse") R_vali[i] = Err["rmse.tst"];
        #if (opts$rgain.type=="made") R_vali[i] = Err["made.tst"];
        S_test[i] = Err[i,"rmse.tst"]
        T_test[i] = Err[i,"rmae.theil.tst"]

        #=============================================
        # PART 5: SOME MORE GRAPHICS
        #=============================================
        if (opts$GD.DEVICE!="non" & !is.null(opts$gr.fctcall)) {
          # execute the graphics command given in text string opts$gr.fctcall
          eval(parse(text=opts$gr.fctcall));
        }
    } # for (i in 1:opts$NRUN)

    #=============================================
    # PART 6: OVERALL EVALUATION
    #=============================================
    if (opts$NRUN>1) {
        # print output averaged over all NRUN runs "for (i)"
        # Expected result: rmse$test & rmse$train should approach same value
        cat1(opts,"\nAverage over all ",opts$NRUN," runs: \n")
        cat1(opts,sprintf("RMAE.trn: (%7.5f +- %7.5f)%%\n", mean(R_train)*100, sd(R_train)*100));
        cat1(opts,sprintf("RMAE.tst: (%7.5f +- %7.5f)%%\n", mean(R_vali)*100, sd(R_vali)*100));
        cat1(opts,sprintf("Theil.train: (%7.2f +- %4.2f)%%\n", mean(T_train), sd(T_train)));
        cat1(opts,sprintf("Theil.test: (%7.2f +- %4.2f)%%\n", mean(T_test), sd(T_test)));
        cat1(opts,sprintf("RMSE.train: (%7.2f +- %4.2f)%%\n", mean(S_train), sd(S_train)));
        cat1(opts,sprintf("RMSE.test: (%7.2f +- %4.2f)%%\n", mean(S_test), sd(S_test)));
    }

    result = list(lastRes = res     # last run, last fold: result from tdmRegress
              	#, opts = res$opts    # deprecated (12/2011), use result$lastRes$opts or Opts(result)
              	, R_train = R_train
              	, R_vali = R_vali
              	, T_train = T_train
              	, T_test = T_test
              	, S_train = S_train
              	, S_test = S_test
              	, Err = Err
              	, predictions = predictions
              	);
    class(result) <- c("TDMregressor", "TDM")     # this causes > result; or > print(result);
                                                  # NOT to print out the whole list (might be very long!!)
                                                  # but to call instead the function  print.TDMregressor
                                                  # (which in turn calls tdmRegressSummary)
   	if (exists("SAVESEED")) assign(".Random.seed", SAVESEED, envir=globalenv()); 		#load the saved RNG status
    result;
} # tdmRegressLoop

######################################################################################
# tdmRegressSummary
#
#'   Print summary output for \code{result} from \code{tdmRegressLoop} and add \code{result$y}.
#'
#'   \code{result$y} is "OOB RMAE" on training set for methods RF or MC.RF.
#'   \code{result$y} is "RMAE" on test set (=validation set) for all other methods.
#'   \code{result$y} is the quantity which the tuner seeks to minimize.
#'
#'   @param result  return value from a prior call to \code{\link{tdmRegressLoop}}, an object of class \code{TDMregressor}.
#'   @param opts    a list from which we need here the following entries
#'     \describe{
#'       \item{\code{NRUN}}{ number of runs (outer loop)}
#'       \item{\code{method}}{}
#'       \item{\code{VERBOSE}}{}
#'       \item{\code{dset}}{ [NULL] if !=NULL, attach it to result}
#'     }
#'   @param dset    [NULL] if not NULL, add this data frame to the return value (may cost a lot of memory!)
#'
#'   @return \code{result},  an object of class \code{TDMregressor}, with \code{result$y}, \code{result$sd.y}
#'          (and optionally also \code{result$dset}) added
#'
#' @seealso   \code{\link{tdmRegress}}, \code{\link{tdmRegressLoop}}, \code{\link{tdmClassifySummary}}
#' @author Wolfgang Konen, FHK, Sep'2010 - Oct'2011
#' @export
######################################################################################
tdmRegressSummary <- function(result,opts,dset=NULL)
{
    res <- result$lastRes;
    cat1Records <- function (nrow_noCV) {
      cat1(opts,ifelse(opts$TST.kind=="cv"
                ,  sprintf("   (on %d records in %d folds)",nrow(res$d_train)+nrow(res$d_test),opts$TST.NFOLD)
                ,  sprintf("   (on %d records)",nrow_noCV)
                ),"\n");
    }
    #print2(opts,res$allRMAE);		   # RMAE for each response variable , but only for lastRes
    y = mean(result$R_vali);       # RMAE, average of opts$NRUN runs
    ytr = mean(result$R_train);
    if (opts$MOD.method %in% c("RF","MC.RF")) {
      cat1(opts,sprintf("\n%sTrain OOB RMAE: %7.3f",ifelse(opts$TST.kind=="cv","CV ",""),ytr));
      cat1(opts,ifelse(opts$NRUN>1,sprintf(" +-%7.3f",sd(result$R_train)),""));
      cat1Records(nrow(res$d_train));
      result$y=ytr;           # the score (to be minimized by SPOT) is "RMAE OOB"
      result$sd.y=sd(result$R_train);
    } else {
      result$y=y;             # the score (to be minimized by SPOT) is "RMAE test set"
      result$sd.y=sd(result$R_vali);
    }
    cat1(opts,sprintf("%s  Vali    RMAE: %7.3f",ifelse(opts$TST.kind=="cv","CV ",""),y));
    cat1(opts,ifelse(opts$NRUN>1,sprintf(" +-%7.3f",sd(result$R_vali)),""));
    cat1Records(nrow(res$d_test));

    if (!is.null(dset)) result$dset=dset;          # might cost a lot of memory

    result;

} # tdmRegressSummary

