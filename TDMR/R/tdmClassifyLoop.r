######################################################################################
# tdmClassifyLoop:
#
#'    Core classification double loop of TDMR returning a \code{\link{TDMclassifier}} object. 
#'    
#'    tdmClassifyLoop contains a double loop (opts$NRUN and CV-folds)
#'    and calls \code{\link{tdmClassify}}. It is called  by all classification R-functions main_*. \cr
#'    It returns an object of class \code{\link{TDMclassifier}}.
#'
#'
#'   @param dset    the data frame containing training and validation data.
#'   @param tset    [NULL] If not NULL, this is the validation data set. If NULL, the validation data 
#'                  set is build from \code{dset} according to the procedure prescribed in \code{opts$TST.*}. 
#'   @param response.variables   name of column which carries the target variable - or -
#'                   vector of names specifying multiple target columns
#'                   (these columns are not used during prediction, only for evaluation)
#'   @param input.variables     vector with names of input columns
#'   @param opts    a list from which we need here the following entries
#'     \describe{
#'       \item{\code{NRUN}}{ number of runs (outer loop)}
#'       \item{\code{TST.SEED}}{ =NULL: get a new random number seed with \code{\link{tdmRandomSeed}}. =any value: set the random number seed
#'               to this value to get reproducible random numbers and thus reproducible training-test-set-selection.
#'               (only relevant in case TST.kind=="cv" or "rand") (see also MOD.SEED in \code{\link{tdmClassify}})}
#'       \item{\code{TST.kind}}{ how to create cvi, handed over to \code{\link{tdmModCreateCVindex}}. If TST.kind="col", then cvi is taken from dset[,opts$TST.col].}
#'       \item{\code{GD.RESTART}}{ [TRUE] =TRUE/FALSE: do/don't restart graphic devices}
#'       \item{\code{GD.DEVICE}}{ ["non"|"win"|"pdf"|"png"]}
#'     }
#'
#'   @return \code{result},  an object of class \code{\link{TDMclassifier}}, this is a list with results, containing
#'       \item{lastRes}{ last run, last fold: result from \code{\link{tdmClassify}}}
#'       \item{C_train}{ classification error on training set}
#'       \item{G_train}{ gain on training set}
#'       \item{R_train}{ relative gain on training set (percentage of max. gain on this set)}
#'       \item{*_vali}{ --- similar, with vali set instead of training set ---   }
#'       \item{*_vali2}{ --- similar, with vali2 set instead of training set ---  }
#'       \item{Err}{ a data frame with as many rows as opts$NRUN and 9 columns corresponding to 
#' 				the nine variables described above }
#'       \item{predictions}{ last run: data frame with dimensions [nrow(dset),length(response.variable)]. In case of CV, all 
#'              CV predictions (for each record in dset), in other cases mixed validation / train set predictions.  }
#'       \item{predProbList}{ the ith element in this list has data on the prediction probabilities of the ith run. See info on 
#'              \code{predProb} in \code{\link{tdmClassify}}. }
#'
#'    Each performance measure \code{C_*, G_*, R_*} is a vector of length \code{opts$NRUN}. To be specific, \code{C_train[i]} is the
#'    classification error on the training set from the \code{i}-th run. This error is \code{mean(res$allEVAL$cerr.trn)}, i.e. the
#'    mean of the classification errors from all response variables when \code{res} is the return value of  \code{\link{tdmClassify}}.
#'    In the case of cross validation, for each performance measure an additional averaging over all folds is done.
#'
#' @seealso   \code{\link{print.TDMclassifier}}, \code{\link{tdmClassify}}, \code{\link{tdmRegress}}, \code{\link{tdmRegressLoop}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), FHK, Sep'2010 - Jun'2012
#' @aliases TDMclassifier 
#' @example demo/demo00-0classif.r
#' @export
######################################################################################
tdmClassifyLoop <- function(dset,response.variables,input.variables,opts,tset=NULL) {
  	if (exists(".Random.seed")) SAVESEED<-.Random.seed	   #save the Random Number Generator RNG status
#print(.Random.seed[1:6])
    if (class(opts)[1] != "tdmOpts")  stop("Class of object opts is not tdmOpts");
    if (is.null(opts$PRE.PCA.numericV)) opts$PRE.PCA.numericV <- input.variables;
    if (opts$NRUN<=0) stop(sprintf("opts$NRUN has to be positive, but it is %d",opts$NRUN));
    if (opts$PRE.PCA!="none" & opts$PRE.SFA!="none") stop("It is not allowed to activate opts$PRE.PCA and opts$PRE.SFA simultaneously.")
 
    #if (opts$READ.TST==TRUE & opts$TST.kind!="col")
    #  warning(sprintf("Are you sure you want opts$READ.TST=TRUE, but opts$TST.kind!='col'? Actual value is opts$TST.kind='%s'.",opts$TST.kind));
    
  	if (!all(response.variables %in% names(dset)))
  	  stop(sprintf("Not all response.variables are in names(dset)!\n  %s\n  response.variables=%s",
  	               "Note that response.variables have to be strings (names of columns in dset), not the columns themselves.",
  	               paste(response.variables,collapse=",")))

    predProbList=list();
    C_train <- C_vali <- C_vali2 <- G_train <- G_vali <- NULL # reserve names (dynamic extension of
    R_train <- R_vali <- R_vali2 <- G_vali2 <- NULL           # these vectors at and of for-i-loop)
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
        } else if (opts$TST.SEED=="algSeed") {    # use the seed from SPOT:
          # opts$ALG.SEED is set in tdmStartSpot to des$SEED[k]+r. This meens that the overall r'th 
          # evaluation of a design point gets the seed spotConfig$alg.seed+r
          newseed=opts$ALG.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed); 
        } else {
          newseed=opts$TST.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed);  # if you want reproducably the same training/test sets,
        }                     # but different for each run i and each repeat (opts$rep)
#print(.Random.seed[1:6])

        cvi <- tdmModCreateCVindex(dset,response.variables,opts,stratified=TRUE);    
        nfold = max(cvi,1);     # Why ',1'? - In the special case where all records belong to train sample (max(cvi)=0), we want 
                                # to have exactly one fold  (one pass through k-loop) and NOT k in { 1, 0 }

#print(cvi[1:20])
#browser();
        #=============================================
        # PART 4: MODELING, EVALUATION, VISUALIZATION
        #=============================================
        alltrn = alltst = NULL;
        predictions <- as.data.frame(dset[,response.variables]);
        names(predictions) <- response.variables;
        predProb=list();
        for (k in 1:nfold) {
            opts$k=k;
            opts$the.nfold = nfold;
            cat1(opts,"\n")
            if (any(names(dset)=="IND.dset")) stop("Name clash in dset, which has already a column IND.dset. Please consider renaming it.")            
            if (is.null(tset)) { 
              d_test  <- dset[cvi==k, ]            # validation set
              d_test  <- tdmBindResponse(d_test , "IND.dset", which(cvi==k));
            } else {    # i.e. if (!is.null(tset))
              # validation set (or test set) is passed in from the caller
              if (opts$TST.kind %in% c("cv"))  #,"col"))
                stop(sprintf("Option opts$TST.kind=\"%s\" together with tset!=NULL is currently not implemented. Consider opts$TST.kind=\"rand\".",opts$TST.kind));
              d_test <- tset;
              d_test  <- tdmBindResponse(d_test , "IND.dset", rep(-1,nrow(tset)));
            }
            d_train <- dset[cvi!=k & cvi>=0, ]   # training set (disregard records with cvi<0)
            d_dis   <- dset[cvi!=k & cvi==-1, ]  # disregard set (needed for active learning)
            d_train <- tdmBindResponse(d_train, "IND.dset", which(cvi!=k & cvi>=0));
            d_dis   <- tdmBindResponse(d_dis  , "IND.dset", which(cvi!=k & cvi==-1));
            ntst=nrow(d_test);
            ntrn=nrow(d_train);
            #if (nrow(d_train)+nrow(d_test)+nrow(d_dis) != nrow(dset))
            #  stop("Something wrong, the union of d_train, d_test and d_dis does not cover dset");
 
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
              pca <- tdmPrePCA.train(d_preproc,opts);                 # see tdmPreprocUtils.r
              d_train <- tdmPrePCA.apply(d_train,pca$pcaList,opts,d_train)$dset;
              d_test <- tdmPrePCA.apply(d_test,pca$pcaList,opts,d_train)$dset;
              d_dis <- tdmPrePCA.apply(d_dis,pca$pcaList,opts,d_train)$dset;
              input.variables <- union(pca$numeric.variables,other.variables);
              #opts$PRE.SFA.numericV <- pca$numeric.variables;        # -- now obsolete, no SFA after PCA allowed --
              if (length(setdiff(input.variables,names(d_train)))>0) 
                  stop("Some elements of input.variables are not columns of d_train");
            }

            # the SFA preprocessing is now in tdmClassify, inside response.variable for-loop
            # (because SFA (for classification) depends on response.variable)

            res <- tdmClassify(d_train,d_test,d_dis,d_preproc,response.variables,input.variables,opts)

            # predProb$Val: bind the different folds of CV (cross validation) together. If no CV, then nfold=1 ->> only results from one vali set.
            predProb$Val = rbind(predProb$Val,res$predProb$Val);
            # predProb$Trn: take only the results from the first fold
            if (k==1) predProb$Trn = res$predProb$Trn;
            
            alltrn = rbind(alltrn,as.data.frame(list(cerr.trn=mean(res$allEVAL$cerr.trn) * ntrn
                                                    ,gain.trn=mean(res$allEVAL$gain.trn) * ntrn
                                                    ,rgain.trn=mean(res$allEVAL$rgain.trn) * ntrn
                                                    ,ntrn=ntrn
                                                    )));
            alltst = rbind(alltst,as.data.frame(list(cerr.tst=mean(res$allEVAL$cerr.tst) * ntst
                                                    ,gain.tst=mean(res$allEVAL$gain.tst) * ntst
                                                    ,rgain.tst=mean(res$allEVAL$rgain.tst) * ntst
                                                    ,cerr.tst2=mean(res$allEVAL$cerr.tst2) * ntst
                                                    ,gain.tst2=mean(res$allEVAL$gain.tst2) * ntst
                                                    ,rgain.tst2=mean(res$allEVAL$rgain.tst2)
                                                    ,ntst=ntst
                                                    )));
            if (is.null(tset)) {
              predictions[cvi==k,response.variables] <- res$d_test[,paste("pred_",response.variables,sep="")];
              predictTest = NULL;
            } else {
              predictTest = res$d_test[,paste("pred_",response.variables,sep="")];
            }
            if (!(opts$TST.kind=="cv")) {
              predictions[cvi!=k & cvi>=0,response.variables] <- res$d_train[1:ntrn,paste("pred_",response.variables,sep="")];
                                                                     # Why '1:ntrn'? - If opts$ncopies>0, then tdmClassify has added
                                                                     # opts$ncopies bootstrap patterns to the ntrn original patterns.
                                                                     # These bootstrap patterns have to be discarded here.
            } 
        }  # for (k in 1:nfold)

        ### -- OLD: Err = colSums(allerr)/nfold     # assumes that each fold has the same (or nearly the same) size
        # each measure is a weighted average over folds, weighted with #cases in each fold:
        Err = rbind(Err,c(colSums(alltrn)/sum(alltrn$ntrn),colSums(alltst)/sum(alltst$ntst)));
        Err[i,"ntrn"]=ifelse(opts$TST.kind=="cv",nrow(dset),nrow(d_train));   # the number of distinct cases 
        Err[i,"ntst"]=ifelse(opts$TST.kind=="cv",nrow(dset),nrow(d_test));    # used for this line of measures

        predProbList[[i]] <- list()
        predProbList[[i]]$Val <- predProb$Val;
        predProbList[[i]]$Trn <- predProb$Trn;

        cat1(opts,"\n",ifelse(opts$TST.kind=="cv","CV",""),"Relative gain on   training set   ",Err[i,"rgain.trn"],"%\n")
        cat1(opts,"",  ifelse(opts$TST.kind=="cv","CV",""),"Relative gain on validation set   ",Err[i,"rgain.tst"],"%\n\n")

        #opts$name.of.prediction <- paste("pred_", response.variable, sep="")
        C_train[i] = Err[i,"cerr.trn"]
        G_train[i] = Err[i,"gain.trn"]
        R_train[i] = Err[i,"rgain.trn"]
        C_vali[i] = Err[i,"cerr.tst"]
        G_vali[i] = Err[i,"gain.tst"]
        R_vali[i] = Err[i,"rgain.tst"]
        C_vali2[i] = Err[i,"cerr.tst2"]
        G_vali2[i] = Err[i,"gain.tst2"]             # for comparision in case of opts$MOD.method="RF": results with
        R_vali2[i] = Err[i,"rgain.tst2"]            # default cutoff 1/n.class instead of opts$CLS.cutoff

        #=============================================
        # PART 5: SOME MORE GRAPHICS
        #=============================================
        if (opts$GD.DEVICE!="non" & !is.null(opts$gr.fctcall)) {
          # execute the graphics command given in text string opts$gr.fctcall
          eval(parse(text=opts$gr.fctcall));
        }
    } # for (i in 1:opts$NRUN)

    #write.table(d_test, file=paste(opts$dir.output, sub(".csv", "", filename), "_predictions.csv", sep=""), quote=FALSE, sep=";", dec=".", row.names=FALSE, col.names=TRUE)

    #=============================================
    # PART 6: OVERALL EVALUATION
    #=============================================
    if (opts$NRUN>1) {
        # print output averaged over all NRUN runs "for (i)"
        cat1(opts,"\nAverage over all ",opts$NRUN," runs: \n")
        cat1(opts,sprintf("cerr$train: (%7.5f +- %7.5f)%%\n", mean(C_train)*100, sd(C_train)*100));
        cat1(opts,sprintf("cerr$vali:  (%7.5f +- %7.5f)%%\n", mean(C_vali)*100, sd(C_vali)*100));
        cat1(opts,sprintf("gain$train: (%7.2f +- %4.2f)\n", mean(G_train), sd(G_train)));
        cat1(opts,sprintf("gain$vali:  (%7.2f +- %4.2f)\n", mean(G_vali), sd(G_vali)));
        cat1(opts,sprintf("rgain.train: %7.3f%%\n", mean(R_train)));
        cat1(opts,sprintf("rgain.vali:  %7.3f%%\n\n", mean(R_vali)));
    }
    result = list(lastRes = res     # last run, last fold: result from tdmClassify
                	#, opts = res$opts   # deprecated (12/2011), use result$lastRes$opts or Opts(result)
                	, C_train = C_train
                	, G_train = G_train
                	, R_train = R_train
                	, C_vali = C_vali
                	, G_vali = G_vali
                	, R_vali = R_vali
                	, C_vali2 = C_vali2
                	, G_vali2 = G_vali2
                	, R_vali2 = R_vali2
                	, Err = Err
                	, predictions = predictions     # predictions on dset from the *last* run: all CV predictions in case of CV, 
                                                  # mixed vali-train predictions else
                  , predictTest = predictTest
                	, predProbList = predProbList   # ith element: a list with predProb data from the ith run 
                	);
    if (!is.null(opts$TST.COL))         # needed when result is input for coTraining  (see ssl_methods.r)
      if (opts$TST.COL %in% names(dset))  result$TST = dset[,opts$TST.COL]

    class(result) <- c("TDMclassifier","TDM")     # this causes > result; or > print(result);
                                                  # NOT to print out the whole list (might be very long!!)
                                                  # but to call instead the function  print.TDMclassifier
                                                  # (which in turn calls tdmClassifySummary)
   	if (exists("SAVESEED")) assign(".Random.seed", SAVESEED, envir=globalenv()); 		#load the saved RNG status
    result;
} # tdmClassifyLoop

######################################################################################
# tdmClassifySummary
#
#'   Print summary output for \code{result} from \code{tdmClassifiyLoop} and add \code{result$y}.
#'
#'   \code{result$y} is "minus OOB rgain" on training set for methods RF or MC.RF.
#'   \code{result$y} is "minus rgain" on test set (=validation set) for all other methods.
#'   \code{result$y} is the quantity which the tuner seeks to minimize.
#'
#'   @param result  return value from a prior call to \code{\link{tdmClassifyLoop}}, an object of class \code{TDMclassifier}.
#'   @param opts    a list from which we need here the following entries
#'     \describe{
#'       \item{\code{NRUN}}{ number of runs (outer loop)}
#'       \item{\code{method}}{}
#'       \item{\code{VERBOSE}}{}
#'       \item{\code{dset}}{ [NULL] if !=NULL, attach it to result}
#'     }
#'   @param dset    [NULL] if not NULL, add this data frame to the return value (may cost a lot of memory!)
#'
#'   @return \code{result},  an object of class \code{TDMclassifier}, with \code{result$y}, \code{result$sd.y}
#'          (and optionally also \code{result$dset}) added
#'
#' @seealso   \code{\link{tdmClassify}}, \code{\link{tdmClassifyLoop}}, \code{\link{print.TDMclassifier}}, \code{\link{tdmRegressSummary}}
#' @author Wolfgang Konen, FHK, Sep'2010 - Oct'2011
#' @export
######################################################################################
tdmClassifySummary <- function(result,opts,dset=NULL)
{
    res <- result$lastRes;
    cat1Records <- function (nrow_noCV) {
      cat1(opts,ifelse(opts$TST.kind=="cv"
                ,  sprintf("   (on %d records in %d folds)",nrow(res$d_train)+nrow(res$d_test),opts$TST.NFOLD)
                ,  sprintf("   (on %d records)",nrow_noCV)
                ),"\n");
    }
    #print2(opts,res$allEVAL);		      # EVAL for each response variable , but only for lastRes
    y = mean(result$R_vali);
    ytr = mean(result$R_train);
    maxScore    = result$G_vali[1]/(result$R_vali[1]/100);
    maxScore.tr = result$G_train[1]/(result$R_train[1]/100);
    z=data.frame(TYPE=c("rgain","meanCA","minCA")
                ,DESC=c("relative gain, i.e. percent of correctly classified records"
                        ,"mean class accuracy, i.e. average over class levels","minimum class accuracy"));
    cat1(opts,sprintf("\nRelative gain is \"%s\"",opts$rgain.type));
    cat1(opts,sprintf(" (%s)", z$DESC[which(z$TYPE==opts$rgain.type)]));
    if (opts$MOD.method %in% c("RF","MC.RF") & opts$RF.OOB==TRUE) {
      cat1(opts,sprintf("\n%sTrain OOB relative gain: %7.3f",ifelse(opts$TST.kind=="cv","CV ",""),ytr));
      cat1(opts,ifelse(opts$NRUN>1,sprintf(" +-%7.3f",sd(result$R_train)),""));
      cat1Records(nrow(res$d_train));
#      cat1(opts,sprintf("   (on %d records)",nrow(res$d_train)));
      result$y=-ytr;           # the score (to be minimized by SPOT) is "minus OOB score"
      result$sd.y=sd(result$R_train);
    } else {
      cat1(opts,"\n");
      result$y=-y;             # the score (to be minimized by SPOT) is "minus test score"
      result$sd.y=sd(result$R_vali);
    }
    cat1(opts,sprintf("%s  Vali    relative gain: %7.3f",ifelse(opts$TST.kind=="cv","CV ",""),y));
    cat1(opts,ifelse(opts$NRUN>1,sprintf(" +-%7.3f",sd(result$R_vali)),""));
    cat1Records(nrow(res$d_test));

    cat1(opts,sprintf("%s Vali2    relative gain (predict always with %s): %7.3f"
                     ,ifelse(opts$TST.kind=="cv","CV ",""), res$allEVAL$test2.string, mean(result$R_vali2)));
    cat1(opts,ifelse(opts$NRUN>1,sprintf(" +-%7.3f",sd(result$R_vali2)),""));
    cat1Records(nrow(res$d_test));

    if (!is.null(dset)) result$dset=dset;          # might cost a lot of memory

    result;

} # tdmClassifySummary

