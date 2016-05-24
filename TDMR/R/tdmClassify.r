# require(randomForest);  # now via direct call 'randomForest::'
# require(e1071);         # svm(), naiveBayes 

######################################################################################
# tdmClassify
# TODO: make code also applicable if d_test is a 0-row data frame.
#
#' Core classification function of TDMR. 
#'
#'      tdmClassify is called by \code{\link{tdmClassifyLoop}} and returns an object of class \code{tdmClass}. \cr
#'      It trains a model on training set \code{d_train} and evaluates it on test set \code{d_test}.
#'      If this function is used for tuning, the test set \code{d_test} plays the role of a validation set.
#'
#' Currently  d_dis is allowed to be a 0-row data frame, but d_train and d_test must have at least one record. \cr
#'
#'   @param d_train     training set
#'   @param d_test      validation set, same columns as training set
#'   @param d_dis       'disregard set', i.e. everything what is neither train nor test. The model is 
#'                      applied to all records in d_dis (needed for active learning, see ssl_methods.r)
#'   @param d_preproc   data used for preprocessing. May be NULL, if no preprocessing is done 
#'                   (opts$PRE.SFA=="none" and opts$PRE.PCA=="none"). If preprocessing is done, 
#'                   then d_preproc is usually all non-validation data.
#'   @param response.variables   name of column which carries the target variable - or - 
#'                   vector of names specifying multiple target columns
#'                   (these columns are not used during prediction, only for evaluation)
#'   @param input.variables     vector with names of input columns 
#'   @param opts     additional parameters [defaults in brackets]
#'     \describe{
#'     \item{\code{SRF.*}}{ several parameters for \code{\link{tdmModSortedRFimport}} }
#'     \item{\code{RF.*}}{ several parameters for RF (Random Forest, defaults are set, if omitted)  }
#'     \item{\code{SVM.*}}{ several parameters for SVM (Support Vector Machines, defaults are set, if omitted)}
#'     \item{\code{filename}}{ }
#'     \item{\code{data.title}}{ }
#'     \item{\code{MOD.method}}{ ["RF"] the main training method
#'                   ["RF"|"MC.RF"|"SVM"|"NB"]: use [Random forest| MetaCost-RF| SVM| Naive Bayes] for the main model}
#'     \item{\code{MOD.SEED}}{ =NULL: get a new random number seed with \code{\link{tdmRandomSeed}} (different RF trainings). \cr
#'                   =any value: set the random number seed to this value (+i) to get reproducible random
#'                   numbers. In this way, the model training part (RF, NNET, ...) gets always a fixed seed
#'                   (see also TST.SEED in \code{\link{tdmClassifyLoop}}) }
#'     \item{\code{CLASSWT}}{ class weights (NULL, if all classes should have the same weight)
#'                   (currently used only by methods RF, MC.RF and by \code{\link{tdmModSortedRFimport}})  }
#'     \item{\code{fct.postproc}}{ [NULL] name of user-def'd function for postprocessing of predicted output  }
#'     \item{\code{GD.DEVICE}}{ if !="non", then make a pairs-plot of the 5 most important variables
#'                   and make a true-false bar plot }
#'     \item{\code{VERBOSE}}{ [2] =2: most printed output, =1: less, =0: no output }
#'     }
#'         
#'   @return  \code{res}, an object of class \code{tdmClass}, this is a list containing
#'       \item{\code{d_train}}{ training set + predicted class column(s) }
#'       \item{\code{d_test}}{ test set + predicted class column(s) }
#'       \item{\code{d_dis}}{  disregard set + predicted class column(s)  }
#'       \item{\code{avgEVAL}}{ list with evaluation measures, averaged over all response variables }
#'       \item{\code{allEVAL}}{ data frame with evaluation measures, one row for each response variable }
#'       \item{\code{lastCmTrain}}{ a list with evaluation info for training set (confusion matrix, gain, class errors, ...)  }
#'       \item{\code{lastCmVali}}{  a list with evaluation info for validation set (confusion matrix, gain, class errors, ...) }
#'       \item{\code{lastModel}}{       the last model built (i.e. for the last response variable) }
#'       \item{\code{lastProbs}}{    a list with three probability matrices (row: records, col: classes) v_train, v_test, v_dis, if the model provides probabilities; NULL else. }
#'       \item{\code{lastPred}}{     name of the colum where the prediction of the last model is appended to the datasets d_train, d_test and d_dis }
#'       \item{\code{predProb}}{     a list with two data frames Trn and Val. They contain at least a column IND.dset (index of each train / validation 
#'              record into data frame dset). If the model has probabilities, then they contain in addition a column for each response
#'              variable with the prediction probabilities. }
#'       \item{\code{opts}}{ parameter list from input, some default values might have been added }
#'
#'    The 9 evaluation measures in avgEVAL and allEVAL are
#'         cerr.*    (misclassification errror),
#'         gain.*    (total gain) and
#'         rgain.*   (relative gain, i.e. total gain divided by max. achievable gain in *)
#'    where * = [trn | tst | tst2 ] stands for [ training set | test set | test set with special treatment ]
#'    and the special treatment is either opts$test2.string = "no postproc" or = "default cutoff".
#'       \cr\cr
#'    The five items \code{lastCmTrain}, \code{lastCmVali}, \code{lastModel}, \code{lastProbs}, \code{lastPred} are 
#'    specific for the *last* model (the one built for the last response variable in the last run and last fold) 
#'
#' @seealso  \code{\link{print.tdmClass}} \code{\link{tdmClassifyLoop}} \code{\link{tdmRegressLoop}}
#' @author Wolfgang Konen, FHK, Sep'2009 - Jun'2012
#' @examples
#' #*# This demo shows a simple data mining process (phase 1 of TDMR) for classification on
#' #*# dataset iris.
#' #*# The data mining process in tdmClassify calls randomForest as the prediction model.
#' #*# It is called opts$NRUN=1 time with one random train-validation set splits.
#' #*# Therefore data frame res$allEval has one row
#' #*#
#' opts=tdmOptsDefaultsSet()                       # set all defaults for data mining process
#' gdObj <- tdmGraAndLogInitialize(opts);          # init graphics and log file
#' 
#' data(iris)
#' response.variables="Species"                    # names, not data (!)
#' input.variables=setdiff(names(iris),"Species")
#' opts$NRUN=1
#' 
#' idx_train = sample(nrow(iris))[1:110]
#' d_train=iris[idx_train,]
#' d_vali=iris[-idx_train,]
#' d_dis=iris[numeric(0),]
#' res <- tdmClassify(d_train,d_vali,d_dis,NULL,response.variables,input.variables,opts)
#' 
#' cat("\n")
#' print(res$allEVAL)
#'
#' @export
######################################################################################
tdmClassify <- function(d_train,d_test,d_dis,d_preproc,response.variables,input.variables,opts)
{
    first <- TRUE; 
    filename <- opts$filename
    saved.input.variables <- input.variables;       # save copy for response.variable loop
    
    opts <- tdmOptsDefaultsSet(opts);               # be sure that all necessary  defaults are set

    if (opts$MOD.method %in% c("RF","MC.RF","SVM")) {
      d_train <- na_input_check(d_train,"d_train",input.variables);
      d_test <- na_input_check(d_test,"d_test",input.variables);
      if (nrow(d_dis)>0) d_dis  <- na_input_check(d_dis,"d_dis",input.variables);
    }

    if (is.null(opts$PRE.SFA.numericV)) opts$PRE.SFA.numericV <- input.variables;
    if (!is.null(opts$RF.mtry)) opts$RF.mtry=min(length(input.variables),opts$RF.mtry);
    # this is to avoid the "invalid mtry" warning (but it does not eliminate all cases)
    if (is.null(opts$i)) opts$i=1;               
    if (is.null(opts$k)) opts$k=1;  
    if (is.null(opts$the.nfold)) opts$the.nfold=1;
    
    if (opts$i==opts$NRUN & opts$k==opts$the.nfold) {
      if (opts$GD.RESTART) {
        tdmGraphicCloseDev(opts); 
        tdmGraphicInit(opts);
      }
    }   
    predProb <- list()
    predProb$Trn <- data.frame(IND.dset=d_train$IND.dset);
    predProb$Val <- data.frame(IND.dset=d_test$IND.dset);
    
    avgEVAL <- allEVAL <- NULL;    
    for (response.variable in response.variables) {   
        input.variables <- saved.input.variables;
        
        # SFA preprocessing, if requested by opts$PRE.SFA, is done *inside* the response.variable-for-loop
        # because SFA training depends on the response variable
        if (opts$PRE.SFA!="none") {
          # a) do SFA on the numeric variables of d_train, if opts$PRE.SFA!="none"
          # b) add monomials of degree 2 for the first opts$PRE.SFA.npc numeric variables
          # c) apply this SFA and monomials to d_test and d_dis in the same way
          other.variables <- setdiff(input.variables,opts$PRE.SFA.numericV);
          sfa <- tdmPreSFA.train(d_preproc,response.variable,opts);                 # see tdmPreprocUtils.r                 
          d_train <- tdmPreSFA.apply(d_train,sfa$sfaList,opts,d_train)$dset;
          d_test <- tdmPreSFA.apply(d_test,sfa$sfaList,opts,d_train)$dset;
          d_dis <- tdmPreSFA.apply(d_dis,sfa$sfaList,opts,d_train)$dset;

          input.variables <- union(sfa$numeric.variables,other.variables);
          if (length(setdiff(input.variables,names(d_train)))>0) 
              stop("Some elements of input.variables are not columns of d_train");
        }

        if (nrow(d_test)==0 & !(opts$MOD.method %in% c("RF","MC.RF"))) 
          warning("No validation data and opts$MOD.method is not based on RF!");
        
        # tdmParaBootstrap is now also inside response.variable for-loop, because the enhancement depends on the 
        # response variable   
        if (opts$ncopies>0) {
            if (length(response.variables)>1) stop("tdmParaBootstrap currently only allowed for tasks with single response variable --> check later TDMR versions");
            # because tdmParaBootstrap is currently for a single response variable and because the code in the two lines below
            # will not work correctly if lines are added to d_train in more than one pass 
            cat1(opts,opts$filename,": Adding", opts$ncopies, "parametric bootstap patterns to training set (tdmParaBootstrap) ...\n");
            d_train <- tdmParaBootstrap(d_train,response.variable,input.variables,opts);
            predProb$Trn <- data.frame(IND.dset=d_train$IND.dset);
        }
        if (!is.factor(d_train[,response.variable]))   {
          warning(paste("Column",response.variable," of d_train is not a factor >> we change it to factor!"));
          d_train[,response.variable] <- as.factor(d_train[,response.variable]);
        }
        lev.resp <- levels(d_train[,response.variable]);
        n.class <- length(lev.resp);
        if (!is.null(opts$CLS.CLASSWT)) {
          if (any(is.na(opts$CLS.CLASSWT))) 
            stop(sprintf("opts$CLS.CLASSWT contains NA's! Consider appropriate lines 'opts$CLS.CLASSWT[i] = ...' in APD file.\n  opts$CLS.CLASSWT = "),
                 paste(opts$CLS.CLASSWT,collapse=", "))
          
          if (n.class !=  length(opts$CLS.CLASSWT)) 
            stop("Length of opts$CLS.CLASSWT differs from the number of levels in response variable.");
          if (is.null(names(opts$CLS.CLASSWT))) 
            names(opts$CLS.CLASSWT) <- lev.resp;    # RF needs a *named* vector for option 'classwt', otherwise it will change (!) this vector 
        }
        if (length(setdiff(names(opts$CLS.CLASSWT),lev.resp))>0) 
          stop("Names in opts$CLS.CLASSWT differ from the levels of the response variable.");

        opts$CLS.cutoff <- tdmModAdjustCutoff(opts$CLS.cutoff,n.class,text="opts$CLS.cutoff");
        opts$SRF.cutoff <- opts$CLS.cutoff; 
        # --- old  and deprecated (too complicated: ---
        # opts$SRF.cutoff <- tdmModAdjustCutoff(opts$SRF.cutoff,n.class,text="opts$SRF.cutoff");
        
        if (is.null(opts$CLS.gainmat)) {
          # default (generic) gain matrix:
          opts$CLS.gainmat <- matrix ( diag(n.class),nrow=n.class,ncol=n.class,
                                  dimnames=list(lev.resp,   # row names    (true levels)
                                                lev.resp    # column names (predicted levels)
                                                ))
          # instead of diag(n.class) an existing opts$CLS.gainmat can specify a user def'd gain matrix, e.g. for n.class=2:
          #     gainmat <- matrix(c(+1,-2,
          #                         0,+1), byrow=TRUE,nrow=2,ncol=2, ...
        }
                
        #=============================================
        # PART 4.1: SUMMARY OF TRAINING DATA
        #=============================================
        #cat1(opts,filename,": Summary of training data ...\n")
        #print(summary(d_train))           # most columns are of numeric type
                                           # -> summary min,max,quantiles...
             
        #=============================================
        # PART 4.2: IMPORTANCE SELECTION (BY USING RF)
        #=============================================
        # determine the importance of all input var's by constructing a test RF
        # --- this step is skipped if SRF.kind=="none", then you use all     ---
        # --- input variables and you do not see the importance of variables ---
        if (opts$SRF.kind!="none") {
          cat1(opts,filename,": Importance check ...\n");
          opts$RF.sampsize <- tdmModAdjustSampsizeC(opts$SRF.samp, d_train, response.variable, opts);
          SRF <- tdmModSortedRFimport(d_train,response.variable,
                                      input.variables,opts)

          input.variables <- as.character(SRF$input.variables);
          opts <- SRF$opts;       # some defaults might have been added, some opts$SRF.* values or list opts$srf might be changed
          SRF$opts <- NULL;       # simplify list result, which will contain both, SRF and opts
        }  else {
          SRF=NULL;
          if (opts$i==1) {
            cat1(opts,filename,": Using all input variables: \n");
            if (opts$VERBOSE>=1) print(input.variables);
          }
        } 
        
        # We set here the random number generator (RNG) seed again such that the subsequent RF training 
        # starts from the same seed, regardless whether opts$SRF.kind=="none" or !="none" (the latter
        # means extra calls to RNG in tdmModSortedRFimport)
        if (is.null(opts$MOD.SEED)) {
          # NEW: when called via SPOT, the RNG might be at (different but) fixed seed in each call.
          #      But if MOD.SEED==NULL we want different seeds (for RF training) to see the variability       
          set.seed(tdmRandomSeed());                                                                  
        } else if (opts$MOD.SEED=="algSeed") {  # use the seed from SPOT:
          # opts$ALG.SEED is set in tdmStartSpot to des$SEED[k]+r. This meens that the overall r'th 
          # evaluation of a design point gets the seed spotConfig$alg.seed+r
          newseed=opts$ALG.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed); 
        } else {
          newseed=opts$MOD.SEED+(opts$i-1)+opts$NRUN*(opts$rep-1);
          set.seed(newseed) # if you want reproducably the same model training,
        }                   # but different for each run i

    
        res.rf <- res.nb <- train.predict <- test.predict <- NULL
        #================================================================
        # PART 4.3: MODELING: TRAIN RANDOM FOREST (OR OTHER METHOD)
        #================================================================
        to.train <- d_train[,c(input.variables,response.variable)]
        to.test <- d_test[,c(input.variables,response.variable)]
        if (opts$MOD.method %in% c("RF","MC.RF"))
          opts$RF.sampsize <- tdmModAdjustSampsizeC(opts$RF.samp, to.train, response.variable, opts);
        train.rf <- function(response.variable,to.train,opts) {
            cat1(opts,opts$filename,": Train RF with sampsize =", opts$RF.sampsize,"...\n")
            testit::assert("Cutoff is bigger than 1",sum(opts$CLS.cutoff)<=1)
            if (!is.null(opts$CLS.CLASSWT)) {
                cat1(opts,"Class weights: ", opts$CLS.CLASSWT,"\n")
                #cwt = opts$CLS.CLASSWT*1.0;  # strange, but necessary: if we omit '*1' then cwt seems to be a copy-by-reference of
                                              # opts$CLS.CLASSWT. After randomForest call  cwt is changed and also opts$CLS.CLASSWT would be changed (!)
                cwt = sprintf("classwt=c(%s)",paste(opts$CLS.CLASSWT,collapse=","))
            } else { cwt="classwt=NULL"; }
            if (!is.null(opts$CLS.cutoff)) cat1(opts,"Cutoff: ", opts$CLS.cutoff,"\n")
            formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
            # we work here with a command text string and eval(parse(...)) to allow for the presence or
            # absence of certain options like "mtry" or "cutoff" which are not allowed to be NULL.
            # The individual "eval(...)" on the following lines are for clarity of res.rf$call 
            # (it should read "..., ntree=400, ..." and not "..., ntree=opts$RF.ntree, ...") 
            # BUT: we cannot use "eval(...)" for the lists cutoff and cwt (and sampsize), therefore we use here
            # the 'non-speaking' variables cwt, opts$CLS.cutoff and add below res.rf$cutoff, res.rf$classwt 
            # for optional later reference or user inspection.
            #
            # TODO: a simpler alternative could be: 
            #     rf.options = paste(rf.options," sampsize=c(",paste(opts$RF.samp,collapse=","),"),",sep="")
            # which writes " sampsize=c(200,300),". TODO: CHECK Validity of this approach
            #
            rf.options = paste(" ntree=",eval(opts$RF.ntree));
            rf.options = paste(rf.options," sampsize=opts$RF.sampsize",sep=",")
            rf.options = paste(rf.options,cwt," na.action=randomForest::na.roughfix"," proximity=FALSE",sep=",")
            #if (!is.null(cwt))  paste(rf.options,paste("classwt=",eval(cwt)),sep=",")    # not run
            if (!is.null(opts$RF.mtry)) rf.options = paste(rf.options,paste(" mtry=",eval(opts$RF.mtry)),sep=",")
            if (!is.null(opts$CLS.cutoff)) rf.options = paste(rf.options," cutoff=opts$CLS.cutoff",sep=",")
            if (!is.null(opts$RF.nodesize)) rf.options = paste(rf.options,paste(" nodesize=",eval(opts$RF.nodesize)),sep=",")
            #dbg_chase_cutoff_bug(formul,to.train,d_train,response.variable,rf.options,opts);
            flush.console();       
            eval(parse(text=paste("res.rf <- randomForest::randomForest( formul, data=to.train,",rf.options,")"))); 
            res.rf$HasVotes = TRUE; 
            res.rf$HasProbs = TRUE; 
            res.rf$cutoff = opts$CLS.cutoff;
            res.rf$classwt = opts$CLS.CLASSWT;
            res.rf$sampsize = opts$RF.sampsize;
            res.rf;
        } 
        train.mc.rf <- function(response.variable,to.train,opts) {
            cat1(opts,opts$filename,": Train tdmMetacostRf ...\n")
            flush.console();
            res.rf <- tdmMetacostRf(response.variable,to.train,opts)
            res.rf$HasVotes = TRUE; 
            res.rf$HasProbs = TRUE; 
            res.rf;
        }
#        train.kSVM <- function(response.variable,to.train,opts) {
#            cat1(opts,filename,": Train kSVM ...\n")
#		        res.rf = list()
#      		  class(res.rf) <- "TDMres_c"
#            formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
#  	        flush.console();
#            res.rf$model <- ksvm(formul
#                    	           , data=to.train
#          						           , kernel="rbfdot"
#	                               , kpar=list(sigma=opts$SVM.sigma)
#                        	       , cost=opts$SVM.cost
#                         	       , type="spoc-svc"  ### type=spoc-svc   #type="C-svc"
#            							       , class.weights=opts$CLS.CLASSWT 
#            							       , na.action=na.omit)
#            res.rf$HasVotes = FALSE;            
#            res.rf$HasProbs = FALSE; 
#	        res.rf;
#        } # train.kSVM
        train.svm <- function(response.variable,to.train,opts) {
            kernelChoices = c("linear","polynomial","radial","sigmoid");
            kernelType = kernelChoices[opts$SVM.kernel];
            cat1(opts,filename,": Train SVM (kernel=",kernelType,") ...\n");
	          #require(e1071)		
            if (!is.null(opts$CLS.CLASSWT)) cat1(opts,"Class weights: ", opts$CLS.CLASSWT,"\n")
            if (!is.null(opts$CLS.cutoff)) cat1(opts,"Cutoff: ", opts$CLS.cutoff,"\n")
            flush.console();
            formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
            res.rf <- e1071::svm(formul 
            							, data=to.train 
            							, kernel=kernelType
                          , gamma=opts$SVM.gamma
                          , coef0=opts$SVM.coef0
                          , degree=opts$SVM.degree
                 	        , cost=opts$SVM.cost
                          , type="C-classification"
            							, class.weights = opts$CLS.CLASSWT
            							, tolerance = opts$SVM.tolerance
            							, na.action=na.omit
            							, probability=TRUE)
            res.rf$HasVotes = FALSE;            
            res.rf$HasProbs = TRUE; 
      			res.rf;
        }# train.svm
        train.ada <- function(response.variable,to.train,opts) {
            coefChoices = c("Breiman","Freund","Zhu");
            coefType = coefChoices[opts$ADA.coeflearn];
            cat1(opts,filename,": Train AdaBoost ( mfinal =",opts$ADA.mfinal,", coeflearn =",coefType,") ...\n");
            flush.console();
            formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
            
            
            #require(adabag)
            if (requireNamespace("adabag", quietly = TRUE)) {
              res.rf <- adabag::boosting(formul
                                         , data=to.train
                                         , mfinal=opts$ADA.mfinal
                                         , coeflearn=coefType
              );
            } else {
              stop("Package adabag is not available on this platform")
            }
            
            res.rf$HasVotes = FALSE;
            res.rf$HasProbs = TRUE;
      			res.rf;
        }# train.ada
        train.naiveBayes <- function(response.variable,to.train,opts) {
            cat1(opts,filename,": Train NB ...\n")
            flush.console();
            formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
            # CAUTION 1: Be careful, that the column ordering in 'naiveBayes(...,data=..)' is 
            # the same as below in 'predict(res.nb,newdata=..)'. Therefore we need here
            #     to.train=d_train[,c(input.variables,response.variable)]
            # instead of a simple 'data=d_train'
            res.rf <- e1071::naiveBayes(formula=formul, data=to.train)
            res.rf$HasVotes = FALSE;                          
            res.rf$HasProbs = FALSE; 
            res.rf;
        }
        train.CMB <- function(opts) {
            res.rf <- list();
            res.rf$HasVotes = FALSE;                          
            res.rf$HasProbs = TRUE;
            res.rf;
        }

        ptm <- proc.time()
        cat1(opts, "Run ",ifelse(opts$the.nfold>1,paste(opts$i,".",opts$k,sep="")           ,opts$i)    ,"/",
                          ifelse(opts$the.nfold>1,paste(opts$NRUN,".",opts$the.nfold,sep=""),opts$NRUN) ,":\n"); 
        res.rf = switch(opts$MOD.method
          ,"RF"  =  train.rf(response.variable,to.train,opts)
          ,"MC.RF"  =  train.mc.rf(response.variable,to.train,opts)
          ,"SVM" =  train.svm(response.variable,to.train,opts)
          ,"ADA" =  train.ada(response.variable,to.train,opts)
#          ,"kSVM" =  train.kSVM(response.variable,to.train,opts)
          ,"NB"  =  train.naiveBayes(response.variable,to.train,opts)
          ,"CMB" = train.CMB(opts)
          ,"INVALID"
          );

        if (res.rf[1]=="INVALID") {
          cat1(opts,sprintf("*** Invalid opts$MOD.method=%s ***\n",opts$MOD.method));
        }
        #print(res.rf)
        cat1(opts,"Proc time: ",(proc.time()-ptm)[1],"\n");       
    
        #=====================================================================
        # PART 4.4: APPLY RANDOM FOREST  (OR OTHER METHOD)
        #=====================================================================
        apply.rf <- function(res.rf,to.train,to.test,to.dis,opts) {        
            cat1(opts,filename,": Apply",opts$MOD.method,"...\n")
            app = list()            
            app$train.predict <- res.rf$predicted
            # (res.rf$predicted is the *OOB-prediction* on the training set)
            # Think about this! Why is it WRONG (or too optimistic) to use here
            #      app$train.predict <- predict(res.rf, newdata=d_train)   
            # as the prediction for the training set?
            
            app$test.predict <- predict(res.rf, newdata=to.test)
            app$test.votes <- predict(res.rf, newdata=to.test, type="vote")
            #app$test.predict <- factor(rep(0,nrow(to.test)),levels=c(0,1)) # dummy for DMC2010 to predict the naive model
            
            app$dis.predict <- predict(res.rf, newdata=to.dis)
            app$dis.votes <- predict(res.rf, newdata=to.dis, type="vote")
            
            app;
        }
        apply.other <- function(res.rf,to.train,to.test,to.dis,opts) {
            predict.TDMres_c <- function(res.rf, newdata) {predict(res.rf$model, newdata=newdata);}#required for S4 models
            cat1(opts,opts$filename,": Apply",opts$MOD.method,"...\n")
            app = list()  
            app$train.predict <- predict(res.rf, newdata=to.train)
            app$test.predict <- predict(res.rf, newdata=to.test)
            if (nrow(to.dis)>0) app$dis.predict <- predict(res.rf, newdata=to.dis)
            app;
        }
        apply.svm <- function(res.rf,to.train,to.test,to.dis,opts) {
            cat1(opts,opts$filename,": Apply",opts$MOD.method,"...\n")
            app = list()  

            app$train.prob <- attr(predict(res.rf, newdata=to.train, probability=TRUE),"probabilities");
            if (is.null(opts$CLS.cutoff)) {
              app$train.predict <- predict(res.rf, newdata=to.train)
            } else {
              cat1(opts,"Cutoff: ", opts$CLS.cutoff,"\n")
  			      # Apply new weighting scheme for training set            
     			    app$train.predict <- apply(app$train.prob/(matrix(1,nrow(to.train),1) %*% opts$CLS.cutoff),1,which.max);
     			    app$train.predict <- colnames(app$train.prob)[app$train.predict];
     			    # -- only for debug: --
   	  		    #if (!check.apply.svm(res.rf,app$train.prob,app$train.predict,to.train,d_train[,response.variable],opts)) browser();  
            }

      			if (nrow(to.test)>0) {
              app$test.prob <- attr(predict(res.rf, newdata=to.test, probability=TRUE),"probabilities");
              if (is.null(opts$CLS.cutoff)) {
                app$test.predict <- predict(res.rf, newdata=to.test)
              } else {
          			# Apply new weighting scheme for test set            
       			    app$test.predict <- apply(app$test.prob/(matrix(1,nrow(to.test),1) %*% opts$CLS.cutoff),1,which.max)
     	  		    app$test.predict <- colnames(app$test.prob)[app$test.predict];
       			    # -- only for debug: --
     	  		    #if (!check.apply.svm(res.rf,app$test.prob,app$test.predict,to.test,d_test[,response.variable],opts)) browser(); 
    			    }
  			    }
  			    
      			if (nrow(to.dis)>0) {
        			app$dis.prob <- attr(predict(res.rf, newdata=to.dis, probability=TRUE),"probabilities");
              if (is.null(opts$CLS.cutoff)) {
                app$dis.predict <- predict(res.rf, newdata=to.dis)
              } else {
          			# Apply new weighting scheme for disregard set (unknown label data, e.g. for active learning)
       			    app$dis.predict <- apply(app$dis.prob/(matrix(1,nrow(to.dis),1) %*% opts$CLS.cutoff),1,which.max)
       			    app$dis.predict <- colnames(app$dis.prob)[app$dis.predict];
    			    }
  			    }
            app;
        }
        apply.ada <- function(res.rf,to.train,to.test,to.dis,opts) {        
            cat1(opts,filename,": Apply",opts$MOD.method,"...\n")
            app = list()                        
            pa <- predict(res.rf, newdata=to.train);
            app$train.predict <- pa$class;
            app$train.votes <- pa$votes;
            app$train.prob <- pa$prob;
            
            pa <- predict(res.rf, newdata=to.test);
            app$test.predict <- pa$class;
            app$test.votes <- pa$votes;
            app$test.prob <- pa$prob;
            
      			if (nrow(to.dis)>0) {
              pa <- predict(res.rf, newdata=to.dis);
              app$dis.predict <- pa$class;
              app$dis.votes <- pa$votes;
              app$dis.prob <- pa$prob;
            }
            
            app;
        }
        apply.CMB <- function(to.train,to.test,to.dis,opts) {
            cat1(opts,opts$filename,": Apply",opts$MOD.method,"...\n")
            app = list()  
  
            app$train.prob <- opts$result1$lastRes$lastProbs$v_train * opts$result2$lastRes$lastProbs$v_train;
      			app$test.prob <- opts$result1$lastRes$lastProbs$v_test * opts$result2$lastRes$lastProbs$v_test;
      			app$dis.prob <- opts$result1$lastRes$lastProbs$v_dis * opts$result2$lastRes$lastProbs$v_dis;

   			    app$train.predict <- colnames(app$train.prob)[apply(app$train.prob,1,which.max)];
   			    app$test.predict <- colnames(app$test.prob)[apply(app$test.prob,1,which.max)];
   			    app$dis.predict <- colnames(app$dis.prob)[apply(app$dis.prob,1,which.max)];   			    
   			    app;
        }
          
        ptm <- proc.time()
        opts$response.variable <- response.variable;
        pred.vars = input.variables;
        if (opts$MOD.method=="ADA") pred.vars = c(pred.vars,response.variable);
        # CAUTION 2: AdaBoost requires both, response.variable and input.variables
        
        #if (opts$MOD.method=="NB") {
          # CAUTION 3: predict.naiveBayes requires that only the input variables
          # enter via 'newdata', so we have to exclude column response.variable
          # Otherwise a strange error message occurs:
          #     Fehler in FUN(1:6[[1L]], ...) : Indizierung ausserhalb der Grenzen
          to.train <- data.frame(a=d_train[,pred.vars]); 
          to.test <- data.frame(a=d_test[,pred.vars]);
          to.dis <- data.frame(a=d_dis[,pred.vars]);
          names(to.train) <- names(to.test) <- names(to.dis) <- pred.vars;
          # Why 'data.frame', 'a=' and 'names...'? - In this way it works also
          # for length(input.variables)==1. If we had instead the simple 
          #     to.train <- d_test[,input.variables]
          # this would lead to an unnamed vector.
        #}

        app = switch(opts$MOD.method
          ,"RF" =,"MC.RF" =  apply.rf(res.rf,to.train,to.test,to.dis,opts)
          ,"NB" =, "kSVM" = apply.other(res.rf,to.train,to.test,to.dis,opts) 
    		  ,"SVM" = apply.svm(res.rf,to.train,to.test,to.dis,opts)
    		  ,"ADA" = apply.ada(res.rf,to.train,to.test,to.dis,opts)
    		  ,"CMB" = apply.CMB(to.train,to.test,to.dis,opts)
          ,"INVALID"
          );
        if (app[1]=="INVALID") {
          cat1(opts,sprintf("*** Invalid opts$MOD.method=%s ***\n",opts$MOD.method));
        }
        train.predict <- app$train.predict;
        test.predict <- app$test.predict;
        dis.predict <- app$dis.predict;

        #
        # column "votes" is beneficial for opts$fct.postproc (e.g. in the case of DMC2010)
        # (it is currently only implemented for opts$MOD.method=="RF", so that is why the other methods have res.rf$HasVotes==FALSE)
        #
        if (res.rf$HasVotes) {  # column "votes" is the fraction of trees voting for the majority class
          name.of.votes <- "votes";
          d_train <- tdmBindResponse(d_train, name.of.votes, res.rf$votes[,1])
          d_test  <- tdmBindResponse(d_test, name.of.votes, app$test.votes[,1])
          d_dis   <- tdmBindResponse(d_dis, name.of.votes, app$dis.votes[,1])
          if (opts$test2.string == "default cutoff") {
            #-- choice 1: test2 is the prediction with default cutoff -------------------------------
            cutoff <- rep(1/n.class, n.class);
            test2.predict <- predict(res.rf, newdata=to.test, cutoff=cutoff)
            test2.votes <- predict(res.rf, newdata=to.test, cutoff=cutoff, type="vote")
          } else { # i.e. opts$test2.string == "no postproc";
            #-- choice 2: test2 is a mere copy of test for which we later skip the postprocessing --- 
            test2.predict <- test.predict;
            test2.votes <- app$test.votes[,1]; 
          }
        }
        if (res.rf$HasProbs) {
          lastProbs = NULL;                                        # lastProbs is used by ssl_methods.r
          if (opts$MOD.method %in% c("RF","MC.RF")) {
            lastProbs = list(  v_train = res.rf$votes                
                              , v_test = app$test.votes
                              , v_dis = app$dis.votes
                              );
          }
          if (opts$MOD.method %in% c("SVM","CMB","ADA")) {
            lastProbs = list(  v_train = app$train.prob             
                              , v_test = app$test.prob
                              , v_dis = app$dis.prob
                              );
          }
        } else {
          lastProbs = NULL;
        } 
        cat1(opts,"Proc time: ",(proc.time()-ptm)[1],"\n");

    ######################################################################################
    # this part is specific for the DMC2010-task and method = "RF" (only needed for human analysis)
    ######################################################################################   
    # just for comparision: with default cutoff (and w/o postprocessing) test2.predict should  
    # be strictly worse in comparision to test.predict:
    VOTE2TARGET.CHK = FALSE
    if ( VOTE2TARGET.CHK ) {
      vres.train <- tdmModVote2Target(res.rf$votes[,1],res.rf$predicted,d_train[,response.variable]);
      print(vres.train);
      #vres.test <- tdmModVote2Target(app$test.votes[,1],test.predict,d_test[,response.variable]);
      #print(vres.test);
      #vres.test2 <- tdmModVote2Target(test2.votes[,1],test2.predict,d_test[,response.variable]);
      #print(vres.test2);
    }
    ######################################################################################
                  
            if (opts$APPLY_TIME) {
              xstart = proc.time()
              for (t in 1:100) {
                dummy <- predict(res.rf, newdata=to.train)
              }
              xend = proc.time()
              cat1(opts,"\nElapsed time for 100x APPLY on train set:\n")
              print(xend-xstart)
              print(nrow(to.train))
            }
    
        #=============================================
        # PART 4.5: POSTPROCESSING
        #=============================================
      	if (!is.null(opts$fct.postproc)) {
      		cat1(opts,filename,": User-defined postprocessing: Applying function",opts$fct.postproc," ...\n") 		
      		train.predict <- eval(parse(text=paste(opts$fct.postproc,"(train.predict,d_train,opts)",sep="")));
      		test.predict <- eval(parse(text=paste(opts$fct.postproc,"(test.predict,d_test,opts)",sep="")));
      		dis.predict <- eval(parse(text=paste(opts$fct.postproc,"(dis.predict,d_dis,opts)",sep="")));
      		if (res.rf$HasVotes) {
            opts2 <- opts; 
            opts2$fct.postproc=NULL; # i.e. no postprocessing for test2
            if (!is.null(opts2$fct.postproc))                                   
          		 test2.predict <- eval(parse(text=paste(opts$fct.postproc,"(test2.predict,d_test,opts)",sep="")));
      		}
      	}
        # bind the predicted class pred_... as last column to the data frames
        name.of.prediction <- paste("pred_", response.variable, sep="")
        name2.of.prediction <- paste("pred2_", response.variable, sep="")
        name.of.probability <- paste("prob_", response.variable, sep="")
        d_train <- tdmBindResponse(d_train, name.of.prediction, train.predict)
        d_test <- tdmBindResponse(d_test, name.of.prediction, test.predict)
        d_dis <- tdmBindResponse(d_dis, name.of.prediction, dis.predict)
        if (res.rf$HasVotes) d_test  <- tdmBindResponse(d_test, name2.of.prediction, test2.predict)
        if (res.rf$HasProbs & !is.null(d_train$IND.dset) & !is.null(d_test$IND.dset)) {
          predProb$Trn <- tdmBindResponse(predProb$Trn, response.variable, d_train[,response.variable]);
          predProb$Trn <- tdmBindResponse(predProb$Trn, name.of.prediction, train.predict)
          predProb$Trn <- tdmBindResponse(predProb$Trn, name.of.probability, as.vector(lastProbs$v_train[,1]))
          predProb$Val <- tdmBindResponse(predProb$Val, response.variable, d_test[,response.variable]);
          predProb$Val <- tdmBindResponse(predProb$Val, name.of.prediction, test.predict)
          predProb$Val <- tdmBindResponse(predProb$Val, name.of.probability, as.vector(lastProbs$v_test[,1]))
        } 

    
        #=============================================
        # PART 4.6: EVAL: CALC CONFUSION MATRIX + GAIN
        #=============================================
        cat1(opts,filename,": Calc confusion matrix + gain ...\n")
        
        cat1(opts,"\nTraining cases (",length(train.predict),"):\n")
        cm.train <- tdmModConfmat(d_train,response.variable,name.of.prediction,opts,predProb$Trn);
        # the contents of cm.train$rgain depends on opts$rgain.type
        
        #cat1(opts,"   --- predicted ---\n")
        if (opts$VERBOSE>=1) print(cm.train$mat)                     # confusion matrix on training set
        cat1(opts,sprintf("total gain: %7.1f (is %7.3f%% of max. gain = %7.1f)\n", 
                          cm.train$gain,cm.train$gain/cm.train$gainmax*100,cm.train$gainmax));

        if (nrow(d_test)>0) {
          cat1(opts,"\nValidation cases (",length(test.predict),"):\n")
          cm.vali <- tdmModConfmat(d_test,response.variable,name.of.prediction,opts,predProb$Val);
          print1(opts,cm.vali$mat)                      # confusion matrix on test set
          print1(opts,cm.vali$gain.vector)
          cat1(opts,sprintf("total gain : %7.1f (is %7.3f%% of max. gain = %7.1f)\n", 
                            cm.vali$gain,cm.vali$gain/cm.vali$gainmax*100,cm.vali$gainmax));
          if (res.rf$HasVotes && opts$test2.show) {                          
            cat1(opts,"\nVali2 cases [",opts$test2.string,"] (",length(test.predict),"):\n")
            cm.vali2 <- tdmModConfmat(d_test,response.variable,name2.of.prediction,opts,predProb$Val);
            print1(opts,cm.vali2$mat)                      # confusion matrix on test set, test2-prediction
            cat1(opts,sprintf("total gain2: %7.1f (is %7.3f%% of max. gain = %7.1f)\n", 
                              cm.vali2$gain,cm.vali2$gain/cm.vali2$gainmax*100,cm.vali2$gainmax));
            # DMC2010:                          
            # we usually get (w/o postprocessing): total.gain: 27.4%, total.gain2: 21.5%   
            # but with postprocessing only marginal diff: total.gain: 28.0% +- 0.3%, total.gain2: 27.7% +- 0.03% 
          }  else {
            cm.vali2 = cm.vali;
          }
        } else {    # i.e. if nrow(d_test)==0
          cm.vali <- cm.train;
          cm.vali$cerr[,] <- cm.vali$mat[,] <- cm.vali$ccase[,] <- cm.vali$gain <- NA;
          cm.vali$gain.vector[,] <- cm.vali$gainmax <- cm.vali$rgain <- NA;
          cm.vali2 <- cm.vali;
          #cm.vali <- cm.vali2 <- NULL;   # OLD, this had problems when constructing as.data.frame(EVALa) below
        }
        
        EVAL =  list(        # a list of eval-quantities which can be summed
                     cerr.trn=cm.train$cerr[1,"Total"]       # misclassification error 
                    ,gain.trn=cm.train$gain
                    ,rgain.trn=cm.train$rgain
                    ,cerr.tst=cm.vali$cerr[1,"Total"]
                    ,gain.tst=cm.vali$gain
                    ,rgain.tst=cm.vali$rgain
                    ,cerr.tst2=cm.vali2$cerr[1,"Total"]      #
                    ,gain.tst2=cm.vali2$gain                 # for comparision in method="RF": results with
                    ,rgain.tst2=cm.vali2$rgain               # default cutoff 1/n.class or with no postproc (see TEST2.MODE)
                    #,gainmax=cm.train$gainmax
                    );
        EVALa = EVAL;
        EVALa$test2.string=opts$test2.string;
        EVALa$resp = response.variable;
        EVALa$i = opts$i;
        EVALa$k = opts$k;         

        if (first) {
            avgEVAL=lapply(EVAL,sum);       # sum here, but with "/length(response.variables)" below effectively a mean
            allEVAL=as.data.frame(EVALa);
            first=FALSE;
        } else {
            avgEVAL<- as.list(mapply(sum,EVAL,avgEVAL));
            allEVAL <- rbind(allEVAL,as.data.frame(EVALa));
        }
        rownames(allEVAL)[nrow(allEVAL)] <- response.variable;
                    
        #=============================================
        # PART 4.7: GRAPHICS
        #=============================================
        if (!(opts$i==opts$NRUN & opts$k==opts$the.nfold)) {
          if (opts$GD.RESTART) {
            tdmGraphicCloseDev(opts); 
            tdmGraphicInit(opts);
          } 
        } 
        if (opts$GD.DEVICE!="non" & length(input.variables)>1) {
            tdmGraphicNewWin(opts);
            opts$gr.points=500;
            p=sample(nrow(d_train));
            ind <- p[1:opts$gr.points];
            # scatter plot for each  pair of the 5 most important input variables, colored with class levels
            # (this can become very large, if #records is high)
            mycolors <-  sample(colors(),n.class,replace=FALSE)    # select n.class random colors
            important.variables=input.variables
            maxi = min(5,length(important.variables));            
            pairs(d_train[ind,important.variables[1:maxi]], cex=0.6, gap=0,
                  col=mycolors[as.numeric(d_train[ind,response.variable])],
                  main=paste(opts$data.title,": The five most important inputs", sep=""))
            tdmGraphicCloseWin(opts);
        }
        if (opts$GD.DEVICE!="non") {
            # bar plot of true/false test cases for all classes
            tdmGraphicNewWin(opts);
            if (nrow(d_test)>0) {
              cmt <- cm.vali$mat;
              setStr = "validation set";
            } else {
              cmt <- cm.train$mat;
              setStr = "train set";
              cat1(opts,"\nNOTE: No validation cases -> we plot the true/false training cases.\n");
            }
            height<-t(matrix(c(diag(cmt),rowSums(cmt)-diag(cmt)),nrow=n.class,ncol=2))
            # 'height' is a 2*n.class matrix containing in its first row the number of
            # correct classifications for each class level (diagonal of confusion matrix)
            # and in its 2nd row the number of wrong classifications for each class level
            # (sum of off-diagonal elements in each row of the confusion matrix)
            barplot(height,beside=TRUE,col=c("blue","red"),legend=c("true","false"),
                    names.arg=colnames(cmt),
                    main=paste("True/false classification on",setStr));
            tdmGraphicCloseWin(opts);
        } # if (opts$GD.DEVICE!="non")
        
    } # for (response.variable)
    avgEVAL <- lapply(avgEVAL, function(x) {x/length(response.variables)} );
    
    #=============================================
    # PART 4.8: WRITE RESULTS ON TEST SET TO FILE
    #=============================================
    dir.output <- paste(dirname(opts$dir.output),basename(opts$dir.output),sep="/")  # remove trailing "/", if it exists
    if (!file.exists(dir.output)) {
      success = dir.create(dir.output);     
      if (!success) stop(sprintf("Could not create dir.output=%s",dir.output));
    }
    #outfile = paste(opts$dir.output,sub(".csv", "", filename), "_predictions.csv", sep="")
    #write.table(d_test, file=outfile, quote=FALSE, sep=";", dec=".", row.names=FALSE, col.names=TRUE)

    if (opts$fileMode) {
      if (!is.null(opts$EVALFILE)) {
        colNames = FALSE
        if (opts$i==1 & opts$k==1) colNames = TRUE
      	write.table(data.frame(allEVAL)
      			, file = paste(opts$dir.output,opts$EVALFILE, sep="")
      			, col.names= colNames
      			, row.names= FALSE
      			, append = !colNames         
      			, sep = ";",
      			, quote = FALSE
      			, eol = "\n"
      	);
     	}
   	}

    res =   list(lastCmTrain=cm.train   # from last response.variable 
                ,lastCmVali=cm.vali     #   "     "     "       "
                ,lastModel = res.rf     #   "     "     "       "
                ,lastPred = name.of.prediction #   "     "       "
                ,lastProbs = lastProbs  # NULL or list with 3 probability matrices (row:records, col: classes) v_train, v_test, v_dis
                ,predProb = predProb    # $Trn: data frame (IND.dset, rv1, pred_rv1, prob_rv1, rv2, ...) for training data
                                        # $Val: data frame (IND.dset, rv1, pred_rv1, prob_rv1, rv2, ...) for validation data
                #,avgEVAL=avgEVAL       # ---deprecated---
                ,allEVAL=allEVAL
                ,d_train=d_train
                ,d_test=d_test 
                ,d_dis=d_dis 
                ,SRF=SRF                  # output from tdmModSortedRFimport or NULL
                ,opts=opts                # some defaults might have been added, list opts$srf might have been added
                );
               
    class(res) <- c("tdmClass","TDM")     # this causes > res; or > print(res);
                                          # NOT to print out the whole list (might be very long!!)
                                          # but to call instead the function  print.tdmClass
    res;
}  # tdmClassify

na_input_check <- function(dfr, name_dfr, input.variables) {
    unknownCols = setdiff(input.variables,names(dfr))   
    if (length(unknownCols)>0) stop("Unknown columns %s in data set %s",paste(unknownCols,collapse=","),name_dfr);
    if (any(is.na(dfr[,input.variables]))) {
      warning(sprintf("There are NA's in the input variable columns of %s. Now replacing them using tdmPreNAroughfix().
        Consider to use this (or a more sophisticated imputation technique) prior to calling tdmClassify.",name_dfr));
      dfr[,input.variables] <- tdmPreNAroughfix(dfr[,input.variables]);      # see tdmPreprocUtils.r
    }
    dfr;
}

# debug routine for apply.svm: 
#   there are rare cases (currcntly detected in benchmark dataset wdbc) where the prediction without probabilites
#   and the prediction with probabilities in case opts$CLS.cutoff=c(0.5,0.5) (disable CUTOFF line in wdbc_02.roi) give different results, 
#   although in theory they shouldn't.
#   The differences occur in those cases, where the probablities are close to 0.5, normally only the one or two cases closest to 0.5.
#   They occur in training and in test cases.
#   The routine check.apply.svm detects them, prints some diagnostics and stops in a browser. 
#   The reason for the differences is still unclear.
#
check.apply.svm <- function(res.rf,appProb,appPredict,to.data,respVar,opts) {
  P.alt <- predict(res.rf, newdata=to.data);
  y=(P.alt==appPredict)
  retVal = TRUE;
  if (any(y==FALSE)) {
    if (length(P.alt) != nrow(appProb)) { cat("NOTE: data frame binding problems 1 ahead!\n"); browser(); }
    if (length(P.alt) != length(respVar)) { cat("NOTE: data frame binding problems 2 ahead!\n"); browser(); }
    x=data.frame(y,mid.diff=abs(0.5-appProb[,1]),P.cl=P.alt,true=respVar,M.prob=appProb[,1]);
    print(x[order(x[,2]),]);
    #browser();     	  
    #retVal = FALSE;
  }		    
  retVal;
}
