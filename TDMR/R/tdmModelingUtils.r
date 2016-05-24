######################################################################################
######################################################################################
#
# HELPER FUNCTIONS FOR MODELING
#
######################################################################################
# require(randomForest);    # now via direct call 'randomForest::'
######################################################################################
#
#
######################################################################################
# tdmModCreateCVindex:
#'   Create and return a training-validation-set index vector. 
#' 
#'   Depending on the value of member TST.kind in list opts, the returned index cvi is
#'   \enumerate{
#'   \item TST.kind="cv": a random cross validation index P([111...222...333...]) - or -
#'   \item TST.kind="rand": a random index with P([00...11...-1-1...]) for training (0), validation (1) and disregard (-1) cases - or -
#'   \item TST.kind="col": the column dset[,opts$TST.COL] contains the training (0), validation (1) and disregard (-1) set division
#'         (and all records with a value <0 in column TST.COL are disregarded).
#'   }
#'   Here P(.) denotes random permutation of the sequence. \cr
#'   The disregard set is optional, i.e. cvi may contain only 0 and 1, if desired. \cr
#'   Special case TST.kind="cv" and TST.NFOLD=1: make *every* record a training record, i.e. index [000...]. \cr
#'   In case TST.kind="rand" and stratified=TRUE a \emph{stratified} sample is drawn, where the strata in the 
#'   training case reflect the rel. frequency of each level of the **1st** response variable
#'   and are ensured to be at least of size 1.
#'
#'  @param dset    the data frame for which cvi is needed
#'  @param response.variables  issue a warning if \code{length(response.variables)>1}. Use the first
#'      response variable for determining strata size.
#'  @param opts    a list from which we need here the following entries
#'    \itemize{
#'      \item TST.kind:  ["cv"|"rand"|"col"]
#'      \item TST.NFOLD: number of CV folds (only relevant in case TST.kind=="cv")
#'      \item TST.COL:   column of dset containing the (0/1/<0) index (only relevant in case TST.kind=="col")
#'                       or NULL if no such column exists
#'      \item TST.valiFrac:  fraction of records to set aside for validation (only relevant in case TST.kind=="rand")
#'      \item TST.trnFrac:  [1-opts$TST.valiFrac] fraction of records to use for training (only relevant in case TST.kind=="rand")
#'    }
#'  @param stratified [F] do stratified sampling for TST.kind="rand" with at least one training 
#'         record for each response variable level (classification)
#'
#'  @return cvi  training-validation-set (0/>0) index vector
#'               (all records with cvi<0, e.g. from column TST.COL, are disregarded)
#'
#' @note Currently stratified sampling in case TST.KIND='rand' does only work correctly for \emph{one} response variable.  
#'    If there are more than one, the right fraction of validation records is taken, but the strata are drawn w.r.t. the 
#'    first response variable. (For multiple response variables we would have to return a list of cvi's or to
#'    call tdmModCreateCVindex for each response variable anew.)
#
######################################################################################
tdmModCreateCVindex <- function(dset,response.variables,opts,stratified=FALSE) {
        L = nrow(dset)
        if (opts$TST.kind=="cv") {
            #===============================================
            # CREATE NFOLD CROSSVALIDATION INDEX
            #===============================================
            cat1(opts,opts$filename,": Creating cross validation index cvi ...\n")
            nfold = opts$TST.NFOLD
            if (opts$TST.NFOLD <= 0)
              stop(sprintf("tdmModCreateCVindex: opts$TST.NFOLD must be > 0. Current value is opts$TST.NFOLD = %d", opts$TST.NFOLD));
            cvi = NULL
            ilow=1
            for (m in 1:nfold) {
                ihigh = floor(L/nfold*m)
                cvi[ilow:ihigh] = m         # i.e. [111..2222...333...]
                ilow = ihigh+1
            }
            cvi = sample(cvi)               # permute CVI randomly
            if (opts$TST.NFOLD == 1)
              cvi = 0*cvi;                  # special case: if NFOLD==1 make every record a training record
              
            # TODO: check that each training set combination (n-1 folds) has at least one record per class level. 
            # Can be done by checking that each fold has less records of a certain level than the overall number of records for this level.
            # If not, permute anew
        } else {
            #=============================================
            # DIVIDE INTO TRAINING SET / VALIDATION SET
            #=============================================
            if (opts$TST.kind=="rand") {   # make a (random) division
              if (is.null(opts$TST.valiFrac))
                stop(sprintf("tdmModCreateCVindex: opts$TST.valiFrac is not defined"));
              if (is.null(opts$TST.trnFrac)) opts$TST.trnFrac = 1-opts$TST.valiFrac;
              if (opts$TST.valiFrac >= 1 | opts$TST.valiFrac < 0)
                stop(sprintf("tdmModCreateCVindex: opts$TST.valiFrac must be in [0,1). Current value is opts$TST.valiFrac = %f", opts$TST.valiFrac));
              if (opts$TST.trnFrac > 1 | opts$TST.trnFrac <= 0)
                stop(sprintf("tdmModCreateCVindex: opts$TST.trnFrac must be in (0,1]. Current value is opts$TST.trnFrac = %f", opts$TST.trnFrac));
              if (opts$TST.valiFrac+opts$TST.trnFrac>1)
                stop(sprintf("tdmModCreateCVindex: opts$TST.valiFrac+opts$TST.trnFrac > 1. Current value is opts$TST.valiFrac = %f"
                            , opts$TST.valiFrac,", opts$TST.trnFrac = %f", opts$TST.trnFrac ));
              if (stratified) {
                # ** NEW 06/2011 ** the division is done by ***stratified*** random sampling (recommended for classification):
                cat1(opts,opts$filename,": Stratified random training-validation-index with opts$TST.valiFrac = ",opts$TST.valiFrac*100,"%\n");
                rv <- dset[,response.variables[1]];
                lrv = length(response.variables);
                if (lrv>1) warning(sprintf("Stratified sampling is only done w.r.t. 1st response variable. It is not guaranteed to work for all %d response variables",lrv))  
                # calculate tfr, the number of training set records for each level of the response variable            
                tfr <- sapply(unique(rv),function(x) { round(opts$TST.trnFrac*length(which(rv==x))) });
                # (this seems complicated, but the simpler command: tfr <- round((1-opts$TST.valiFrac)*table(rv));
                # does not work, because 'table' orders the levels alphabetically but 'unique' (or 'strata') below 
                # requires them in the order they appear in column rv.) 
                #
                tfr[tfr<1] <- 1;      # ensure that there is at least one record for each class 
                cvi <- rep(-1,L);
                urv <- unique(rv);
                for (i in 1:length(urv))  cvi[ sample(which(rv==urv[i]), tfr[i]) ] <- 0;
                #
                #--- OLD and slow ---: 
                # the code below with strata from package sampling does the same as the two lines above,
                # but it is prohibitively slow when dset gets larger(50000 rows and more)
                #require(sampling);
                #s2=strata(dset,c(response.variables[1]),size=tfr, method="srswor");
                #cvi[s2$ID_unit] <- 0;          # set the training records   
                
                # select opts$TST.valiFrac of the non-training data as validation data:
                idxCviMinus <- which(cvi==-1);
                p <- sample(length(idxCviMinus));
                vfr <- opts$TST.valiFrac*L;     # index where the validation data ends
                if (trunc(vfr)>0) cvi[idxCviMinus[p[1:vfr]]] <- 1;
                
              } else {  # i.e. stratified=FALSE
                # simple random sampling (recommended for regression):
                p <- sample(L)                  # random permutation of indices 1:L  
                # calculate tfr, the record where the validation set starts (opts$TST.valiFrac)
                tfr <- opts$TST.trnFrac*L;      # index where the training data ends
                vfr <- (1-opts$TST.valiFrac)*L; # index where the validation data starts
                cat1(opts,opts$filename,": Random training-validation-index with opts$TST.valiFrac = ",opts$TST.valiFrac*100,"%\n")
                cvi <- rep(-1,L);        
                if (trunc(tfr)==0) stop(sprintf("No training data! opts$TST.trnFrac=%7.2f",opts$TST.trnFrac));
                cvi[p[1:tfr]] <- 0;                         # training set index ( opts$TST.trnFrac  percent of the data) 
                if (trunc(vfr)<L) cvi[p[(vfr+1):L]] <- 1;   # validation set index ( opts$TST.valiFrac  percent of the data)
              }
            } # opts$TST.kind=="rand"
            else           # i.e. opts$TST.kind=="col"
            {              # take the validation-training-division as delivered in dset[,opts$TST.COL]
              if (is.null(opts$TST.COL))
                stop(sprintf("tdmModCreateCVindex: opts$TST.COL is NULL, but opts$TST.kind=='col'"));
              if (!(opts$TST.COL %in% names(dset)))
                stop(sprintf("tdmModCreateCVindex: Data frame dset does not contain a column opts$TST.COL named \"%s\"", opts$TST.COL));
              cat1(opts,opts$filename,": Using training-validation-index from column",opts$TST.COL,"\n");
              cvi <- dset[,opts$TST.COL];
              if (!is.null(opts$TST.trnFrac)) {
                cat1(opts,opts$filename,": Selecting only ",opts$TST.trnFrac*100,"%% of the available training data for the actual training \n");
                cat1(opts,opts$filename,": (If you are in unbiased run and want all available training data, set tdm$TST.trnFrac to NULL)\n");
                # deselect 1-opts$TST.trnFrac of  all training data (set their cvi to -1):
                idxCviNull <- which(cvi==0);
                L0 = length(idxCviNull);
                p <- sample(L0);
                vfr <- opts$TST.trnFrac*L0;     # index where the training data to use ends
                if (trunc(vfr)>0) cvi[idxCviNull[p[(vfr+1):L0]]] <- -1;
              }
            }
        }

        cvi;
}

######################################################################################
# adjust sampsize (RF) in case of classification: 
#   - if samp==NULL, set to 3000
#   - if samp is scalar, clip it to total number of training records
#   - if samp is vector of length n.class, clip each element to number of training 
#     records for that class level
######################################################################################
tdmModAdjustSampsizeC <- function(samp, to.model, response.variable, opts) {
    if (is.null(samp)) samp=3000;
    if (length(samp)==1) newsamp = min(nrow(to.model),samp)
    else {
      newsamp <- as.vector(summary(to.model[,response.variable]));
      n.class <- length(newsamp);
      #browser()
      if (length(samp)!=n.class)
        stop(sprintf("Wrong length of sampsize 'samp'! It must be either a scalar or a vector of length n.class=%d.\n  samp = ",n.class),
             sprintf("%d ",samp));
      newsamp <- pmin(newsamp,samp);    # clip any element of samp which is larger than the respective 
                                        # number of training records (newsamp)
    }
    
    if (any(is.na(newsamp))) 
      stop(sprintf("sampsize 'samp' contains NA's! Consider appropriate lines 'opts$RF.samp[i] = ...' in APD file.\n  samp = "),
           paste(samp,collapse=", "))
    if (any(samp!=newsamp)) cat1(opts,"Clipping sampsize to ",newsamp,"\n");
    newsamp;
}

######################################################################################
# adjust sampsize (RF) in case of regression: 
#   - if samp==NULL, set to 3000
#   - clip it to total number of training records
######################################################################################
tdmModAdjustSampsizeR <- function(samp, to.model, response.variable, opts) {
  if (is.null(samp)) samp=3000;
  if (length(samp)==1) newsamp = min(nrow(to.model),samp)
  else {
    stop(sprintf("Wrong length of sampsize! It must be a scalar, but it is samp = "),
         sprintf("%d ",samp));
  }
  
  if (any(samp!=newsamp)) cat1(opts,"Clipping sampsize to ",newsamp,"\n");
  newsamp;
}

######################################################################################
# adjust cutoff (optional parameter for classification with Random Forest or SVM)
#   a) if length(cutoff)=n.class-1,  add cutoff[n.class] as the remainder to 1.
#   b) if cutoff[w] < 0 for exactly one w, then set cutoff[w] as the remainder to 1.
#      (If the others have a sum >= 1, then issue a warning and scale them to sum 0.9.)
#   c) if sum(cutoff) != 1, scale it to 1
#   d) if sum(cutoff) > 1 by a small amount (<=1e-8), reduce it to 1 or smaller.
######################################################################################
tdmModAdjustCutoff <- function(cutoff,n.class,text="cutoff")
{
    if (!is.null(cutoff)) {
      if (any(is.na(cutoff))) 
        stop(sprintf("cutoff '%s' contains NA's! Consider appropriate lines '%s[i] = ...' in APD file.\n  %s = ",text,text,text),
             paste(cutoff,collapse=", "))
      
      if (length(cutoff)==n.class-1) cutoff=c(cutoff,-1);   # assume that the last element is to be adjusted as remainder to 1
      if (length(cutoff)!=n.class) stop(sprintf("length(cutoff) differs from n.class=%d. cutoff = ",n.class),sprintf("%7.3f ",cutoff));
      w = which(cutoff<0);
      if (length(w)>1) stop("cutoff has more than one element < 0.\n cutoff = ",sprintf("%7.3f ",cutoff));
      if (length(w)==1) {
          s = sum(cutoff[-w]);
          if (s>=1) {
            warning(sprintf("One element of %s is not specified (<0), but the others have a sum>=1, there is no remainder.  ",text),
                    sprintf("%s = c(%s).\n  ",text,paste(cutoff,collapse=",")),
                    "Reducing the sum of the others to 0.9.")
            browser()
            cutoff=cutoff/s*0.9;
          }
          cutoff[w] = 1-sum(cutoff[-w]);
      }
      if (abs(sum(cutoff)-1)>1e-8) cutoff = cutoff/sum(cutoff);
      #
      # Due to rounding inaccurarcies, it can happen that sum(cutoff)=1+ 2.2e-16,
      # i.e. slightly larger than 1. randomForest does not allow this,
      # sum(cutoff) may not be larger than 1. The following code enforces this
      # by modifying cutoff as little as possible:
      eps = sum(cutoff)-1;
      if (eps>0) {
        if (eps>1e-8) stop("Something wrong with eps");   # this should not happen
        eps = max(1.5*eps,1e-10);
        cutoff = cutoff-eps/n.class;
        if (sum(cutoff)>1) stop("Something wrong, sum(cutoff) is still >1 ");
      }
    }
    
    cutoff;
}

#--- the old, now deprecated version ----------------
######################################################################################
# adjust cutoff (optional parameter for Random Forest and SVM)
#   a) if length(cutoff)=n.class-1,  add cutoff[n.class] as the remainder to 1.
#   b) if sum(cutoff) != 1, scale it to 1
#   c) if sum(cutoff) > 1 by a small amount (<=1e-8), reduce it to 1 or smaller.
######################################################################################
tdmModAdjustCutoff.OLD <- function(cutoff,n.class)
{
    if (!is.null(cutoff)) {
      if (length(cutoff)!=n.class) {
        if (length(cutoff)!=n.class-1)
          stop("Length of cutoff is not n.class-1!");
        cutoff[n.class] = 1-sum(cutoff[1:(n.class-1)]);
      }
      if (abs(sum(cutoff)-1)>1e-8) cutoff = cutoff/sum(cutoff);
      #
      # Due to rounding inaccurarcies, it can happen that sum(cutoff)=1+ 2.2e-16,
      # i.e. slightly larger than 1. randomForest does not allow this,
      # sum(cutoff) may not be larger than 1. The following code enforces this
      # by modifying cutoff as little as possible:
      eps = sum(cutoff)-1;
      if (eps>0) {
        if (eps>1e-8) stop("Something wrong with eps");
        eps = max(1.5*eps,1e-10);
        cutoff = cutoff-eps/n.class;
        if (sum(cutoff)>1) stop("Something wrong, sum(cutoff) is still >1 ");
      }
    }

    cutoff;
}

######################################################################################
# tdmModSortedRFimport:
#
#'   Sort the input variables decreasingly by their RF-importance.
#'
#'       Build a Random Forest using \code{importance=TRUE}. Usually the RF is smaller (50 trees), to speed up computation.
#'       Use na.roughfix for missing value replacement.
#'       Decide which input variables to keep and return them in SRF$input.variables
#'
#'   @param d_train   training set
#'   @param response.variable   the target column from \code{d_train} to use for the RF-model
#'   @param input.variables   the input columns from \code{d_train} to use for the RF-model
#'   @param opts options, here we use the elements [defaults in brackets]:
#'    \itemize{
#'     \item SRF.kind:  \cr
#'          ="xperc": keep a certain importance percentage, starting from the most important variable \cr
#'          ="ndrop": drop a certain number of least important variables \cr
#'          ="nkeep": keep a certain number of most important variables \cr
#'          ="none": do not call \code{\link{tdmModSortedRFimport}} at all (see tdmRegress.r and tdmClassify.r)
#'     \item SRF.ndrop:   [0] how many variables to drop (if SRF.kind=="ndrop")
#'     \item SRF.XPerc:   [0.95] if >=0, keep that importance percentage, starting with the most
#'                   important variables (if SRF.kind=="xperc")
#'     \item SRF.calc:    [TRUE] =TRUE: calculate importance & save on SRF.file, =F: load from SRF.file
#'                   (SRF.file = Output/<filename>.SRF.<response.variable>.Rdata)
#'     \item SRF.ntree:   [50] number of RF trees
#'     \item SRF.verbose: [2]
#'     \item SRF.maxS:    [40] how many variables to show in plot
#'     \item SRF.minlsi:  [1] a lower bound for the length of SRF$input.variables
#'     \item RF.sampsize: sampsize for RF, set prior to calling this func via tdmModAdjustSampsize(opts$SRF.samp,...)
#'     \item GD.DEVICE:   if !="non", then make a bar plot on current graphic device
#'     \item CLS.CLASSWT: class weight vector to use in random forest training
#'    }
#' @return \code{SRF},    a list with the following elements:
#'     \item{input.variables}{   the vector of input variables which remain after importance
#'                    processing. These are sorted by decreasing importance.}
#'     \item{s_input}{all input.variables sorted by decreasing (**NEW**) importance}
#'     \item{s_imp1}{ the importance for s_input}
#'     \item{s_dropped}{   vector with name of dropped variables}
#'     \item{lsd}{    length of s_dropped}
#'     \item{perc}{   the percentage of total importance which is in the dropped variables}
#'     \item{opts}{   some defaults might have been added}
#'
#' @author Wolfgang Konen, Patrick Koch \email{wolfgang.konen@@fh-koeln.de}
#' @export
######################################################################################
tdmModSortedRFimport <- function(d_train, response.variable, input.variables, opts)
{
    # require(randomForest);    # now via direct call 'randomForest::'
    opts <- tdmOptsDefaultsSet(opts);

    ptm <- proc.time();
    filename <- opts$filename;
#--- this is now in tdmEnvTMakeNew.r ---
#    dir.output <- paste(dirname(opts$dir.output),basename(opts$dir.output),sep="/")  # remove trailing "/", if it exists
#    if (!file.exists(dir.output)) {
#      success = dir.create(dir.output);
#      if (!success) stop(sprintf("Could not create dir.output=%s",dir.output));
#    }
#    SRF.file <- paste(paste(dir.output,filename,sep="/"),"SRF",response.variable,"Rdata",sep=".")
    

    ############################################################
    # PART 1: Calculate SRF importance (or load it from opts)
    ############################################################
    if (opts$SRF.calc==TRUE) {
      if (opts$SRF.kind=="xperc") {
        if (opts$SRF.XPerc<0) stop(sprintf("opts$SRF.XPerc < 0 : %f",opts$SRF.XPerc));
        if (opts$SRF.XPerc>1) stop(sprintf("opts$SRF.XPerc > 1 : %f",opts$SRF.XPerc));
      }
      
      formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
      to.model <- d_train[,c(response.variable,input.variables)]
      cat1(opts,filename,": Train RF (importance, sampsize=", opts$RF.sampsize,") ...\n")
      if (!is.null(opts$CLS.CLASSWT)) {
        cat1(opts,"Class weights: ", opts$CLS.CLASSWT,"\n")
        #cwt = opts$CLS.CLASSWT*1;# strange, but necessary: if we omit '*1' then cwt seems to be a copy-by-reference of
                                  # opts$CLS.CLASSWT. After randomForest call  cwt is changed and also opts$CLS.CLASSWT would be changed (!)
                                  # (Another way to hinder RF to change the parameter 'classwt' on output is to make this vector a *named*
                                  # vector. This is actually the case, because in tdmClassify we decorate opts$CLS.CLASSWT with the 
                                  # level names of response.variable)
        cwt = sprintf("classwt=c(%s)",paste(opts$CLS.CLASSWT,collapse=","))
      } else { cwt="classwt=NULL"; }
      if (!is.null(opts$SRF.cutoff)) {
        cat1(opts,"Cutoff: ", opts$SRF.cutoff,"\n");
        #cutoff = sprintf("cutoff=c(%s)",paste(opts$SRF.cutoff,collapse=","))  # DON*T!! may lead to "Incorrect cutoff" error due to rounding inaccuracies
      } else { cutoff=NULL; }
      
       res.SRF = switch(opts$SRF.method
        ,"RFimp" = fsRfImportance(opts,formul,to.model,cwt)
        #,"lasso" = fsLasso(opts)
        ,"INVALID"
        ); 
      if (res.SRF[1]=="INVALID1") {
        stop(sprintf("*** Invalid FS method=%s ***\n",opts$SRF.method));
      } 
      if (is.null(res.SRF)) 
        stop(sprintf("FS method %s did not return a suitable result 'res.SRF'",opts$SRF.method));
      
      imp1 <- res.SRF$imp1;
      
      s_input <- input.variables[order(imp1,decreasing=TRUE)]  # input.variables sorted by increasing importance
      s_imp1 <- imp1[order(imp1,decreasing=TRUE)]
      s_imp2 <- pmax(s_imp1,0);
      # It can happen that unimportant variables have a negative importance s_imp1. This is only a statistical
      # artefact, the importance is compatible with 0 (see notes_RF.doc for details). Since cumsum below does not  
      # allow negative values, we make a copy s_imp2 with all negative values clipped to 0.
      
      cat1(opts,filename,": Saving SRF (sorted RF) importance info on opts ...\n");
      opts$srf[[response.variable]] = data.frame(s_input=s_input
                                                ,s_imp1=s_imp1
                                                ,s_imp2=s_imp2 
                                                );
    }
    else {   # i.e. if (opts$SRF.calc==F)
      cat1(opts,filename,": Loading sorted RF importance info from opts ...\n")
      #if (!file.exists(SRF.file)) stop(sprintf("tdmModSortedRFimport: SRF.file=%s does not exist"));
      #load(file=SRF.file);
      if (!any(names(opts$srf)==response.variable))
        stop(sprintf("list opts$srf contains no data frame named %s. Consider opts$SRF.calc==TRUE.",response.variable));
      s_input=opts$srf[[response.variable]]$s_input;
      s_imp1 =opts$srf[[response.variable]]$s_imp1;
      s_imp2 <- pmax(s_imp1,0);
    } # if (opts$SRF.calc)
    if (is.null(s_input) | is.null(s_imp1)) 
      stop("opts$SRF.s_input or opts$SRF.s_imp1 are not available. Consider opts$SRF.calc=TRUE.");
    
    ############################################################
    # PART 2: Select important variables
    ############################################################
    if (opts$SRF.kind=="ndrop") {
        # remove the SRF.ndrop input variables which have the lowest importance
        # (but keep at least one):
        opts$SRF.ndrop <- min(opts$SRF.ndrop,length(s_input)-1)
        lsi <- length(s_input)-opts$SRF.ndrop;
    } else {    # i.e. opts$SRF.kind=="xperc" or "nkeep"
        if (opts$SRF.kind=="nkeep") {
          opts$SRF.nkeep <- lsi <- min(opts$SRF.nkeep,length(s_input));
        } else {    # i.e. opts$SRF.kind=="xperc"
          # keep the minimal number of those most important input variables which
          # sum up together to more than opts$SRF.XPerc * (sum of all importances)
          w = which(cumsum(s_imp2)/sum(s_imp2)<=opts$SRF.XPerc)
          lsi <- max(w)+1;
          lsi <- min(lsi,length(s_imp1)); # for the case SRF.XPerc==1: lsi may not be larger than #input.variables
        }
    }
    # don't let lsi = length(input.variables) become smaller than minlsi
    minlsi <- min(length(s_input),opts$SRF.minlsi);
    lsi <- max(lsi,minlsi);
    input.variables <- as.character(s_input[1:lsi])

    s_dropped <- setdiff(s_input,input.variables);
    lsd <- length(s_dropped);
    # what percentage of total importance is in the dropped variables:
    perc <- (1-sum(s_imp2[1:(length(s_imp2)-lsd)])/sum(s_imp2))*100;
    if (opts$SRF.kind=="xperc" & (perc/100>1-opts$SRF.XPerc | perc<0)) {
        stop(sprintf("Something wrong with perc: %f",perc));
    }

    if (opts$SRF.verbose>=2) {
        lmi <- length(s_input);
        cat(sprintf("Variables sorted by importance (%d %s):\n", lmi,
                    ifelse(lmi>100,", we show first 100 only","")));
        print(s_input[1:min(lmi,100)]);
        if (opts$SRF.kind=="ndrop" | (opts$SRF.kind=="xperc" & opts$SRF.XPerc>0.5)) {
          cat(sprintf("Dropped columns (%d [=%5.1f%% of total importance]%s):\n", lsd, perc,
                      ifelse(lsd>100,", we show first 100 only","")));
          if (lsd>0) print(s_dropped[1:min(lsd,100)]);
        } else { # i.e. opts$SRF.kind=="nkeep" or (..=="xperc" & SRF.XPerc <= 0.5)
          cat(sprintf("Kept columns (%d [=%5.1f%% of total importance]%s):\n", lsi, 100-perc,
                      ifelse(lsi>100,", we show first 100 only","")));
          if (lsi>0) print(s_input[1:min(lsi,100)]);
        }
        cat(sprintf("Proc time: %5.2f\n",(proc.time()-ptm)[1]));
    }

    ############################################################
    # PART 3: Plot the importance
    ############################################################
    if (opts$GD.DEVICE!="non") {
      maxS <- min(length(s_imp1),opts$SRF.maxS)
      tdmGraphicNewWin(opts)
      oldmar = par()$mar
      if (opts$GD.DEVICE!="rstudio") {  # some barplots cause an "Error in plot.new() : figure margins too large"
                                        # on device RStudioGD >> we omit barplot if in RStudio
        par(mar=c(10,3,2,2)+0.1)
        barplot(t(s_imp1[1:maxS]), names.arg=s_input[1:maxS], las=3, cex.lab=0.75, col=2, main=paste(filename, response.variable, sep=" : "))
        par(mar=oldmar)   
      }
      tdmGraphicCloseWin(opts);
      if (opts$SRF.calc==TRUE & opts$SRF.method=="RFimp") {
        tdmGraphicNewWin(opts)
        randomForest::varImpPlot(res.SRF,n.var=maxS)    # NEW: plot both MeanDecreaseAccuracy and MeanDecreaseGini
        tdmGraphicCloseWin(opts);
      }
    }


    SRF = list(input.variables=input.variables
              , s_input=s_input
              , s_imp1=s_imp1
              , s_dropped=s_dropped   # vector with name of dropped variables
              , lsd=lsd               # length of s_dropped
      	      , perc=perc             # the percentage of total importance which is in the dropped variables
      	      , opts=opts             # some elements might have been added
              );
}

# fsRfImportance: helper for tdmModSortedRFimport
fsRfImportance <- function(opts,formul, to.model, cwt) {
  # require(randomForest);    # for 'na.roughfix', now via direct call 'randomForest::'
  testit::assert("Cutoff is bigger than 1",sum(opts$SRF.cutoff)<=1)
  
  # we work here with a command text string and eval(...) to allow for the presence or
  # absence of certain options like "mtry" or "cutoff" which are not allowed to be NULL:
  rf.options = "ntree=opts$SRF.ntree, importance=TRUE"; # NEW: only with 'importance=TRUE' we see MeanDecreaseAccuracy,
  # otherwise we get only MeanDecreaseGini (faster, but sometimes unreliable)
  rf.options = paste(rf.options,"sampsize=opts$RF.sampsize",sep=",")
  rf.options = paste(rf.options,cwt,"na.action=randomForest::na.roughfix","proximity=F",sep=",")
  if (!is.null(opts$SRF.mtry)) rf.options = paste(rf.options,"mtry=opts$SRF.mtry",sep=",")
  if (!is.null(opts$SRF.cutoff)) rf.options = paste(rf.options,"cutoff=opts$SRF.cutoff",sep=",")
  #
  #if (!is.null(opts$SRF.cutoff)) rf.options = paste(rf.options,cutoff,sep=",")  # DON'T!! this may lead to "Incorrect cutoff" error due to rounding!
  #
  #dbg_chase_cutoff_bug(formul,to.model,d_train,response.variable,rf.options,opts);
  #print(.Random.seed[1:6]);
  flush.console();
  res.SRF <- NULL;      
  #browser()
  eval(parse(text=paste("res.SRF <- randomForest::randomForest( formul, data=to.model,",rf.options,")")));
  # select MeanDecreaseAccuracy-importance (NEW 05/11, together with switch 'importance=TRUE' above)
  res.SRF$imp1 <- randomForest::importance(res.SRF, type=1, scale=opts$SRF.scale);      
  #print(.Random.seed[1:6]);
  
  # only for debug or in-depth-analysis:
  #analyzeImportance(res.SRF,input.variables,opts);
  
  res.SRF;
}

#      fsLasso <- function(opts) {
#        # *** does not yet work in all cases !!! ***
#        require(lasso2);
#        formul <- formula(paste("as.numeric(",response.variable,") ~ ."))   # use all possible input variables
#        res.SRF <- NULL;
#        #browser()
#        eval(parse(text=paste("res.SRF <- l1ce( formul, data=to.model, bound=0.5)")));
#        res.SRF$imp1 <- as.matrix(res.SRF$coefficients[-1]);
#        
#        res.SRF;
#      }

# helper for fsRfImportance 
# (only for debug or in-depth-analysis of the importance)
analyzeImportance <- function(res.SRF,input.variables,opts)  {
    #
    # 1st finding: variables which appear never as split variable have necessarily importance=0
    #
    splitVar=NULL;              # the union of all variables which appear ever as a split in a tree of the forest 
    splitLst=NULL;              # the list of all split variables (variables may appear more than once)
    rn=rownames(res.SRF$imp1);
    for (k in 1:res.SRF$ntree) splitVar=union(splitVar,randomForest::getTree(res.SRF,k=k,labelVar=T)$"split var")
    splitVar = sort(splitVar); # remove the <NA>
    noSplitVar = sort(setdiff(input.variables,splitVar))
    for (k in 1:res.SRF$ntree) splitLst=c(splitLst,as.character(randomForest::getTree(res.SRF,k=k,labelVar=T)$"split var"))
    splitLst=sort(splitLst);
    impZeroVar = sort(rn[res.SRF$imp1==0])
    cat(sprintf("%d variables have importance exactly zero\n",length(impZeroVar)));
    cat(sprintf("%d variables never appear in a split\n",length(noSplitVar)));
    if (length(setdiff(noSplitVar,impZeroVar))==0) cat("All variables never appearing in a split have importance exactly zero\n")
    cat(sprintf("setdiff(impZeroVar,noSplitVar) is: ")); print(setdiff(impZeroVar,noSplitVar))
    
    #
    # 2nd finding: variables with negative importance have a sligthly smaller split count than those with positive importance
    #         Most variables with negative importance have a very small raw importance (%IncMSE), it is the division by 
    #         importanceSD which makes the negative value big.
    #
    impNegVar = sort(rn[res.SRF$imp1<0]);
    impPosVar = setdiff(splitVar,impNegVar);
    cat("Split count for variables with negative importance (we show first 6 only):\n")
    for (v in impNegVar[1:min(length(impNegVar),6)]) cat(sprintf("%s: %d  %e\n", v,length(which(splitLst==v)),res.SRF$imp1[rownames(res.SRF$imp1)==v]))
    avg=0;  for (v in impNegVar) avg = avg + length(which(splitLst==v)); 
    cat("Avg. count:", avg/length(impNegVar),"\n");
    cat("Split count for variables with positive importance (we show first 6 only):\n")
    for (v in impPosVar[1:min(length(impPosVar),6)])  cat(sprintf("%s: %d  %f\n", v,length(which(splitLst==v)),res.SRF$imp1[rownames(res.SRF$imp1)==v]))
    avg=0;  for (v in impPosVar) avg = avg + length(which(splitLst==v)); 
    cat("Avg. count:", avg/length(impPosVar),"\n");
    #
    par(mfcol=c(1,2))
    ind1 = order(res.SRF$imp1)
    plot(res.SRF$imp1[ind1],ylim=c(-2,3))
    points(res.SRF$importance[ind1,1]*15,col="green")
    ind2 = order(res.SRF$importance[,1])
    plot(res.SRF$importance[ind2,1]*10,col="green",col.lab="green",ylim=c(-2,3))
    points(res.SRF$imp1[ind2])
    par(mfcol=c(1,1))
    
    #
    # 3rd finding: variables with |importance|=1.01... (exactly the same value) have mostly a split count of exactly 1.
    #         The quotient  %IncMSE/importanceSD   can be shown to be exactly +-sqrt(N/(N-1)) in this case.
    #
    N=opts$SRF.ntree;
    impOneVar =  sort(rn[abs(res.SRF$imp1)==sqrt(N/(N-1))]);     # all variables with |importance| = 1.010153 in case N=50
    cat(sprintf("Split count for variables with |importance| =%9.6f (we show first 6 only):\n",sqrt(N/(N-1)) ));
    for (v in impOneVar[1:min(length(impOneVar),6)])  cat(sprintf("%s: %d  %f\n", v,length(which(splitLst==v)),res.SRF$imp1[rownames(res.SRF$imp1)==v]))
    avg=0;  for (v in impOneVar) avg = avg + length(which(splitLst==v)); 
    cat("Avg. count:", avg/length(impOneVar),"\n");
#browser()        
}

######################################################################################
# tdmModVote2Target
#
#'     Analyze how the vote fraction corresponds to reliability of prediction.
#' 
#' This function analyzes whether in different vote bins the trained RF makes
#' predictions with different reliability. Only for RF-prediction in case of 
#' binary (0/1) classification. \cr\cr Expected result: The larger the fraction
#' of trees voting for class 0 is, the smaller is the percentage of true class-1-
#' cases in this vote bin.
#' This function is somewhat specialized for the DMC2010-task.
#'
#' @param  vote0     vector: which fraction of trees votes for class 0?
#' @param  pred      vector: the predicted class for each record (0/1)
#' @param  target    vector: the true class for each vector (0/1)
#' 
#' @return a data frame with columns 
#'    \item{vcut}{  vote cut v}
#'    \item{count}{ number of cases with vote fraction in [v[i-1],v[i]]}
#'    \item{pred0}{ fraction of 0-predictions}
#'    \item{pCorr}{ fraction of correct predictions}
#'    \item{pR}{    fraction of true 1-cases}
#'
#' @author Wolfgang Konen \email{wolfgang.konen@@fh-koeln.de}
#' @export
######################################################################################
tdmModVote2Target <- function(vote0,pred,target) {
    tvot <- data.frame(v0=vote0,pred=pred,targ=target);
    vcut <- c(0.4,0.5,0.6,0.65,0.70,0.73,0.75,0.80,0.85,0.90,0.95,1.0)
    R <- nrow(tvot)
    tres <- data.frame()
    i=1; vprev=0;
    for (v in vcut) {
      ind <- which(vprev< tvot$v0 & tvot$v0 <= v)
      I <- length(ind);
      tres[i,"vcut"] <- v;
      tres[i,"count"] <- I;
      tres[i,"frac"] <- I/R;
      tres[i,"pred0"] <- length(which(tvot$pred[ind]==0))/I;
      tres[i,"pCorr"] <- length(which(tvot$targ[ind]==tvot$pred[ind]))/I;
      tres[i,"pR"] <- length(which(tvot$targ[ind]==1))/I;
      i=i+1;  vprev=v;
    }
    tres;
}

######################################################################################
# tdmModConfmat
#
#'     Calculate confusion matrix, gain and RGain measure. 
#'
#'  @param  d 		      data frame
#'  @param  colreal     name of column in d which contains the real class
#'  @param  colpred 		name of column in d which contains the predicted class
#'  @param  opts        a list from which we use the elements: \itemize{
#'     \item \code{gainmat}:   the gain matrix for each possible outcome, same size as \code{cm$mat} (see below). \cr
#'               \code{gainmat[R1,P2]} is the gain associated with a record of real class R1 which we
#'               predict as class P2. (gain matrix = - cost matrix)
#'     \item \code{rgain.type}: one out of \{"rgain" | "meanCA" | "minCA" | "arROC" | "arLIFT" | "arPRE" \},
#'               affects output \code{cm$mat} and \code{cm$rgain}, see below.
#'    }
#'  @param  predProb    if not NULL, a data frame with as many rows as data frame \code{d}, containing columns 
#'               (index, true label, predicted label, prediction score). Is only needed for \code{opts$rgain.type=="ar*"}.
#'
#'  @return \code{cm}, a list containing: 
#'     \item{mat}{ matrix with real class levels as rows, predicted class levels columns.  \cr
#'               \code{mat[R1,P2]} is the number of records with real class R1
#'               predicted as class P2, if opts$rgain.type=="rgain".
#'               If opts$rgain.type=="meanCA" or "minCA", then show this number as percentage
#'               of "records with real class R1" (percentage of each row).
#'               CAUTION: If colpred contains NA's, those cases are missing in mat (!)
#'               (but the class errors are correct as long as there are no NA's in colreal)}
#'     \item{cerr}{ class error rates, vector of size nlevels(colreal)+1.  \cr
#'               \code{cerr[X]} is the misclassification rate for real class X. \cr
#'               \code{cerr["Total"]} is the total classification error rate.}
#'     \item{gain}{ the total gain (sum of pointwise product \code{opts$gainmat*cm$mat})  }
#'     \item{gain.vector}{ gain.vector[X] is the gain attributed to real class label X.
#'               gain.vector["Total"] is again the total gain.}
#'     \item{gainmax}{    the maximum achievable gain, assuming perfect prediction}
#'     \item{rgain}{      ratio gain/gainmax in percent, if \code{opts$rgain.type=="rgain"};  \cr
#'               mean class accuracy percentage (i.e. mean(diag(cm$mat)), if \code{opts$rgain.type=="meanCA"}; \cr
#'               min class accuracy percentage (i.e. min(diag(cm$mat)), if \code{opts$rgain.type=="minCA"}; \cr
#'               area under ROC curve (a number in [0,1]), if \code{opts$rgain.type=="arROC"}; \cr
#'               area between lift curve and horizontal line 1.0, if \code{opts$rgain.type=="arLIFT"}; \cr
#'               area under precision-recall curve (a number in [0,1]), if \code{opts$rgain.type=="arPRE"}; \cr
#'               } 
#' 
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), Patrick Koch 
#' @seealso  \code{\link{tdmClassify}}    \code{\link{tdmROCRbase}}
#' @export
######################################################################################
tdmModConfmat <- function(d,colreal,colpred,opts,predProb=NULL)
{
    if (is.null(opts$rgain.type)) opts$rgain.type="rgain";
    if (!(opts$rgain.type %in% c("rgain","meanCA","minCA","arROC","arLIFT","arPRE")))
      stop(sprintf("Invalid option opts$rgain.type=\"%s\"\n  Allowed values are \"rgain\",\"meanCA\",\"minCA\",\"arROC\",\"arLIFT\",\"arPRE\"",opts$rgain.type));
    AREA = opts$rgain.type %in% c("arROC","arLIFT","arPRE");
    if (AREA & is.null(predProb))
      stop(sprintf("Can not calculate the measure rgain for opts$rgain.type==%s, if argument predProb is NULL",opts$rgain.type));

    gainmat = opts$CLS.gainmat;
    col.real <- factor(d[,colreal],levels=colnames(gainmat));
    col.pred <- factor(d[,colpred],levels=colnames(gainmat));
    #cmat = matrix(0,nrow=nlevels(col.real),ncol=nlevels(col.pred),
    #              dimnames=list(levels(col.real),levels(col.pred)))
    #-- cmat now simpler to obtain with function table (see below)
    ccase = matrix(0,nrow=1,ncol=nlevels(col.real)+1,
                  dimnames=list("class cases:",c(levels(col.real),"Total")))

    # class cases
    for (rc in levels(col.real)) {
        ccase[1,rc] <-  length(which(col.real==rc))
        #d1 <- d[col.real==rc, ]      # all records with real class rc
        #ccase[1,rc] <- dim(d1)[1]
        #for (pc in levels(col.pred)) {
        #    cmat[rc,pc] <- length(which(d1[,colpred]==pc))
        #}
    }
    ccase[1,"Total"] <- sum(ccase)      # all cases, including those with NA in col.pred
                                        # (but excluding those with NA in col.real)
    # confusion matrix
    cmat <- table(actual=col.real,predicted=col.pred)
                           # each case with NA in colpred is in *no* cell of cmat

    # class errors
    cerr = matrix(0,nrow=1,ncol=nlevels(col.pred)+1,
                  dimnames=list("class errors:",c(levels(col.pred),"Total")))
    for (rc in levels(col.real)) {
        cerr[1,rc] <- 1-cmat[rc,rc]/ccase[1,rc]
    }
    cerr[1,"Total"] <- 1-sum(diag(cmat))/ccase[1,"Total"]

    gain = sum(gainmat*cmat)
    gain.vector = matrix(0,nrow=1,ncol=nlevels(col.pred)+1,
                  dimnames=list("gain.vector",c(levels(col.pred),"Total")))
    for (rc in levels(col.real)) {
        gain.vector[1,rc] <- sum(gainmat[rc,] * cmat[rc,])
    }
    gain.vector[1,"Total"] <- sum(gain.vector)

    gainmax = sum(apply(gainmat,1,max)*rowSums(cmat))
    # the 1st term is a vector containing the max gain for each row (true class)
    # the 2nd term is a vector containing the #records for each true class

    rgain = gain/gainmax*100;
    nacase = length(which(is.na(col.pred)));
    if (sum(cmat)+nacase!=ccase[1,"Total"])
        stop("tdmModConfmat: Something wrong in NA-counting!");

    if (opts$rgain.type %in% c("meanCA","minCA")) {
      for (rc in levels(col.real)) {
          cmat[rc,] <- cmat[rc,]/ccase[1,rc];
      }
      rgain=switch(opts$rgain.type
            , "meanCA" = mean(diag(cmat))
            , "minCA" = min(diag(cmat))
            );
    }
    if (AREA) {
      perf=switch(opts$rgain.type
            , "arROC" = tdmROCR_calc(predProb,"tpr","fpr")
            , "arLIFT" = tdmROCR_calc(predProb,"lift","rpp")
            , "arPRE" = tdmROCR_calc(predProb,"prec","rec")
            );
      rgain=switch(opts$rgain.type
            , "arROC" =, "arPRE" = tdmROCR_area(perf,"ROC")
            , "arLIFT" = tdmROCR_area(perf,"lift")
            );
    }

    cm = list( mat=cmat
              ,cerr=cerr
              ,ccase=ccase
              ,nacase=nacase
              ,gain=gain
              ,gain.vector=gain.vector
              ,gainmax=gainmax
              ,rgain=rgain
              );

    return(cm)
}

# debug info to chase the "Incorrect cutoff specified bug" in training randomForest
dbg_chase_cutoff_bug <- function(formul,to.model,d_train,response.variable,rf.options,opts) {
  print2(opts,c(Targets=NA,table(to.model[,response.variable])))
  dir.Rdata <- paste(dirname(opts$dir.Rdata),basename(opts$dir.Rdata),sep="/")  # remove trailing "/", if it exists
  if (!file.exists(dir.Rdata)) {
    success = dir.create(dir.Rdata);     
    if (!success) stop(sprintf("Could not create dir.Rdata=%s",dir.Rdata));
  }
  save(formul,to.model,response.variable,rf.options,opts,file=paste(dir.Rdata,"rf_input_dbg.Rdata",sep="/"));
  cat1(opts,"RF-debug-data saved to", paste(dir.Rdata,"rf_input_dbg.Rdata",sep="/"),"\n");
  if (!is.null(opts$SRF.cutoff))
    if (length(levels(d_train[,response.variable]))!=length(opts$SRF.cutoff)) {
      warning("Cutoff problems ahead!!",immediate.=TRUE)
      browser()
    }
    if (sum(opts$SRF.cutoff)-1>0) {
      warning(sprintf("Sorry, but sum(opts$SRF.cutoff) is larger than 1 by %f",sum(opts$SRF.cutoff)-1),immediate.=TRUE)
      browser()
    }
}


