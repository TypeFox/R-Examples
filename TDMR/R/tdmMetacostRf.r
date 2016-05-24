######################################################################################
# tdmMetacostRf: 
#     Function to wrap the MetaCost algorithm around RF (Random Forest)
# INPUT:
#     response.variable     name of target column
#     to.model              data frame with training data
#     opts                  list with options, here we need:
#         RF.ntree          number of trees for RF
#         RF.sampsize       sample size for RF (one number, NOT a vector)
#         RF.mtry           mtry for RF
#         CLASSWT           class weights 
#         cutoff            [default 1/n.class]
#         gainmat           the gain matrix, G[j,i] is the gain associated with 
#                           predicting "class==i" for a record of true class j.
#                           The average gain is to be maximized.
# OUTPUT:
#     res.rf      the random forest, trained to maximize the gain. 
#                 This object can be utilized in the standard way
#                           predict(res.rf, newdata=to.test)
#                 to perform gain-maximizing predictions on new test data.
#
# DETAILS: 
#   see [Domingos99] Pedro Domingos. MetaCost: A General Method for Making Classifiers Cost-Sensitive. 
#   In: Proceedings of the Fifth International Conference on Knowledge Discovery and Data Mining (KDD-99), 1999. 
#
# Author: Wolfgang Konen, FHK, May'2010
#
######################################################################################
tdmMetacostRf <- function(response.variable,to.model,opts)
{
        if (!is.null(opts$CLS.CLASSWT)) {
            cat1(opts,"Class weights: ", opts$CLS.CLASSWT,"\n")
            cwt = opts$CLS.CLASSWT*1.0;    # strange, but necessary: if we omit '*1' then cwt seems to be a copy-by-reference of
                                      # opts$CLS.CLASSWT. After randomForest call  cwt is changed and also opts$CLS.CLASSWT would be changed (!)
        } else { cwt=NULL; }
        if (!is.null(opts$CLS.cutoff)) cat1(opts,"Cutoff: ", opts$CLS.cutoff,"\n")
        formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
        orig.target <- as.factor(to.model[,response.variable]);
        
        # we work here with a command text string and eval(parse(...)) to allow for the presence or
        # absence of certain options like "mtry" or "cutoff" which are not allowed to be NULL.
        # The individual "eval(...)" on the following lines are for clarity of res.rf$call 
        # (it should read "..., ntree=400, ..." and not "..., ntree=opts$RF.ntree, ...") 
        # BUT: we cannot use "eval(...)" for the lists cutoff and cwt, therefore we use here
        # the 'non-speaking' variables cwt, opts$CLS.cutoff and add below res.rf$cutoff, res.rf$classwt 
        # for optional later reference or user inspection.
        rf.options = paste(" ntree=",eval(opts$RF.ntree));
        rf.options = paste(rf.options," sampsize=opts$RF.sampsize",sep=",")
        rf.options = paste(rf.options," classwt=cwt"," na.action=randomForest::na.roughfix"," proximity=F",sep=",")
        #if (!is.null(cwt))  paste(rf.options,paste("classwt=",eval(cwt)),sep=",")    # not run
        if (!is.null(opts$RF.mtry)) rf.options = paste(rf.options,paste(" mtry=",eval(opts$RF.mtry)),sep=",")
        if (!is.null(opts$CLS.cutoff)) rf.options = paste(rf.options," cutoff=opts$CLS.cutoff",sep=",")
        if (!is.null(opts$RF.nodesize)) rf.options = paste(rf.options,paste(" nodesize=",eval(opts$RF.nodesize)),sep=",")
        
        DBG=F
        if (DBG) {          # debug info to chase the "Incorrect cutoff specified bug": 
          print2(opts,c(Targets=NA,table(to.model[,response.variable])))      # debug info: 
          dir.Rdata <- paste(dirname(opts$dir.Rdata),basename(opts$dir.Rdata),sep="/")  # remove trailing "/", if it exists
          if (!file.exists(dir.Rdata)) {
            success = dir.create(dir.Rdata);     
            if (!success) stop(sprintf("Could not create dir.Rdata=%s",dir.Rdata));
          }
          save(formul,to.model,response.variable,rf.options,opts,file=paste(dir.Rdata,"rf_input_dbg.Rdata",sep=""));
          cat1(opts,"RF-debug-data saved to", paste(dir.Rdata,"rf_input_dbg.Rdata",sep=""),"\n");
          if (!is.null(opts$CLS.cutoff))
            if (length(levels(to.model[,response.variable]))!=length(opts$CLS.cutoff)) {
              warning("Cutoff problems ahead!!",immediate.=TRUE)
              browser()
            }
        }
        flush.console();          
	      
	      res.rf <- NULL;          # just to make "R CMD check" happy 
        # Train an initial RF to label the data in such a way that the structural gain is maximized 
        eval(parse(text=paste("res.rf <- randomForest::randomForest( formul, data=to.model,",rf.options,")")))

        Pjx = res.rf$votes        # we take the fractions of OOB-votes as an unbiased estimator for P(j|x)
                                  # Pjx is a matrix with #columns = #target labels and #rows = #records in to.model
    
        if (any(is.na(Pjx))) {    # if there are records which were never OOB (should happen very rarely, if ntree is high enough)
            # take the voting on all trees as surrogate:
            p.votes <- predict(res.rf, newdata=to.model, type="vote")
            ind = which(is.na(Pjx))
            Pjx[ind] = p.votes[ind]
        }
        
        # assign to each record that class level which maximizes the structural gain 
        # (minimizes the structural risk):
        mc.target <- factor(apply(Pjx %*% opts$CLS.gainmat, 1,  function(x) { colnames(Pjx)[which.max(x)] } )
                           , levels=levels(orig.target))
                    # this line implements the formula 
                    #       arg max_i ( sum_j ( P(j|x) * G[j,i] ) )
                    # where we return instead of the maximizing i its corresponding column name, i.e. the target label.
                    # (Pjx %*% opts$CLS.gainmat is a matrix with #rows = #records and the second argument "1" indicates 
                    # that apply should loop over the 1st dimension of this matrix, i.e. the rows.) 

        # assure that each leval of orig.target appears at least once in mc.target. If this were not the case, then 
        # the call to randomForest would crash with the cryptic error message "Incorrect cutoff specified" (!)                    
        for (cl in levels(orig.target)) {
          if (!any(mc.target==cl)) {
            w = which(orig.target==cl);
            mc.target[w[1]]=cl;         # a simple fix: make the first record of class cl in orig.target also cl in mc.target
          }                             # >> we guarantee that each cl appears at least once as target
        }
        to.model[,response.variable] <- mc.target;    # Only a local change. In the calling context, the 
                                                      # column response.variable remains untouched.
        # Train a new RF (model M) on the relabeled data. M will maximize the average gain, when evaluated on 
        # the original response.variable.
        
        DBG=F
        if (DBG) {          # debug info to chase the "Incorrect cutoff specified bug": 
          #if (!is.null(opts$CLS.cutoff)) cat2(opts,"Cutoff: ", opts$CLS.cutoff,"\n")
          print2(opts,c(Targets=NA,table(to.model[,response.variable])))      # debug info: 
          dir.Rdata <- paste(dirname(opts$dir.Rdata),basename(opts$dir.Rdata),sep="/")  # remove trailing "/", if it exists
          if (!file.exists(dir.Rdata)) {
            success = dir.create(dir.Rdata);     
            if (!success) stop(sprintf("Could not create dir.Rdata=%s",dir.Rdata));
          }
          save(formul,to.model,response.variable,rf.options,opts,file=paste(dir.Rdata,"rf_input_dbg.Rdata",sep=""));
          cat1(opts,"RF-debug-data saved to", paste(dir.Rdata,"rf_input_dbg.Rdata",sep=""),"\n");
          if (!is.null(opts$CLS.cutoff))
            if (length(levels(to.model[,response.variable]))!=length(opts$CLS.cutoff)) {
              warning("Cutoff problems ahead!!",immediate.=TRUE)
              browser()
            }
          flush.console();
        }

        eval(parse(text=paste("res.rf <- randomForest::randomForest( formul, data=to.model,",rf.options,")")))

        res.rf;
}

# ----- this function should not be used, it is only a preliminary test --------------
# ----- (instead use:    pred <- predict(res.rf, newdata=to.test)     ) --------------
applyMetacostRf <- function(app,res.rf,to.test,opts) 
{
        response.variable=opts$response.variable;
        formul <- formula(paste(response.variable, "~ ."))   # use all possible input variables
        # we work here with a command text string and eval(...) to allow for the presence or
        # absence of certain options like "mtry" which are not allowed to be NULL:
        rf.options = "ntree=opts$RF.ntree";
        rf.options = paste(rf.options,"sampsize=opts$RF.sampsize",sep=",")
        rf.options = paste(rf.options,"na.action=randomForest::na.roughfix","proximity=F",sep=",")
        if (!is.null(opts$RF.mtry)) rf.options = paste(rf.options,"mtry=opts$RF.mtry",sep=",")

        app$test.votes <- predict(res.rf, newdata=to.test, type="vote")
        Pjx = app$test.votes
        mc.target <- as.factor(apply(Pjx %*% opts$CLS.gainmat, 1,  function(x) { colnames(Pjx)[which(x==max(x))]}))
        
        # re-train
        to.test[,response.variable] <- mc.target
        eval(parse(text=paste("res.rf <- randomForest::randomForest( formul, data=to.test,",rf.options,")")))

        app$test.predict <- predict(res.rf, newdata=to.test)
        app$test.votes <- predict(res.rf, newdata=to.test, type="vote")
        app;
}
