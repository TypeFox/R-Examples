######################################################################################
#tdmOptsDefaultsSet:
#
#'   Default values for list \code{opts}.  
#' 
#'   Set up and return a list \code{opts} with default settings. The list \code{opts} 
#'   contains all DM-related settings which are needed by main_<TASK>.     
#'   \cr\cr
#'   For better readability, most elements of  \code{opts} are arranged in groups:
#'     \tabular{ll}{
#'      \code{dir.*} \tab  path-related settings  \cr
#'      \code{READ.*} \tab  data-reading-related settings  \cr
#'      \code{TST.*} \tab  resampling-related settings (training, validation and test set, CV)  \cr    
#'      \code{PRE.*} \tab  preprocessing parameters \cr
#'      \code{SRF.*} \tab  several parameters for \code{\link{tdmModSortedRFimport}}   \cr
#'      \code{MOD.*} \tab  general settings for models and model building  \cr
#'      \code{RF.*} \tab  several parameters for model RF (Random Forest)    \cr
#'      \code{SVM.*} \tab  several parameters for model SVM (Support Vector Machines)  \cr
#'      \code{ADA.*} \tab  several parameters for model ADA (AdaBoost)  \cr
#'      \code{CLS.*} \tab  classification-related settings  \cr
#'      \code{GD.*} \tab  settings for the graphic devices  \cr
#'     }
#'
#'   The path-related settings are relative to \code{dir(tdm$mainFile)}, if it is def'd, else relative to the current dir. \cr
#'   Finally, the function \code{\link{tdmOptsDefaultsFill}(opts)} is called to fill in further details, depending on the current 
#'   settings of \code{opts}.
#'
#' @param opts    (optional) the options already set
#' @param path    ["./"] the path string to precede all settings \code{opts$dir.*}. 
#'                If \code{path=="./"} then nothing is preceeded to \code{opts$dir.*}.   
#'
#' @return a list \code{opts}, with defaults set for all options relevant for a DM task, 
#'    containing the following elements
#' 			\item{dir.txt}{[<path>/data] where to find .txt/.csv files} 
#' 			\item{dir.data}{[<path>/data] where to find other data files, including .Rdata  } 
#' 			\item{dir.Rdata}{[<path>/Rdata] -- deprecated, use opts$dir.data -- } 
#' 			\item{dir.output}{[<path>/Output] where to put output files} 
#' 			\item{filename}{["default.txt"] the task data} 
#' 			\item{filetest}{[NULL] the test data, only relevant for READ.TST=T}
#' 			\item{fileMode}{[TRUE] if =T, write opts$EVALFILE=*_train_eval.csv, *_train.csv.SRF.*.RData file and *_train.log file} 
#'   		\item{logFile}{[TRUE] if =T, and if opts$fileMode=T, write log file to *_train.log } 
#' 			\item{data.title}{["Default Data"] title for plots} 
#' 			\item{READ.TXT}{[T] =T: read data from .csv and save as .Rdata, =F: read from .Rdata}                                                   
#' 			\item{READ.NROW}{[-1] read this amount of rows or -1 for 'read all rows'} 
#' 			\item{READ.TST}{[F] =T: read unseen test data from opts$filetest (usually you will do this only for the final model and only with TST.kind="col")} 
#' 			\item{READ.CMD}{["\code{\link{tdmReadCmd}(filename,opts)}"] 
#'                      the command to be passed into \code{\link{tdmReadData}}. It has to contain the placeholder 'filename'. The default 
#'                      in brackets implies 'read.csv(file=paste(opts$dir.txt, filename, sep=""), nrow=opts$READ.NROW)' which includes the 
#'                      further settings header=T, sep="," and dec="."   } 
#' 			\item{READ.INI}{[TRUE] read the task data initially, i.e. prior to tuning, using \code{\link{tdmReadData}} .  
#'                      If =FALSE, the data are read anew in each pass through main_TASK, i.e. in each tuning step (deprecated). } 
#' 			\item{TST.kind}{["rand"] one of the choices from \{"cv","rand","col"\}, see \code{\link{tdmModCreateCVindex}} for details  } 
#' 			\item{TST.COL}{["TST.COL"] name of column with train/test/disregard-flag} 
#' 			\item{TST.NFOLD}{[3] number of CV-folds (only for TST.kind=="cv")} 
#' 			\item{TST.valiFrac}{[0.1] set this fraction of the train-validation data aside for validation (only for TST.kind=="rand")} 
#' 			\item{TST.testFrac}{[0.1] set prior to tuning this fraction of data aside for testing (if tdm$umode=="SP_T" and opts$READ.INI==TRUE)
#'                      or set this fraction of data aside for testing after tuning (if tdm$umode=="RSUB" or =="CV") } 
#'      \item{TST.trnFrac}{[NULL] train set fraction, if NULL then \code{\link{tdmModCreateCVindex}} will set it to 1 - opts$TST.valiFrac. }
#' 			\item{TST.SEED}{[NULL] a seed for the random test set selection (\code{\link{tdmRandomSeed}}) and random validation set selection. 
#'            (\code{\link{tdmClassifyLoop}}). If NULL, use \code{\link{tdmRandomSeed}}. } 
#' 			\item{PRE.PCA}{["none" (default)|"linear"] PCA preprocessing: [don't | do normal PCA (prcomp) ] } 
#' 			\item{PRE.PCA.REPLACE}{[T] =T: replace with the PCA columns the original numerical columns, =F: add the PCA columns  } 
#' 			\item{PRE.PCA.npc}{[0] if >0: add monomials of degree 2 from the first PRE.PCA.npc columns (PCs) (only active, if opts$PRE.PCA!="none")} 
#' 			\item{PRE.SFA}{["none" (default)|"2nd"] SFA preprocessing (see package \code{\link[rSFA]{rSFA}}: [don't | do ormal SFA with 2nd degree expansion ] } 
#' 			\item{PRE.SFA.REPLACE}{[F] =T: replace the original numerical columns with the SFA columns; =F: add the SFA columns } 
#' 			\item{PRE.SFA.npc}{[0] if >0: add monomials of degree 2 from the first PRE.SFA.npc columns (only acitve, if opts$PRE.SFA!="none") } 
#' 			\item{PRE.SFA.PPRANGE}{[11] number of inputs after SFA preprocessing, only those inputs enter into SFA expansion } 
#' 			\item{PRE.SFA.ODIM}{[5] number of SFA output dimensions (slowest signals) to return } 
#' 			\item{PRE.SFA.doPB}{[T] =F|T: don't | do parametric bootstrap for SFA in case of marginal training data } 
#' 			\item{PRE.SFA.fctPB}{[sfaPBootstrap] the function to call in case of parametric bootstrap, see \code{\link[rSFA]{sfaPBootstrap}} 
#'                      in package \code{\link[rSFA]{rSFA}} for its interface description } 
#'      \item{PRE.allNonVali}{[F] if =T, then use all non-validation data in the training-validation set for PCA or SFA preprocessing.
#'                      If =F, use only the training set for PCA or SFA processing (only relevant if opts$PRE.PCA!="none" or opts$PRE.SFA!="none").  } 
#' 			\item{PRE.Xpgroup}{[0.99] bind the fraction 1-PRE.Xpgroup in column OTHER (see \code{\link{tdmPreGroupLevels}})  } 
#' 			\item{PRE.MaxLevel}{[32] bind the N-32+1 least frequent cases in column OTHER (see \code{\link{tdmPreGroupLevels}})  } 
#' 			\item{SRF.kind}{["xperc" (default) |"ndrop" |"nkeep" |"none" ] the method used for feature selection, see \code{\link{tdmModSortedRFimport}}  } 
#'      \item{SRF.ndrop}{   [0] how many variables to drop (only relevant if SRF.kind=="ndrop")  }
#'      \item{SRF.nkeep}{ [NULL] how many variables to keep, NULL="keep all" (only relevant if SRF.kind=="nkeep") }
#'      \item{SRF.XPerc}{  [0.95] if >=0, keep that importance percentage, starting with the most important variables (if SRF.kind=="xperc")  }
#'      \item{SRF.calc}{   [T] =T: calculate importance & save on SRF.file, =F: load from srfFile
#'                      (srfFile = Output/<confFile>.SRF.Rdata) }
#'      \item{SRF.ntree}{  [50] number of RF trees }
#'      \item{SRF.samp}{    sampsize for RF in importance estimation. 
#'                          See RF.samp for further info on sampsize. }
#'      \item{SRF.verbose}{ [2] }
#'      \item{SRF.maxS}{    [40] how many variables to show in plot }
#'      \item{SRF.minlsi}{  [1] a lower bound for the length of SRF$input.variables  }
#'      \item{SRF.method}{ ["RFimp"] }
#'      \item{SRF.scale}{  [TRUE] option 'scale' for call importance() in \code{\link{tdmModSortedRFimport}}   }
#' 			\item{MOD.SEED}{[NULL] a seed for the random model initialization (if model is non-deterministic). If NULL, use \code{\link{tdmRandomSeed}}. } 
#' 			\item{MOD.method}{["RF" (default) |"MC.RF" |"SVM" |"NB" ]: use [RF | MetaCost-RF | SVM | Naive Bayes ] in \code{\link{tdmClassify}}  \cr
#'                      ["RF" (default) |"SVM" |"LM" ]: use [RF | SVM | linear model ] in \code{\link{tdmRegress}}  } 
#' 			\item{RF.ntree}{[500] } 
#' 			\item{RF.samp}{[NULL] sampsize for RF in model training. If RF.samp is a scalar, then it specifies the 
#'   		                      total size of the sample. For classification, it can also be a vector of length n.class 
#'     	                      (= # of levels in response variable), then it specifies the size of each strata. The sum 
#'                             of the vector is the total sample size. If NULL, RF.samp will be replaced by 3000 later
#'                             in tdmModAdjustSampsize*.  } 
#' 			\item{RF.mtry}{[NULL] } 
#' 			\item{RF.nodesize}{[1] } 
#' 			\item{RF.OOB}{[TRUE] if =T, return OOB-training set error as tuning measure; if =F, return validation set error } 
#' 			\item{RF.p.all}{[FALSE]  }
#' 			\item{SVM.kernel}{[3] =1: linear, =2: polynomial, =3: RBF, =4: sigmoid}
#' 			\item{SVM.epsilon}{[0.005] needed only for regression}
#' 			\item{SVM.gamma}{[0.005] }
#' 			\item{SVM.coef0}{[0.0] (needed only for opts$SVM.kernel=="polynomial" or =="sigmoid")}
#' 			\item{SVM.degree}{[3] (needed only for opts$SVM.kernel=="polynomial")}
#' 			\item{SVM.tolerance}{[0.008] }
#' 			\item{ADA.coeflearn}{[1] =1: "Breiman", =2: "Freund", =3: "Zhu" as value for boosting(...,coeflearn,...) (AdaBoost)  }
#' 			\item{ADA.mfinal}{[10] number of trees in AdaBoost = mfinal boosting(...,mfinal,...)  }
#' 			\item{ADA.rpart.minsplit}{[20] minimum number of observations in a node in order for a split to be attempted  }
#' 			\item{CLS.cutoff}{[NULL] vote fractions for the classes (vector of length n.class = # of levels in response variable). 
#'   		                The class i with maximum ratio (\% votes)/CLS.cutoff[i] wins. 
#'                      If NULL, then each class gets the cutoff 1/n.class (i.e. majority vote wins). 
#'                      The smaller CLS.cutoff[i], the more likely class i will win. }
#'      \item{CLS.CLASSWT}{ [NULL] class weights for the n.class classes, e.g. \cr
#'                                   c(A=10,B=20) for a 2-class problem with classes A and B         \cr
#'                      (the higher, the more costly is a misclassification of that real class). It should be a named vector with the same
#'                      length and names as the levels of the response variable. If no names are given, the levels of the response variables 
#'                      in lexicographical order will be attached in \code{\link{tdmClassify}}. CLS.CLASSWT=NULL for no weights.  }
#' 			\item{CLS.gainmat}{[NULL] (n.class x n.class) gain matrix. If NULL, CLS.gainmat will be set to unit matrix in \code{\link{tdmClassify}} }
#' 			\item{rgain.type}{["rgain" (default) |"meanCA" |"minCA" ] in case of \code{\link{tdmClassify}}: For classification, the measure 
#'                      \code{Rgain} returned from \code{\link{tdmClassifyLoop}} in \code{result$R_*} is
#'                      [relative gain (i.e. gain/gainmax) | mean class accuracy | minimum class accuracy ]. 
#'                      The goal is to maximize  \code{Rgain}. \cr 
#'                      For binary classification there are the additional measures [ "arROC" | "arLIFT" | "arPRE" ], see 
#'                      \code{\link{tdmModConfmat}}. \cr
#'                      For regression, the goal is to minimize \code{result$R_*} returned from \code{\link{tdmRegress}}. In this case, possible values are 
#'                      \code{rgain.type} = ["rmae" (default) |"rmse"  | "made" ] which stands for 
#'                      [ relative mean absolute error | root mean squared error | mean absolute deviation ].  } 
#' 			\item{ncopies}{[0] if >0, activate \code{\link{tdmParaBootstrap}} in \code{\link{tdmClassify}}  } 
#'      \item{fct.postproc}{[NULL] name of a function with signature \code{(pred, dframe, opts)} where \code{pred} is the prediction of the model on the 
#'                      data frame \code{dframe} and \code{opts} is this list. This function may do some postprocessing on \code{pred}  and
#'                      it returns a (potentially modified) \code{pred}. This function will be called in \code{\link{tdmClassify}} if it is not \code{NULL}.  }
#' 			\item{GD.DEVICE}{["win"] ="win": all graphics to (several) windows (\code{windows} or \code{X11} in package \code{grDevices}) \cr
#'   		                ="rstudio": same as "win", but all graphics go to the RStudio device \cr
#'                      ="pdf": all graphics to one multi-page PDF \cr
#'                      ="png": all graphics in separate PNG files in \code{opts$GD.PNGDIR} \cr
#'                      ="non": no graphics at all \cr
#'                      This concerns the TDMR graphics, not the SPOT (or other tuner) graphics. 
#'                      If running R from RStudio (if there is a device with name "RStudioGD")
#'                      then the default "win" is changed to "rstudio" automatically.    } 
#' 			\item{GD.RESTART}{[T] =T: restart the graphics device (i.e. close all 'old' windows or re-open 
#'                      multi-page pdf) in each call to \code{\link{tdmClassify}} or \code{\link{tdmRegress}}, resp. \cr
#'                      =F: leave all windows open (suitable for calls from SPOT) or write more pages in same pdf. } 
#' 			\item{GD.CLOSE}{[T] =T: close graphics device "png", "pdf" at the end of main_*.r (suitable for main_*.r solo) or \cr
#'                      =F: do not close (suitable for call from tdmStartSpot, where all windows should remain open)  } 
#' 			\item{NRUN}{[2] how many runs with different train & test samples  - or - how many CV-runs, if \code{opts$TST.kind}="cv"  } 
#'   		\item{APPLY_TIME}{[FALSE]   } 
#'   		\item{test2.show}{[FALSE]   } 
#'     	\item{test2.string}{["default cutoff"]   } 
#' 			\item{VERBOSE}{[2] =2: print much output, =1: less, =0: none} 
#'
#' @note  The variables opts$PRE.PCA.numericV and opts$PRE.SFA.numericV (string vectors of numeric input columns to be used for PCA or SFA) 
#'      are not set by \code{\link{tdmOptsDefaultsSet}} or \code{\link{tdmOptsDefaultsFill}}. Either they are supplied by the user or, 
#'      if NULL, TDMR will set them to \code{input.variables} in \code{\link{tdmClassifyLoop}}, assuming that all columns are numeric. 
#----- now obsolete: ----
#      If PCA is done, its output \code{pca$numeric.variables} will overwrite \code{opts$PRE.SFA.numericV} (because the numeric variables 
#      after PCA become the input for SFA).
#'
#' @seealso  \code{\link{tdmOptsDefaultsFill}} \code{\link{tdmDefaultsFill}}
#' @author Wolfgang Konen, FHK, 2011 - 2013
#' @export
######################################################################################
tdmOptsDefaultsSet <- function(opts=NULL, path="./") {
  if (is.null(opts)) {
      opts = list()

      opts$dir.data <- paste("./", "data/", sep="")
      opts$dir.txt  <- paste("./", "data/", sep="")
      opts$dir.Rdata <- paste("./", "Rdata/", sep="")
      opts$dir.output <- paste("./", "Output/", sep="")
      opts$filename = "default.txt"
      opts$data.title <- "Default Data"
      opts$fileMode = TRUE    # =T: write opts$EVALFILE=*_train_eval.csv, *_train.csv.SRF.*.RData file and *_train.log file
      opts$logFile = TRUE
      
      opts$READ.TXT = TRUE    # =T: read data from .csv and save as .Rdata, =F: read from .Rdata
      opts$READ.NROW = -1     # [-1] read this amount of rows or -1 for 'read all rows' 
      opts$READ.TST = FALSE   # =T: read unseen test data (do this only for the final model and only with TST.kind="col")
                              # and fill column dset[,opts$TST.COL] accordingly (set it to 1 for those test records)
                              # =F: set a part of the train data aside as test data (as prescribed by TST.kind)
      opts$READ.CMD = "tdmReadCmd(filename,opts)"; 
      opts$READ.INI = TRUE;   # read in the task data initially, i.e. prior to tuning
      opts$TST.kind <- "rand" # ["cv"|"rand"|"col"] see tdmModCreateCVindex in tdmModelingUtils.r
      opts$TST.COL ="TST.COL";# column with train/test/disregard-flag
      opts$TST.NFOLD =  3     # number of CV-folds (only for TST.kind=="cv")
      opts$TST.valiFrac = 0.10    # set this fraction of data aside for validation (only for TST.kind=="rand")
      opts$TST.SEED = NULL    # [NULL] a seed for the random test set selection
      opts$MOD.SEED = NULL    # [NULL] a seed for the random model initialization (if model is non-deterministic)
      
      opts$PRE.PCA = "none"   # ["none"|"linear"|"kernel"] PCA preprocessing: [don't | normal pca (prcomp) | kernel pca (kernlab) ]
                              # (so far we have problems with kernlab in this app --> currently commented out in tdmPreprocUtils.r)
      opts$PRE.knum = 0;      # [0] if >0 and if PRE.PCA="kernel", take only a subset of PRE.knum records from dset                         
      opts$PRE.PCA.REPLACE=T; # [T] =T: replace the original numerical columns with the PCA columns; =F: add the PCA columns
      opts$PRE.PCA.npc <- 0;  # [0] if >0: add monomials of degree 2 for the first PRE.PCA.npc columns (PCs)
      opts$PRE.SFA = "none"   # ["none"|"2nd"] SFA preprocessing: [don't | normal SFA with 2nd degree expansion ]
      opts$PRE.SFA.REPLACE=F; # [F] =T: replace the original numerical columns with the SFA columns; =F: add the SFA columns
      opts$PRE.SFA.npc <- 0;  # [0] if >0: add monomials of degree 2 for the first PRE.SFA.npc columns 
      opts$PRE.SFA.PPRANGE=11; 
      opts$PRE.SFA.ODIM=5;
      opts$PRE.SFA.doPB=TRUE;
      opts$PRE.SFA.fctPB=rSFA::sfaPBootstrap;
      opts$PRE.allNonVali=F;  # [F] =T: use all non-validation data for PCA or SFA preprocessing (only relevant if opts$PRE.PCA!="none" or opts$PRE.SFA!="none")
      opts$PRE.Xpgroup=0.99;
      opts$PRE.MaxLevel=32;
    
      opts$SRF.kind = "xperc" # ["xperc"|"ndrop"|"nkeep"|"none"] see tdmModSortedRFimport in tdmModelingUtils.r
      opts$SRF.ndrop = 10;    # 0..n: how many variables (those with lowest importance) to drop (only for SRF.kind=="ndrop")
      opts$SRF.ntree=50;      #
      opts$SRF.verbose=2;     #
      opts$SRF.maxS=40;       #
      opts$SRF.minlsi=1;      #
      opts$SRF.XPerc = 0.95;  # 0.0..1.0: how much of the overall importance to keep
      opts$SRF.scale = TRUE;  # option 'scale' for call importance() in tdmModSortedRFimport 
      opts$SRF.calc = TRUE    # =T: calculate SRF variable ranking
                              # =F: reload previously calculated and stored variable ranking
                                    
      opts$MOD.method="RF";   # ["RF"|"MC.RF"|"SVM"|"NB"]: use [RF| MetaCost-RF| SVM| Naive Bayes] in tdmClassify
                              # ["RF"|"SVM"|"LM"]: use [RF| SVM| linear model] in tdmRegress

      opts$RF.ntree = 500
      opts$RF.samp = NULL # 1000
      opts$RF.mtry = NULL
      opts$RF.nodesize = 1
      opts$RF.OOB = TRUE;     # if =T, return OOB-training set error as tuning measure; if =F, return test set error
      opts$RF.p.all=FALSE;
      opts$SVM.kernel=3;      # [3] =1: linear, =2: polynomial, =3: RBF, =4: sigmoid
      opts$SVM.cost=1.0;
      opts$SVM.epsilon=0.005; # needed only for regression
      opts$SVM.gamma=0.005;
      opts$SVM.coef0=0.0;     # needed only for opts$SVM.kernel=="polynomial" or =="sigmoid"
      opts$SVM.tolerance=0.008;
      opts$SVM.degree=3;      # needed only for opts$SVM.kernel=="polynomial"
      opts$ADA.coeflearn=1;   # [1] =1: "Breiman", =2: "Freund", =3: "Zhu" as value for boosting(...,coeflearn,...) (AdaBoost)
      opts$ADA.mfinal=10;     # [10] number of trees in AdaBoost = mfinal boosting(...,mfinal,...) 
      opts$ADA.rpart.minsplit=20; # minimum number of observations in a node in order for a split to be attempted  
      opts$CLS.cutoff = NULL  # [NULL] vote fractions for the n.class classes. The class i with
                              # maximum ratio (% votes)/RF.cutoff[i] wins. If NULL, then each
                              # class gets the cutoff 1/n.class (i.e. majority vote wins).
      opts$CLS.CLASSWT = NULL # class weights for the n.class classes, 
      opts$CLS.gainmat <- NULL# if [NULL], opts$CLS.gainmat will be set to unit matrix in \code{\link{tdmClassify}}
      opts$rgain.type="rgain" # ["rgain" (def.)|"meanCA"|"minCA"] for tdmClassify: the measure result$R_* will contain 
                              # relative gain (i.e. gain/gainmax), mean or minimum class accuracy. 
                              # In each case the goal is to maximize the measure.
                              # ["rmae" (default) |"rmse" ] for regression (tdmRegress). Here the goal is to minimize result$R_*.
      opts$ncopies = 0;       # if >0, activate tdmParaBootstrap in tdmClassify                       

      opts$DO.POSTPROC = FALSE; # --- deprecated, use opts$fct.postproc=NULL ---
      opts$fct.postproc=NULL;
      opts$DO.GRAPHICS=T      # --- deprecated, use opts$GD.DEVICE="non" ---
      opts$GD.DEVICE="win"    # ="pdf": all graphics to one multi-page PDF
                              # ="win": all graphics to (several) windows (X11)
                              # ="rstudio": all graphics to the *one* RStudio graphics device
                              # ="non": no graphics at all
                              # This concerns the TDMR graphics, not the SPOT graphics
      opts$GD.RESTART=TRUE    # [T] =T: restart the graphics device (i.e. close all 'old' windows
                              # or re-open multi-page pdf) in each call to tdmClassfiy or tdmRegress, resp.
                              # =F: leave all windows open (suitable for calls from SPOT)
                              # or write more pages in same pdf
      opts$GD.CLOSE=TRUE      # [T] =T: close graphics device "png", "pdf" at the end of main_*.r
                              # (suitable for main_*.r solo) or =F: do not close (suitable for
                              # call from tdmStartSpot, where all windows should remain open)
      opts$NRUN =  2          # how many runs with different train & test samples  - or -
                              # how many CV-runs, if opts$TST.kind=="cv"
      opts$rep=1;             # the number of the repeat (1,...,spotConfig$max.repeats) in case of repeated evocations from tuner
                              # (needed only internally as private storage for the RNG, see tdmClassify.r
      opts$APPLY_TIME=FALSE;                              
      opts$test2.show <- FALSE;
      opts$test2.string <- "default cutoff";
      opts$VERBOSE=2;

      opts$srf = list();  # private storage to pass results back/forth from tdmModSortedRFimport
                          # opts$srf[["resp"]] holds a data frame with results from/for SRF with response variable "resp".
  }  #if (is.null(opts))
  
# DON'T!! - saving the absolute path in elements of opts gives problems when envT is saved on .RData and later loaded
#           (perhaps on another computer where the absolute path does not exist) --> this made the example in unbiasedRun.r
#           not runnable on CRAN computers, because demoSonar.RData had the path 
#           dir.data = "C:/Users/wolfgang/Documents/R/win-library/2.15/TDMR/demo02sonar/./data/"
#
#  if (!is.null(path)) 
#    if (!(path=="./" | length(grep(path,opts$dir.data))==1)) {
#        #preceede the following strings with 'path' *only*, if path is not "./" and  
#        #if opts$dir.data not contains 'path' (first pass through tdmOptsDefaulsSet)
#        opts$dir.data <- paste(path,opts$dir.data,sep="")
#        opts$dir.txt  <- paste(path,opts$dir.txt,sep="")
#        opts$dir.Rdata <- paste(path,opts$dir.Rdata,sep="")
#        opts$dir.output <- paste(path,opts$dir.output,sep="")
#    }
    

  class(opts) <- c("tdmOpts","TDM")     

  tdmOptsDefaultsFill(opts);
}

######################################################################################
# tdmOptsDefaultsFill:
#
#' Fill the current \code{opts}. 
#'
#'   There is no need to call this function directly, use \code{\link{tdmOptsDefaultsSet}}(opts) 
#'   instead, which calls tdmOptsDefaultsFill in the end. \cr
#'   \code{tdmOptsDefaultsFill} fills the current \code{opts} with further default 
#'   parameters if they are not yet defined. The defaults may depend on previously 
#'   defined elements of \code{opts} (e.g. \code{opts$filename}). 
#'
#'   What is the difference between \code{\link{tdmOptsDefaultsSet}} and \code{\link{tdmOptsDefaultsFill}}? 
#'   \code{tdmOptsDefaultsSet} is for all parameters that do NOT depend on previously def'd elements of \code{opts}.
#'   \code{tdmOptsDefaultsFill} is used to fill in further \code{opts} elements, if not yet defined, depending on 
#'   previous settings (e. g. opts$LOGFILE is derived from opts$filename).
#'
#' @param opts    the options 
#  -- deprecated -- @param suffix  the suffix of \code{opts$filename}. If NULL, take opts$filesuffix (which, if also NULL, is  
#  --            --            inferred from opts$filename)
#' @return \code{opts},  the extended options, where additional elements, if they are not yet def'd,  are set as: 
#' 			\item{TST.COL}{["TST.COL"] } 
#' 			\item{PDFFILE}{["*_pic.pdf"] file for multipage graphics in case \code{opts$GD.DEVICE}="pdf" } 
#' 			\item{GD.PNGDIR}{["PNG*"] directory for .png files in case \code{opts$GD.DEVICE}="png" } 
#' 			\item{LOGFILE}{["*.log"] where to log the output } 
#' 			\item{EVALFILE}{["*_eval.csv"] file with evaluation results allEVAL } 
#' 			\item{SRF.samp}{sample size for SRF, derived from \code{RF.samp} } 
# 			\item{SRF.cutoff}{[opts$CLS.cutoff] } 
#'      	\item{rgain.string}{ one out of c("RGain","MeanCA","MinCA","RMAE","RMSE","MADE","AreaROC","AreaLift","AreaPrRe"), 
#'                          depending on \code{opts$rgain.type} }
#'
#' Here * is the stripped part of \code{opts$filename} (w/o suffix).
#'
#' All files and directories in the above settings are relative to dir  \code{opts$dir.output}.
#'
#' @seealso  \code{\link{tdmOptsDefaultsSet}}
#' @author Wolfgang Konen, FHK, Mar'2011 - May'2012
#' @keywords internal
#' @export
######################################################################################
tdmOptsDefaultsFill <- function(opts) {  #,suffix=NULL) {
    if (is.null(opts)) opts = tdmOptsDefaultsSet();
    if (class(opts)[1] != "tdmOpts")  stop("Class of object opts is not tdmClass");

    filename = opts$filename; 
    ssp = strsplit(filename,".",fixed=TRUE);
    if (length(ssp[[1]])==1) {                  
      filename = paste(filename,".suf",sep=""); # filename has no suffix --> add a dummy suffix ".suf"
      ssp = strsplit(filename,".",fixed=TRUE);
    }
    filesuffix = paste(".",tail(unlist(ssp),1),sep="");

    if (is.null(opts$TST.COL)) opts$TST.COL="TST.COL";
    # bug fix 05/12: removed the is.null-check from  opts$PDFFILE, $GD.PNGDIR, $LOGFILE, $EVALFILE
    # (otherwise it might happen that opts$LOGFILE = "default.log"):
    opts$PDFFILE=sub(filesuffix,"_pic.pdf",filename)
    opts$GD.PNGDIR=paste("PNG",sub(filesuffix,"",filename),"/",sep="");
    opts$LOGFILE=sub(filesuffix,".log",filename)
    opts$EVALFILE=sub(filesuffix,"_eval.csv",filename)      # contains evaluation results allEVAL
    
    if (is.null(opts$rgain.type)) opts$rgain.type="rgain";
    #if (is.null(opts$rgain.string)) {
      rgainTypeVals = c("rgain","meanCA","minCA","rmae","rmse","made","arROC","arLIFT","arPRE");
      rgainStringVals=c("RGain","MeanCA","MinCA","RMAE","RMSE","MADE","AreaROC","AreaLift","AreaPrRe");
      ind = which(opts$rgain.type==rgainTypeVals);
      if (length(ind)==0) stop(sprintf("Could not find opts$rgain.type=%s in c(%s)"
                                        ,opts$rgain.type,"'rgain','meanCA','minCA','rmae','rmse','made','arROC','arLIFT','arPRE'"));
      opts$rgain.string = rgainStringVals[ind];
    #}

    if (!is.null(opts$DO.GRAPHICS) && is.null(opts$GD.DEVICE))
      if (opts$DO.GRAPHICS==F) opts$GD.DEVICE="non";        # DO.GRAPHICS is now (01/2011) deprecated

    if (!is.null(opts$GRAPHDEV) && is.null(opts$GD.DEVICE))
      opts$GD.DEVICE=opts$GRAPHDEV;                         # opts$GRAPHDEV is now (12/2011) deprecated

    #if("rstudio" %in% names(sessionInfo()$loadedOnly)) {      # are we running R from RStudio?
    #  if (opts$GD.DEVICE=="win") opts$GD.DEVICE="rstudio";    # --- if yes, we have only one graphic device
    #}    
    # unfortunately the above if-clause does no longer (with RStudio0.97, R3.0.2) detect  
    # whether we are running R from RStudio. Use instead:
    if(Sys.getenv("RSTUDIO") =="1")  {
      if (opts$GD.DEVICE=="win") opts$GD.DEVICE="rstudio";
    }
    
    # code which was previously in tdmClassify. Now we put it here and call tdmOptsDefaultsFill from tdmClassify
    # (cleaner code, less places where opts-values are set)
    #
    if (is.null(opts$MOD.method)) opts$MOD.method="RF";
    if (is.null(opts$RF.p.all)) opts$RF.p.all=FALSE;
    if (is.null(opts$RF.OOB)) opts$RF.OOB=TRUE;
    if (is.null(opts$APPLY_TIME)) opts$APPLY_TIME=FALSE;
    if (is.null(opts$DO.POSTPROC)) opts$DO.POSTPROC=FALSE;    # --- deprecated, use opts$fct.postproc=NULL ---
    if (is.null(opts$GD.RESTART)) opts$GD.RESTART=TRUE;  
    if (is.null(opts$VERBOSE)) opts$VERBOSE=2;
    if (is.null(opts$test2.show)) opts$test2.show <- FALSE;
    if (is.null(opts$test2.string)) opts$test2.string <- "default cutoff";
    
    # code which was previously in tdmModSortedRFimport. Now we put it here and call tdmOptsDefaultsFill from tdmModSortedRFimport
    # (cleaner code, less places where opts-values are set)
    #
    if (is.null(opts$SRF.samp)) opts$SRF.samp=min(opts$RF.samp,3000);  # new 06/2011
    if (is.null(opts$SRF.cutoff)) opts$SRF.cutoff=opts$CLS.cutoff;     # new 01/2011, might change the default behaviour
    if (is.null(opts$SRF.kind)) opts$SRF.kind="xperc";
    if (is.null(opts$SRF.XPerc)) opts$SRF.XPerc=0.95;
    if (is.null(opts$SRF.calc)) opts$SRF.calc=TRUE;
    if (is.null(opts$SRF.ntree)) opts$SRF.ntree=50;
    if (is.null(opts$SRF.verbose)) opts$SRF.verbose=2;
    if (is.null(opts$SRF.maxS)) opts$SRF.maxS=40;
    if (is.null(opts$SRF.minlsi)) opts$SRF.minlsi=1;
    if (is.null(opts$SRF.method)) opts$SRF.method="RFimp";
    if (is.null(opts$SRF.scale)) opts$SRF.scale = TRUE;  # option 'scale' for call importance() in tdmModSortedRFimport 

    opts;
}
