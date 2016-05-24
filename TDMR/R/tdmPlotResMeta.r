######################################################################################
# tdmPlotResMeta:
#
#' Interactive plots of RES data frames and their metamodels.
#'
#' Makes interactive plots for any of the result data frames contained in \code{envT} together with fitted metamodels.
#' \code{tdmPlotResMeta} creates a \code{\link[twiddler]{twiddle}} interface which allows to select
#'    \itemize{
#'      \item tuner:      one of the tuners
#'      \item nExperim:   a knob to select one of the experiments in \code{envT} (only if  \code{envT$tdm$nExperim > 1})
#'      \item reportFunc: one of \{"spotReport3d", "spotReportContour"\}, see \code{\link{SPOT}} 
#'      \item modelFit:   one of \{"spotPredictGausspr", "spotPredictRandomForest", "spotPredictForrester", "spotPredictMlegp"\}, see \code{\link{SPOT}} 
#'      \item skipLastSteps: [0] skip the skipLastSteps sequential steps from the selected RES data frame.  
#'      \item skipIncomplete: [FALSE] if this checkbox is checked (i.e. skipIncomplete==TRUE), then all CONFIGs having *less* than
#'            seq.design.maxRepeats repeats are skipped 
#'      \item nSkip:      [0] skip the nSkip worst CONFIGs
#'      \item y10Exp:     [0] multiply in the RES data frame the column Y (objective function) by the factor 10^y10Exp. This affects the 
#'            coloring of the metamodel in case "spotReport3d": A surface in the Y-range [1e-3, 1e-2] would have only a single color, 
#'            but a surface in the range [10,100] will have a richer color scheme.
#'     } 
#'
#' In case of "spotReport3d", the plot shows the metamodel as colored surface and the true observations of the objective function as 
#' black points. Each point stands for a CONFIG, i.e. it is an aggregation (usually: mean) of all repeats for a design point configuration.
#' The options \code{nskip} and \code{skipIncomplete} allow to suppress certain CONFIGs which may be outliers and thus hinder the view 
#' on the 'interesting' region. Note however, that this suppression also affects the metamodel fit, so use these options with care.
#' 
#' In case of "spotReportContour", the plot shows a contour plot of the metamodel with a black point indicating the best CONFIG found. 
#'
#'   If \code{envT$tunerVal$meta.compare} is TRUE, the quality of the metamodel is evaluated with the following procedure:
#'   A RES data frame \code{newRes} with N new design CONFIGs (not used during metamodel fit) is taken and their true target
#'   value (evaluation of target function) is compared with the metamodel's value. The mean absolute deviation   \cr
#'        sum( | metamodel_i - newCONFIG_i | ) /N     \cr
#'   is printed and returned in \code{envT$tunerVal$meta.newMAD}.    \cr
#'   The RES data frame \code{newRes} is either taken from \code{envT$tunerVal$meta.newRes} (if this element is not NULL) or it is 
#'   taken from another member of envT$resGrid (preferably from a 'lhd' tuner result, because it is evenly distributed). \cr
#'   If \code{envT$tunerVal$meta.compare} is NULL, it is set to FALSE.
#' 
#' @param envT   environment with results as returned from \code{\link{tdmBigLoop}} or as loaded from an appropriate .RData file. 
#'      We use here especially \code{envT$resGrid}, \code{envT$tdm} and \code{envT$tunerVal} 
#'      (the \code{spotConfig} from the last tuning experiment).
#'
#'
#' @note Side Effects:
#'   In case of "spotReport3d", one or several RGL plot windows are created which can be manipulated interactively. 
#'   A certain RGL window \code{n} can be selected with \code{rgl.set(n)}.  
#'   An interactively manipulated RGL window can be saved with \code{rgl.snapshot("myplot.png")}.
#'
#'   In case of "spotPredictMlegp" there is a new element sC$seq.mlegp.min.nugget=1 set in the source code which causes
#'   the MLEGP fit to become much smoother in the presence of a noisy target function (nugget effect).
#'
#' @examples
#'    \dontrun{ 
#'      ##
#'      ## Read previous tuning results 'envT' from demo02sonar/demoSonar.RData 
#'      ## (relative to the TDMR package directory). 
#'      ## Then, tdmPlotResMeta lets you explore interactively the RES data frame(s):
#'      load(paste(find.package("TDMR"), "demo02sonar","demoSonar.RData",sep="/"));
#'      tdmPlotResMeta(envT);
#'    }
#'
#' @seealso   \code{\link{tdmBigLoop}}
#' @author Wolfgang Konen
#' @export
######################################################################################
tdmPlotResMeta <- function(envT) {
  #require(twiddler);     # now specific call with 'twiddler::'

  ################################################################################
  showMeta <- function( theTuner,nExp,confFile,reportFunc,modelFit,
                        nSkip,skipIncomplete,skipLastSteps,y_10Exp,xAxis,yAxis) {
      iGrid = getInd(confFile,nExp,theTuner,envT);
      res = envT$resGrid[[iGrid]];     
      sC <- envT$tunerVal;
      sC$spot.fileMode = FALSE;
      names(sC$alg.roi) <- c("lower","upper","type");
      sC$alg.currentResult = res;  
      
      sC <- configureAxis(sC,xAxis,yAxis,iGrid,envT);
      
      if (sC$configureAxisSuccess==TRUE)  {
        #
        # Preprocess RES data
        #
        sC <- preprocResData(sC,res,nSkip,skipIncomplete,skipLastSteps);
        
        # 
        # Prepare and execute the spot(...,"rep") call
        #
        spotConfFile="NULL";  # "NULL" means: don't read any file, take all values from sC=envT$tunerVal
        sC$seq.predictionModel.func = modelFit; # override the setting in envT$tunerVal
        sC$report.func <- reportFunc;           # override the setting in envT$tunerVal
        sC$seq.modelFit <- NULL;                # needed to force spot to use sC$seq.predictionModel.func
        #sC$seq.mlegp.min.nugget = 1;     # NEW, for "spotPredictMlegp", experimental!
        
        sC$report.main <- paste(theTuner,", nExp=",as.character(nExp),sep="");      # plot title
        sC$alg.currentResult[1] = sC$alg.currentResult[1]*10^y_10Exp;
        sC <- spot(spotConfFile,"rep",NA,sC);

        #
        # Compare in certeain new design points the metamodel's Y with the true Y --> print sC$meta.newMAD, the 
        # mean absolute deviation  | model - true |  of those points
        #
        if (is.null(sC$meta.compare)) sC$meta.compare=FALSE;
        if (sC$meta.compare==TRUE) {
          newTuner = ifelse("lhd" %in% envT$tdm$tuneMethod, "lhd", theTuner);
          if (newTuner!="lhd") warning(sprintf("There is no tuner 'lhd' in envT. Using tuner %s instead (might have uneven distribution).",theTuner));
          if (newTuner==theTuner & envT$tdm$nExperim==1) warning(paste("There is only one result data frame in envT."
                                                                      ,"Design points for model fit and for test are the same!"));
          newNexp = nExp%%envT$tdm$nExperim + 1;
          newRes= envT$resGrid[[iGrid]]; 
          sC <- reportTestMeta(sC,newRes); 
          print(head(cbind(sC$newDes,sC$newY,sC$meta.Y),10));
          cat("sC$meta.newMAD = ",sC$meta.newMAD,"\n");
   	      flush.console();
 	      }
        envT$tunerVal$meta.newMAD <- sC$meta.newMAD; # return meta.newMAD in envT (if envT is environment, not list)
      } # if (sC$configureAxisSuccess==TRUE)
      return(sC);
  } # function showMeta
  
  ################################################################################
  buildTwidCmd <- function(envT) {
    vars = getUnionOfVars(envT);  #rownames(envT$tunerVal$alg.roi);
    tve = length(vars);
    tne = envT$tdm$nExperim;
    tce = length(envT$runList);
    res = envT$resGrid[[1]];
    #maxstep = 0;              # use this setting to disable slider 'skipLastSteps'
    maxstep = max(res$STEP);  # experimental, yields slider 'skipLastSteps'. 
                              # CAUTION: maxstep only valid, if envT$resGrid[[1]] contains as many STEPs as all other RES data frames in envT
    twiddleCmd <- paste("twiddler::twiddle(showMeta(tuner,nExp",sep="");
    if (tne==1)   twiddleCmd <- paste(twiddleCmd,"=1",sep="");
    twiddleCmd <- paste(twiddleCmd,",confFile",sep="");
    if (tce==1)   twiddleCmd <- paste(twiddleCmd,"=envT$runList[1]",sep="");
    twiddleCmd <- paste(twiddleCmd,",reportFunc,modelFit,nSkip,skipIncomplete,skipLastSteps,y_10Exp,xAxis",sep="");
    if (tve==2)   twiddleCmd <- paste(twiddleCmd,"=vars[1]",sep="");
    twiddleCmd <- paste(twiddleCmd,",yAxis",sep="");
    if (tve==2)   twiddleCmd <- paste(twiddleCmd,"=vars[2]",sep="");
    twiddleCmd <- paste(twiddleCmd,"), eval=FALSE",sep="");
                  # eval=FALSE triggers two buttons "EVAL" and "CLOSE"  and inhibits auto-evaluation in twiddle
    tm =  envT$tdm$tuneMethod;
    if (length(tm)>1) {
      twiddleCmd <- paste(twiddleCmd,", tuner=combo(\"",tm[1],"\"",sep="");
      for (i in 2:length(tm)) twiddleCmd <- paste(twiddleCmd,",\"",tm[i],"\"",sep="");
      twiddleCmd <- paste(twiddleCmd,",label=as.character(\"tuner\"))",sep="");
      # the lines above construct a command which looks e.g. like
      #       twiddle(showMeta(tuner), eval=FALSE, tuner=combo("spot","lhd",label=as.character("tuner")))   
      # if tm = c("spot","lhd"). We need to do it this (complicated) way, including eval(parse(...))
      # below, to get the right text strings into the combo boxes.
      # --- It would be nicer, if we could just say ' twiddle(..., tuner=combo(tm,label="tuner"),...)
      # --- if tm is a string vector like c("spot","lhd")
    } else { # i.e. only one tuneMethod
      twiddleCmd <- paste(twiddleCmd,", tuner=combo(\"",tm[1],"\",\"",tm[1],"\",label=as.character(\"tuner\"))",sep="");
      # A bit awkward, but combo needs at least two entries --> we put twice the first and only entry into the combo box.
      # It does also NOT work to set 'tuner <<- "spot"', because then a knob-interface is created 
      # (unclear why, it must have to do s.th. with "spot" being of type string, because a setting 'tuner <<- 2' would 
      # work in twiddler, but is of course meaningless in the context of TDMR)
    }
    if (tne>1) {
      twiddleCmd <- paste(twiddleCmd,", nExp=knob(c(1,",tne,"), res=1, label=\"nExper\")",sep="");
    } 
    
    if (tce>1) {
      twiddleCmd <- paste(twiddleCmd,", confFile=combo(\"",envT$runList[1],"\"",sep="");
      for (i in 2:tce) twiddleCmd <- paste(twiddleCmd,",\"",envT$runList[i],"\"",sep="");
      twiddleCmd <- paste(twiddleCmd,",label=as.character(\"confFile\"))",sep="");
    } 
    
    # if there are more than 2 design variables, build two combo boxes to select two of them: 
    if (tve>2) {
      twiddleCmd <- paste(twiddleCmd,", xAxis=combo(\"<none>\"",sep="");
      for (i in 1:length(vars)) twiddleCmd <- paste(twiddleCmd,",\"",vars[i],"\"",sep="");
      twiddleCmd <- paste(twiddleCmd,",label=as.character(\"xAxis\"))",sep="");
      twiddleCmd <- paste(twiddleCmd,", yAxis=combo(\"<none>\"",sep="");
      for (i in 1:length(vars)) twiddleCmd <- paste(twiddleCmd,",\"",vars[i],"\"",sep="");
      twiddleCmd <- paste(twiddleCmd,",label=as.character(\"yAxis\"))",sep="");
    } else {
      if (length(vars)<2) stop("Need at least two design variables!");
    }
      
    twiddleCmd <- paste(twiddleCmd,", reportFunc = combo(\"spotReport3d\",\"spotReportContour\",label=\"reportFunc\")",sep="");
    twiddleCmd <- paste(twiddleCmd,", modelFit = combo(\"spotPredictGausspr\",\"spotPredictRandomForest\",\"spotPredictForrester\",\"spotPredictMlegp\",label=\"modelFit\")",sep="");
    twiddleCmd <- paste(twiddleCmd,", nSkip=knob(c(0,",min(10,nrow(res)),"), res=1, label=\"nSkip\")",sep="");
    if (maxstep>0)
      twiddleCmd <- paste(twiddleCmd,", skipLastSteps=knob(c(0,",maxstep,"), res=1, label=\"skipLastSteps\")",sep="");
    twiddleCmd <- paste(twiddleCmd,", y_10Exp=knob(c(0,3), res=1, label=\"y_10Exp\")",sep="");
    twiddleCmd <- paste(twiddleCmd,", skipIncomplete=toggle(default=FALSE,label=\"Skip incomplete CONFIGs\")",sep="");
    twiddleCmd <- paste(twiddleCmd,")",sep="");
    #cat(twiddleCmd);
    eval(parse(text=twiddleCmd));
  } # function buildTwidCmd
  
  buildTwidCmd(envT);
    
  cat("tdmPlotResMeta finished\n");
  return
}

################################################################################
# configureAxis: 
#   private helper fct for showMeta
configureAxis <- function(sC,xAxis,yAxis,iGrid,envT) {     
  # NOTE: package tcltk must be in the 'Imports' list in DESCRIPTION (generated via tdmGeneralUtils.r) 
      sC$configureAxisSuccess=FALSE;
      
      # envT$roiGrid is only available since TDMR 0.4.0. If it is there, 
      # it is safer to use its ROI instead of sC$alg.roi
      if (!is.null(envT$roiGrid)) if (length(envT$roiGrid)>0) 
        sC$alg.roi = envT$roiGrid[[iGrid]];
        
      # If there are more than 2 design variables or if at least one of xAxis or yAxis has value "<none>", 
      # then SPOT will open a 2nd twiddle-interface to select two of the design variables
      sC$report.interactive=TRUE;
      if (nrow(sC$alg.roi)==2) {   # only 2 design variables
        sC$report.aIndex=1; 
        sC$report.bIndex=2; 
        sC$report.interactive=FALSE;
      }
      if (xAxis!="<none>" & yAxis!="<none>") {
        vars = setdiff(names(envT$bstGrid[[iGrid]]),c(sC$alg.resultColumn,"COUNT","CONFIG")); # OLD: rownames(sC$alg.roi);
        sC$report.aIndex=which(vars==xAxis); 
        sC$report.bIndex=which(vars==yAxis); 
        if (length(sC$report.aIndex)==0) {
          ttx <- tcltk::tktoplevel(); 
          tcltk::tkpack(l1<-tcltk::tklabel(ttx,text=sprintf("There is no tuning param xAxis=%s in envT$bstGrid[[%d]]",xAxis,iGrid)), 
                 l2<-tcltk::tkbutton(ttx,text="OK",command=function()tcltk::tkdestroy(ttx)));
          sC$configureAxisSuccess=FALSE;
          return(sC);
        }
        if (length(sC$report.bIndex)==0) {
          tty <- tcltk::tktoplevel(); 
          tcltk::tkpack(l1<-tcltk::tklabel(tty,text=sprintf("There is no tuning param yAxis=%s in envT$bstGrid[[%d]]",yAxis,iGrid)), 
                 l2<-tcltk::tkbutton(tty,text="OK",command=function()tcltk::tkdestroy(tty)));
          sC$configureAxisSuccess=FALSE;
          return(sC);
        }
        sC$report.interactive=FALSE;
      }
      sC$configureAxisSuccess=TRUE;
      sC;
}        

################################################################################
# preprocResData: 
#   private helper fct for showMeta
preprocResData <- function(sC,res,nSkip,skipIncomplete,skipLastSteps) {
      mergedData <- spotPrepareData(sC);
      cat("Raw number of CONFIGs:",length(mergedData$CONFIG),"\n");
      if (skipLastSteps>0) {
        maxSTEP = max(0,max(res$STEP)-skipLastSteps);  # Why 'max(0,...)'? - We have to ensure that at least one STEP remains.
        configSkip = which(mergedData$STEP>maxSTEP)
        if (length(configSkip)>0) {
          cat("Skipping",length(configSkip),"CONFIGs in the last steps\n") 
          res <- res[!(res$CONFIG %in% configSkip),]
          sC$alg.currentResult = res;  
          mergedData <- spotPrepareData(sC);
        }
      }
      if (skipIncomplete) {
        configSkip = which(mergedData$count<max(mergedData$count))
        if (length(configSkip)>0) {
          cat("Skipping",length(configSkip),"incomplete CONFIGs \n") 
 	        res <- res[!(res$CONFIG %in% configSkip),]
          sC$alg.currentResult = res;  
          mergedData <- spotPrepareData(sC);
        }
      }
      if (nSkip>0) {
        ind=order(-mergedData$mergedY);
        if (length(ind)>nSkip) {
          cat("Skipping",nSkip,"CONFIGs with highest (worst) Y\n") 
          configSkip <- mergedData$CONFIG[head(ind,nSkip)]
 	        res <- res[!(res$CONFIG %in% configSkip),]
          sC$alg.currentResult = res;  
          mergedData <- spotPrepareData(sC);
        }
      }
      cat("Shown number of CONFIGs:",length(mergedData$CONFIG),"\n");
      flush.console();
      
      sC;     # with sC$alg.currentResult potentially modified
} 	      

################################################################################
# getUnionOfVars: 
#   private helper fct for forming the choices for xAxis, yAxis
getUnionOfVars <- function(envT) {
  vars = NULL;
  for (i in 1:length(envT$bstGrid)) {
    vars=union(vars, setdiff(names(envT$bstGrid[[i]]),c(envT$tunerVal$alg.resultColumn,"COUNT","CONFIG")) );
  }
  vars;
}

################################################################################
# getInd: 
#   private helper fct for picking the right element from envT$resGrid
getInd <- function(confFile,nExp,theTuner,envT) {
  nTuner = length(envT$tdm$tuneMethod);
  indTuner = which(envT$tdm$tuneMethod==theTuner);
  if (length(indTuner)==0) stop(paste("Could not find tuner ",theTuner,"in envT$tdm$tuneMethod"));
  nConf = which(envT$runList==confFile);
  if (length(nConf)==0) stop(paste("Could not find conf file ",confFile,"in envT$runList"));
  if (nExp<1 | nExp>envT$tdm$nExperim) stop(paste("nExp is not in range {1,...,",envT$tdm$nExperim,"}",sep=""));
  ind = indTuner + nTuner*((nExp-1) + envT$tdm$nExperim*(nConf-1));
}

################################################################################
# testMeta: 
#   test function for reportMeta
testMeta <- function(envT) {
        confFile = envT$runList[1];
        res = envT$resGrid[[getInd(confFile,1,"lhd",envT)]];        # RES for metamodel fit 
        newRes = envT$resGrid[[getInd(confFile,1,"lhd",envT)]];     # RES for new (test) design points
    
        sC <- envT$tunerVal;
        names(sC$alg.roi) <- c("lower","upper","type");
        sC$alg.currentResult = res;  

        sC$seq.predictionModel.func = "spotPredictMlegp" # "spotPredictForrester" # "spotPredictGausspr" #
        sC$seq.mlegp.min.nugget = 0.5;     # NEW, experimental!

        sC$spot.fileMode=F;
        sC <- reportTestMeta(sC,newRes); 
        print(head(cbind(sC$newDes,sC$newY,sC$meta.Y),15));
        cat("sC$meta.newMAD = ",sC$meta.newMAD,"\n");
        return(sC);
}

################################################################################
# reportTestMeta: 
#   Helper function for testMeta and for showMeta in tdmPlotResMeta.
#   Given a spotConfig with its RES data frame in spotConfig$alg.currentResult,
#   fit a metamodel to this RES data frame.
#   Simultaneously, calculate the CONFIGs of the new RES data frame newRes and the metamodel's response
# 
#   Return: a modified spotConfig with
#       spotConfig$newDes       the new CONFIGs (from newRes)
#       spotConfig$newY         the true response of the new CONFIGs (aggregated over repeats)
#       spotConfig$meta.Y       the metamodel's response to these new CONFIGs
#       spotConfig$meat.newMAD  the mean abs deviation  | meta.Y - newY |
#       spotConfig$seq.modelFit the metamodel
#
reportTestMeta <- function(spotConfig,newRes) {		
	if (is.null(spotConfig$alg.currentResult)) stop("Need a data frame spotConfig$alg.currentResult to perform!");
	
  sCnew = spotConfig; 
  sCnew$alg.currentResult = newRes;
 	mergedData <- spotPrepareData(sCnew);
  spotConfig$newDes = mergedData$x[,];       # the CONFIGs of newRes become the largeDesign for spotConfig
  spotConfig$newY = mergedData$mergedY[];
	rawB <- spotGetRawDataMatrixB(spotConfig);
	mergedData <- spotPrepareData(spotConfig)
	mergedB <- spotGetMergedDataMatrixB(mergedData, spotConfig);	
	
	spotConfig1 <- eval(call(spotConfig$seq.predictionModel.func
                          , rawB
                          , mergedB           # fit the predictionModel to these CONFIGs
                          , spotConfig$newDes # evaluate the predictionModel at these CONFIGs
                          , spotConfig));
                          
  spotConfig$meta.Y <- spotConfig1$seq.largeDesignY[[1]]; 
  spotConfig$seq.modelFit <- spotConfig1$seq.modelFit;
  spotConfig$meta.newMAD = mean(abs(spotConfig$meta.Y-spotConfig$newY))                              

	return(spotConfig)	
}

#tdmPlotResMeta(envT);

