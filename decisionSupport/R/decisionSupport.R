#
# file: decisionSupport.R
#
# This file is part of the R-package decisionSupport
#
# Authors:
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF)
#	http://www.worldagroforestry.org
#
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include estimate_read_csv_old.R
#' @include individualEvpiSimulation.R
NULL
##############################################################################################
# decisionSupport(inputFilePath, outputPath, welfareFunction, numberOfModelRuns,
#                     randomMethod="calculate",	functionSyntax="data.frameNames",
#                     write_table=TRUE, 
#                     oldInputStandard=FALSE)
##############################################################################################
#' Welfare Decision and Value of Information Analysis wrapper function.
#'
#' This function performs a Welfare Decision Analysis via a Monte Carlo simulation from input files
#' and analyses the value of different information about the input variables. This value of
#' information analysis can be done via combined PLSR - VIP analysis or via IndividualEVPI
#' calculation. Results are saved as plots and tables.
#' @param inputFilePath Path to input csv file, which gives the input \code{\link{estimate}}.
#' @param outputPath Path where the result plots and tables are saved.
#' @param welfareFunction The welfare function.
#' @param numberOfModelRuns The number of running the welfare model for the underlying Monte Carlo 
#' simulation.
#' @param write_table \code{logical}: If the full Monte Carlo simulation results and PLSR results should be
#'  written to file.
#' @param plsrVipAnalysis \code{logical}: If PLSR-VIP analysis shall be performed.
#' @param individualEvpiNames \code{character vector}: names of variables, which for the 
#'   IndividualEVPI shall be obtained via Monte Carlo simulation. If \code{=NULL} (the default), no 
#'   IndividualEVPI is calculated; if \code{="all"}, the IndividualEVPI is calculated for all
#'   variables. \dfn{Note:} depending on \code{numberOfModelRuns} and the complexity of
#'   \code{welfare} this might take a long time.
#' @param sortEvpiAlong \code{character}: result name along which the summary of the IndividualEVPI
#'    shall be sorted. Only relevant if \code{sortEvpiAlong!=NULL}. 
#' @param oldInputStandard \code{logical}: If the old input standard should be used
#' 	(\code{\link{estimate_read_csv_old}}).
# @param verbosity \code{integer}: if \code{0} the function is silent; the larger the value the
#   more verbose is output information.
#' @inheritParams welfareDecisionAnalysis
#' @details
#'  This function integrates the most important features of  
#' 	\link[=decisionSupport-package]{this package} into a single function. It is wrapped arround the functions 
#' 	\code{\link{welfareDecisionAnalysis}}, \code{\link{plsr.mcSimulation}}, 
#' 	\code{\link[chillR:VIP]{VIP}} and \code{\link{individualEvpiSimulation}}.
#'   \subsection{Combined PLSR - VIP Analysis}{
#'   The combined Partial Least Squares Regression (PLSR) and Variables Importance in Projection 
#'   (VIP) analysis is implemented via: \code{\link{plsr.mcSimulation}} and 
#'   \code{\link[chillR:VIP]{VIP}}.
#'   }
#'   \subsection{IndividualEVPI Calculation}{
#'   Implementation: \code{\link{individualEvpiSimulation}}
#'   }
#' @seealso \code{\link{mcSimulation}}, \code{\link{estimate}}, \code{\link{estimate_read_csv}}, 
#' 	\code{\link{plsr.mcSimulation}}, \code{\link[chillR:VIP]{VIP}}, 
#' 	\code{\link{welfareDecisionAnalysis}}, \code{\link{individualEvpiSimulation}}, 
#' 	\code{\link{decisionSupport-package}}
#' @export
decisionSupport <- function(inputFilePath, outputPath, welfareFunction, numberOfModelRuns,
                            randomMethod="calculate",	
                            functionSyntax="data.frameNames",
                            relativeTolerance=0.05,
                            write_table=TRUE, 
                            plsrVipAnalysis=TRUE,
                            individualEvpiNames=NULL,
                            sortEvpiAlong=if(individualEvpiNames) individualEvpiNames[[1]] else NULL,
                            oldInputStandard=FALSE,
                            verbosity=1){
  # Read estimate from file:
  if(!oldInputStandard){
    #	print("newInputStandard")
    estimateObject<-estimate_read_csv(fileName=inputFilePath)
  }else{
    #	print("oldInputStandard")
    estimateObject<-estimate_read_csv_old(fileName=inputFilePath)
  }
  if(verbosity > 0)
    cat("Estimate read from file: ",inputFilePath, "\n")
  if(verbosity > 1)
    print(estimateObject)
  # Create the output directory if necessary
  if ( !file.exists(outputPath) )
      dir.create(outputPath, recursive=TRUE)
  # Run Monte Carlo simulation for the Welfare Decision Analysis:
  if(verbosity > 0)
    cat("Performing Monte Carlo Simulation for the Welfare Decision Analysis:\n")
  welfareDecisionResults<-welfareDecisionAnalysis(estimate=estimateObject,
                                                  welfare=welfareFunction,
                                                  numberOfModelRuns=numberOfModelRuns,
                                                  randomMethod=randomMethod,
                                                  functionSyntax=functionSyntax,
                                                  relativeTolerance=relativeTolerance,
                                                  verbosity=verbosity)
  if(verbosity > 0)
    cat("Monte Carlo Simulation for Welfare Decision Analysis done.\n")
  if(verbosity > 1){
    print(summary(welfareDecisionResults$mcResult))
    print(summary(welfareDecisionResults))
  }
  if (write_table){
    mcSimulationResultsFilePath<-file.path(outputPath, "mcSimulationResults.csv")
    write.csv(data.frame(welfareDecisionResults$mcResult$x,welfareDecisionResults$mcResult$y),
              mcSimulationResultsFilePath)
    if (verbosity > 0)
      cat("Monte Carlo results written into file: ", mcSimulationResultsFilePath, "\n")
  }
  # Write histogram of results to png files:
  for(i in names(welfareDecisionResults$mcResult$y)) {
    png(file.path(outputPath, paste(i, "_distribution.png",sep="")), width=1000, height=500)
    par(mar=c(5.1,5.1,4.1,2.1))
    hist(welfareDecisionResults$mcResult, lwd=3, cex.lab=2 ,cex.axis=2, prob=FALSE, resultName=i)
    dev.off()
  }
  # Write the summary of the resulting distributions to file:
  mcSummary<-summary(welfareDecisionResults$mcResult, digits=2)
  mcSimulationSummaryFilePath<-file.path(outputPath,"mcSummary.csv")
  write.csv(mcSummary$summary,mcSimulationSummaryFilePath)
  if (verbosity > 0)
    cat("Monte Carlo summary written into file: ", mcSimulationSummaryFilePath, "\n")
  # Write the summary of the Welfare Decision Analysis to file:
  welfareDecisionSummary<-summary(welfareDecisionResults, digits=2)
  welfareDecisionSummaryFilePath<-file.path(outputPath,"welfareDecisionSummary.csv")
  write.csv(welfareDecisionSummary$summary,welfareDecisionSummaryFilePath)
  if (verbosity > 0)
    cat("Welfare Decision Analysis summary written into file: ", welfareDecisionSummaryFilePath, "\n")
  # Partial lest squares analysis:
  if(plsrVipAnalysis){
    if (!requireNamespace("chillR", quietly = TRUE)) 
      stop("Package \"chillR\" needed. Please install it.",
           call. = FALSE)
    for(i in names(welfareDecisionResults$mcResult$y)) {
      # Run PLSR on the Monte Carlo result:
      plsrResult<-plsr.mcSimulation(welfareDecisionResults$mcResult, resultName=i)
      # Run VIP analysis on the PLSR result:
      vipResult<-chillR::VIP(plsrResult)["Comp 1",]
      # Prepare plotting of PLSR-VIP results:
      coef<-plsrResult$coefficients[,,1]
      color_bars<-chillR::color_bar_maker(vipResult, coef, threshold=0.8, 
                                          col1="RED", col2="DARK GREEN", col3="DARK GREY")
      vipResult_select<-vipResult[order(vipResult,decreasing=TRUE)[1:50]]
      if(length(vipResult)<50) vipResult_select<-vipResult_select[1:length(vipResult)]
      col_select<-color_bars[order(vipResult,decreasing=TRUE)[1:50]]
      if(length(vipResult)<50) col_select<-col_select[1:length(vipResult)]
      # Plot PLSR-VIP results to file:
      png(file.path(outputPath, paste(i, "_PLS_VIP.png",sep="")),height=1400,width=1400)
      par(mar=c(5.1,55,4.1,2.1))
      barplot(rev(vipResult_select),horiz=TRUE,las=1,col=rev(col_select),cex.names=2,cex.axis=1,main="VIP for most important variables",cex.main=2,axes=FALSE)
      axis(side=1,cex.axis=2,lwd=5,padj=0.7)
      abline(v=0.8,lwd=3)
      dev.off()
      
      if (write_table){
        vipPlsResultTable<-cbind(vipResult,coef)
        colnames(vipPlsResultTable)<-c("VIP","Coefficient")
        #         if (write_table) write.csv(pls_tab,paste(result_path,resultnames[ress],"_pls_results.csv",sep=""))
        write.csv(vipPlsResultTable,file.path(outputPath, paste(i,"_pls_results.csv",sep="")))
      }
    }
    if (verbosity > 0)
      cat("VIP PLSR results written into directory: ", outputPath, "\n")
  }
  # Individual EVPI analysis:
  if (!is.null(individualEvpiNames)){
    if( !is.character(individualEvpiNames) )
      stop("\"individualEvpiNames\" must be a character vector.")
    else if( any(!individualEvpiNames %in% row.names(estimateObject)) ){
      if( individualEvpiNames!="all") 
        stop("\"individualEvpiNames\" must be a subset of variables defined in file 
              \",inputFilePath\" or =\"all\".")
      #if (individualEvpiNames=="all" )
      else
        individualEvpiNames=row.names(estimateObject)
    }
    if(verbosity > 0)
      cat("Performing Monte Carlo Simulation(s) for the IndividualEVPI Calculation(s):\n")
    individualEvpiResults<-individualEvpiSimulation(welfare=welfareFunction, 
                                                    currentEstimate=estimateObject, 
                                                    perfectProspectiveNames=individualEvpiNames,
                                                    numberOfModelRuns=numberOfModelRuns,
                                                    randomMethod=randomMethod,
                                                    functionSyntax=functionSyntax,
                                                    relativeTolerance=relativeTolerance,
                                                    verbosity=verbosity)
    
    if (verbosity > 1)
      print(sort(summary(individualEvpiResults)))
    individualEvpiSimulationSummaryFilePath<-file.path(outputPath,"individualEvpiSummary.csv")
    write.csv(sort(summary(individualEvpiResults), along=sortEvpiAlong)$summary$evi,
              individualEvpiSimulationSummaryFilePath)
    if (verbosity > 0)
      cat("IndividualEVPI simulation summary written into file: ", 
          individualEvpiSimulationSummaryFilePath, "\n")
  }
  
}
