#' @title reproduceTableWithSourceDataByCiolkowski
#' @description Function reproduces Table, which shows the effect sizes reported by Ciolkowski identifying the type of design used in each study.
#' @author Lech Madeyski
#' @export reproduceTableWithSourceDataByCiolkowski
#' @examples
#' reproduceTableWithSourceDataByCiolkowski()
reproduceTableWithSourceDataByCiolkowski <- function(){
  reproducer::printXTable(data = reproducer::Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR, selectedColumns = 'c("Study", "Ref.", "Control", "Within-subjects", "Cross-over", "d_ByCiolkowski", "d_ByOriginalAuthors")', tableType = "latex", alignCells = "c(rep('l',6), rep('d{1.4}',1), rep('l',1))", digits = c(0,0,0,0,0,0,4,2), caption = "Source Data", label = "tab:SourceDataUsedByCiolkowski", fontSize = "footnotesize", captionPlacement = "top", alignHeader = c( "\\multicolumn{1}{l}{Study}", "\\multicolumn{1}{l}{Ref.}", "\\multicolumn{1}{l}{Control}", "\\multicolumn{1}{p{1.2cm}}{Within-subjects}", "\\multicolumn{1}{p{1.0cm}}{Cross-over}", "\\multicolumn{1}{p{1.9cm}}{$d$ (by Ciol-kowski)}", "\\multicolumn{1}{p{2.1cm}}{$d$ (by origi-nal authors)}" ) )
}

#' @title reproduceTableWithEffectSizesBasedOnMeanDifferences()
#' @description Function reproduces Table, which shows the effect sizes based on mean differences.
#' @author Lech Madeyski
#' @export reproduceTableWithEffectSizesBasedOnMeanDifferences
#' @examples
#' reproduceTableWithEffectSizesBasedOnMeanDifferences()
reproduceTableWithEffectSizesBasedOnMeanDifferences <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  rownames(dataMK)<-dataMK$Study
  printXTable(data=dataMK, selectedColumns='c("Study", "Ref.", "ExpDesign", "Teams", "n_PBR", "n_C", "M_PBR", "M_C", "Diff", "Inc", "V_C", "V_D")', tableType="latex", alignCells="c(rep('l',4), rep('d{2.0}',3), rep('d{1.3}',2), rep('d{2.3}',1), rep('d{3.0}',1), rep('d{1.3}',2))", digits=c(0,0,0,0,0,0,0,3,3,3,0,3,3), caption="Effect Sizes based on Mean Differences", label="tab:MDEffectSize", fontSize="footnotesize", captionPlacement="top", alignHeader=c("\\multicolumn{1}{l}{Study}", "\\multicolumn{1}{l}{Ref.}", "\\multicolumn{1}{l}{ExpDesign}", "\\multicolumn{1}{l}{$Teams$}", "\\multicolumn{1}{l}{$n_{PBR}$}", "\\multicolumn{1}{l}{$n_C$}", "\\multicolumn{1}{l}{$M_{PBR}$}", "\\multicolumn{1}{l}{$M_C$}", "\\multicolumn{1}{l}{$\\overline{\\mathit{diff}}$}", "\\multicolumn{1}{l}{$\\%Inc$}", "\\multicolumn{1}{l}{$V_C$}", "\\multicolumn{1}{l}{$V_D$}"))
}

#' @title reproduceForestPlotRandomEffects()
#' @description Function reproduces Forest Plot of a Random-Effects Meta-analysis of Mean Differences.
#' @author Lech Madeyski
#' @export reproduceForestPlotRandomEffects
#' @examples
#' reproduceForestPlotRandomEffects()
reproduceForestPlotRandomEffects <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  res<-metafor::rma(yi=dataMK$Diff,vi=dataMK$V_D,data=dataMK,method="REML")
  metafor::forest(res, xlab="Overall Effect Size", slab=reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR$Study)
  graphics::text(-1,19," Study")
  graphics::text(0.9,19, "Effect Size   [95% CL]")
}

#' @title reproduceMixedEffectsAnalysisWithExperimentalDesignModerator()
#' @description Function reproduces Mixed-Effects Analysis with Experimental Design as a Moderator.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsAnalysisWithExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsAnalysisWithExperimentalDesignModerator()
reproduceMixedEffectsAnalysisWithExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_D, data=dataMK, method="REML", mods=~factor(ExpDesign)-1)
  resWithExpDesignMod
}

#' @title reproduceMixedEffectsForestPlotWithExperimentalDesignModerator()
#' @description Function reproduces Forest Plot of a Mixed Effects Meta-analysis of Mean Differences with Experimental Design as a Moderator Variable.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsForestPlotWithExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsForestPlotWithExperimentalDesignModerator()
reproduceMixedEffectsForestPlotWithExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_D, data=dataMK, method="REML", mods=~factor(ExpDesign)-1)
  metafor::forest(resWithExpDesignMod,xlab="Overall Effect Size", slab=reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR$Study)
  graphics::text(-1,19," Study")
  graphics::text(0.9,19, "Effect Size   [95% CL]")
}


#' @title reproduceTableWithPossibleModeratingFactors()
#' @description Function reproduces Table with possible moderating factors.
#' @author Lech Madeyski
#' @export reproduceTableWithPossibleModeratingFactors
#' @examples
#' reproduceTableWithPossibleModeratingFactors()
reproduceTableWithPossibleModeratingFactors <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  rownames(dataMK) <- dataMK$Study
  printXTable(data=dataMK, selectedColumns='c("Study", "Ref.", "ControlType", "ParticipantsType", "TeamType", "ArtefactType", "AssociatedWithBasili")', tableType="latex", alignCells="c(rep('c',8))", digits=0, caption="Possible Moderating Factors", label="tab:Moderators", fontSize="footnotesize", captionPlacement="top", alignHeader=c("\\multicolumn{1}{c}{Study}", "\\multicolumn{1}{c}{Ref.}", "\\multicolumn{1}{c}{Control}","\\multicolumn{1}{c}{Participants}", "\\multicolumn{1}{c}{TeamType}", "\\multicolumn{1}{c}{ArtefactType}", "\\multicolumn{1}{c}{AssociatedWithBasili}"))
}


#' @title reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator()
#' @description Function reproduces Mixed-Effects Analysis using Subject Specific Estimated Variance with Experimental Design as a Moderator.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator()
reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignUsingSubjectSpecificVarianceMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_Alt, data=dataMK,method="REML", mods=~factor(ExpDesign)-1)
  resWithExpDesignUsingSubjectSpecificVarianceMod
}


