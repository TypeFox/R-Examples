#' A simple multiple-group ACE model with the \pkg{lavaan} package.
#' 
#' @export
#' @description This function uses the \pkg{lavaan} package to estimate a univariate ACE model, using multiple groups.  
#' Each group has a unique value of \code{R} (i.e., the \emph{R}elatedness coefficient).
#' @usage AceLavaanGroup(dsClean, estimateA=TRUE, estimateC=TRUE, printOutput=FALSE)
#' 
#' @param dsClean The \code{data.frame} containing complete cases for the \code{R} groups to be included in the estimation.
#' @param estimateA Should the \emph{A} variance component be estimated?  A^2 represents the proportion of variability due to a shared genetic influence.
#' @param estimateC Should the \emph{C} variance component be estimated?  C^2 represents the proportion of variability due to a shared environmental influence.
#' @param printOutput Indicates if the estimated parameters and fit statistics are printed to the console.
#' 
#' @details The variance component for \emph{E} is always estimated, while the \emph{A} and \emph{C} estimates can be fixed to zero (when \code{estimateA} and/or \emph{estimateC} are set to \code{FALSE}).
#' @return An \code{AceEstimate} object.
#' @references The \pkg{lavaan} package is developed by Yves Rosseel at Ghent University.  
#' Three good starting points are the package home page (\url{http://lavaan.ugent.be/}), the documentation (\url{http://cran.r-project.org/package=lavaan}) 
#' and the JSS paper.
#' 
#' Rosseel, Yves (2012), \href{http://www.jstatsoft.org/v48/i02/}{lavaan: An R Package for Structural Equation Modeling}. \emph{Journal of Statistical Software, 48}, (2), 1-36.
#' @author Will Beasley
#' @note Currently, the variables in \code{dsClean} must be named \code{O1}, \code{O2} and \code{R}; the letter `O' stands for \emph{O}utcome.  This may not be as restrictive as it initially seems, because \code{dsClean} is intented to be produced by \code{CleanSemAceDataset}.  If this is too restrictive for your uses, we'd like to here about it (\emph{please email wibeasley at hotmail period com}).
#' @seealso \code{\link{CleanSemAceDataset}}.  Further ACE model details are discussed in our package's \href{http://cran.r-project.org/package=NlsyLinks}{vignettes}.
#' @keywords ACE
#' @examples
#' library(NlsyLinks) #Load the package into the current R session.
#' dsLinks <- Links79PairExpanded #Start with the built-in data.frame in NlsyLinks
#' dsLinks <- dsLinks[dsLinks$RelationshipPath=='Gen2Siblings', ] #Use only Gen2 Siblings (NLSY79-C)
#' 
#' oName_S1 <- "MathStandardized_S1" #Stands for Outcome1
#' oName_S2 <- "MathStandardized_S2" #Stands for Outcome2
#' 
#' dsGroupSummary <- RGroupSummary(dsLinks, oName_S1, oName_S2)
#' dsClean <- CleanSemAceDataset(dsDirty=dsLinks, dsGroupSummary, oName_S1, oName_S2)
#' 
#' ace <- AceLavaanGroup(dsClean)
#' ace
#' 
#' #Should produce:
#' # [1] "Results of ACE estimation: [show]"
#' #     ASquared     CSquared     ESquared    CaseCount 
#' #    0.6681874    0.1181227    0.2136900 8390.0000000 
#' 
#' library(lavaan) #Load the package to access methods of the lavaan class.
#' GetDetails(ace)
#' 
#' #Exmaine fit stats like Chi-Squared, RMSEA, CFI, etc.
#' fitMeasures(GetDetails(ace)) #The function 'fitMeasures' is defined in the lavaan package.
#' 
#' #Examine low-level details like each group's individual parameter estimates and standard errors.
#' summary(GetDetails(ace))
#' 
#' #Extract low-level details. This may be useful when programming simulations.
#' inspect(GetDetails(ace), what="converged") #The lavaan package defines 'inspect'.
#' inspect(GetDetails(ace), what="coef")



AceLavaanGroup <-
function( dsClean, estimateA=TRUE, estimateC=TRUE, printOutput=FALSE) {
  # library(lavaan)
  # library(stringr)
 
  rLevels <- base::sort(base::unique(dsClean$R))
  #These five lines enumerate the path coefficient labels to be inserted into the model statement.
  #rString <- stringr::str_c(rLevels, collapse=", ") #The output is typically "1, 0.5, 0.375, 0.25"
  rString <- base::paste(rLevels, collapse=", ") #The output is typically "1, 0.5, 0.375, 0.25"
  # aString <- str_c(rep("a", length(rLevels)), collapse=",") #The output is typically "a,a,a,a"
  # cString <- str_c(rep("c", length(rLevels)), collapse=",") #The output is typically "c,c,c,c"
  # eString <- str_c(rep("e", length(rLevels)), collapse=",") #The output is typically "e,e,e,e"
  # intString <- str_c(rep("int", length(rLevels)), collapse=",") #The output is typically "int,int,int,int"
  groupCount <- base::length(rLevels)
  
  #Generate the model statements that will exist in all models (ie, ACE, AE, AC, & E)
  modelBase <- base::paste("
    O1 ~~ 0 * O2                          #The manifest variables are uncorrelated.
    O1 + O2 ~ rep('int',", groupCount, ") * 1           #The manifest variables are fed the same intercept (for all groups).
    
    E1 =~ rep('e',", groupCount, ") * O1  #Declare the contributions of E to Subject1 (for all groups).
    E2 =~ rep('e',", groupCount, ") * O2  #Declare the contributions of E to Subject2 (for all groups).
    
    E1 ~~ 0 * E2                          #The Es are uncorrelated
    E1 ~~ 1 * E1                          #The Es have a variance of 1
    E2 ~~ 1 * E2
    e2 := e * e                           #Declare e^2 for easy point and variance estimation.
  ")
  

  #Generate the model statements that will exist in all models estimating the A component.
  modelA <- base::paste("
    A1 =~ rep('a',", groupCount,") * O1   #Declare the contributions of A to Subject1 (for all groups).
    A2 =~ rep('a',", groupCount,") * O2   #Declare the contributions of A to Subject2 (for all groups).
    A1 ~~ c(", rString, ") * A2           #Declare the genetic relatedness between Subject1 and Subject2. This coefficient differs for all groups.
    
    A1 ~~ 1 * A1                          #The As have a variance of 1
    A2 ~~ 1 * A2
    a2 := a * a                           #Declare a^2 for easy point and variance estimation.
  ")

  #Generate the model statements that will exist in all models estimating the C component.
  modelC <- base::paste("
    C1 =~ rep('c',", groupCount,") * O1   #Declare the contributions of C to Subject1 (for all groups).
    C2 =~ rep('c',", groupCount,") * O2   #Declare the contributions of C to Subject2 (for all groups).
    
    C1 ~~ 1 * C2                          #The Cs are perfectly correlated. !!Note this restricts the sample to immediate families!!
    C1 ~~ 1 * C1                          #The Cs have a variance of 1
    C2 ~~ 1 * C2
    c2 := c * c                           #Declare c^2 for easy point and variance estimation.
  ")
  
  #If the A/C component's excluded: (1) overwrite the code and (2) declare a2/c2, so the parameter value can be retrieved later without if statements.
  if( !estimateA ) modelA <- "\n a2 := 0 \n"
  if( !estimateC ) modelC <- "\n c2 := 0 \n"
  
  model <- base::paste(modelBase, modelA, modelC) #Assemble the three parts of the model.
  
  #Run the model and review the results
  fit <- lavaan::lavaan(model, data=dsClean, group="GroupID", missing="listwise", information="observed")
  if( printOutput ) lavaan::summary(fit)
  #lavaanify(model) #lavaanify(modelA)
  #parseModelString(modelA)
  #lavaanNames(modelA)
  #parTable(fit)
  #parameterEstimates(fit)
  #inspect(fit)
  #names(fitMeasures(fit)) #str(fitMeasures(fit)[c("chisq", "df")])
  if( printOutput ) print(paste("Chi Square: ", lavaan::fitMeasures(fit)[["chisq"]])) #Print the Chi Square value
  
  #Extract the UNSCALED ACE components.
  est <- lavaan::parameterEstimates(fit)
  a2 <- est[est$label=="a2", "est"]
  c2 <- est[est$label=="c2", "est"]
  e2 <- est[est$label=="e2", "est"]
  
  #variogram-like diagnostics
  # plot(dsGroupSummary$R, dsGroupSummary$Covariance, pch=(4-3*dsGroupSummary$Included))
  # plot(dsGroupSummary$R, dsGroupSummary$Correlation, pch=(4-3*dsGroupSummary$Included))
  # text(dsGroupSummary$PairCount, x=dsGroupSummary$R, y=dsGroupSummary$Correlation)
  # 
  # ggplot(data=dsGroupSummary, aes(x=R, y=Correlation, label=PairCount, color=Included)) + #substitute(N == b, list(b=dsGroupSummary$PairCount)) )) + #geom_line()
  #   layer(geom="line") + layer(geom="text") + scale_x_reverse(limits=c(1,0)) + scale_y_continuous(limits=c(0,1))
  # ggplot(data=dsGroupSummary, aes(x=R, y=Covariance, label=PairCount, color=Included)) + #substitute(N == b, list(b=dsGroupSummary$PairCount)) )) + #geom_line()
  #   layer(geom="line") + layer(geom="text") + scale_x_reverse(limits=c(1,0))
  
  components <- base::as.numeric(base::cbind(a2, c2, e2)[1,] / (a2 + c2 + e2)) #The 'as.numeric' gets rid of the vector labels.
  if( printOutput ) base::print(components) #Print the unity-SCALED ace components.
  caseCount <- base::nrow(dsClean)
  details <- base::list(lavaan=fit)
  #print(paste("R Levels excluded:",  stringr::str_c(rLevelsToExclude, collapse=", "), "; R Levels retained:", rString)) #Print the dropped & retained groups.
  ace <- NlsyLinks::CreateAceEstimate(
    aSquared  = components[1], 
    cSquared  = components[2], 
    eSquared  = components[3], 
    caseCount = caseCount, 
    details   = details
  )
  return( ace )
}
