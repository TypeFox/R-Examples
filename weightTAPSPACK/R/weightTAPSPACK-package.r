#' weightTAPSPACK
#'
#' The weightTAPSPACK subsets The American Panel Survey (TAPS) data by outcome and covariates, models the attrition rates, imputes data for attrited individuals, and finds weights for analysis. 
#' @name weightTAPSPACK
#' @docType package
#' @author David Carlson \email{carlson.david@@wustl.edu}, Michelle Torres \email{smtorres@@wustl}, and Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{subsetTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @details
#' This package is meant to subset The American Panel Survey (TAPS) data by outcome and by covariate variables of interest through the function \code{\link{weightTAPS}}.
#' The subsetting process accounts for respondents attriting from at least one of the waves under analysis, as well as for outcome non-response.
#' The variables of interest must be entered exactly as named in the TAPS dataframe. See \url{http://taps.wustl.edu/data-archive} or use the \code{\link{variablesTAPS}} function
#' to explore the names of the variables by wave. It is important to revise the particular features of each of the variables of interest.
#' 
#' 
#' It is strongly suggested that the outcome variables be entered starting with the earliest wave for easier interpretation of the attrition rates.
#' Other arguments are listed in the help file of the \code{weightTAPS()} function, and must be considered based on the user's needs.
#' The function can be run in interactive mode by simply running \code{weightTAPS()}. The user must answer the questions based on her needs. 
#' 
#' 
#' The function \code{weightTAPS()} should be assigned to an object in order to conduct the analysis of TAPS.
#' \code{weightTAPS()} returns a subset of the complete TAPS dataset that includes only the outcome variable and covariates specified by the user, a set of standard demographics
#' and a new variable with the corresponding weight for each respondent.
#' 
#' 
#' It also retains the respondents that gave an answer to the outcome variable of interest through the waves specified by the user.
#' Respondents that attritted or did not provide an answer to the outcome variable for any of the waves under analysis are removed from the subset data.
#' Missing values in sociodemographic variables are imputed for the respondents that remain in the sample, in order to compute their proper weight.
#' The missing values observed in the covariates of the remaining respondents are imputed through the method selected by the user.
#' Once the TAPS data is subset, the function calculates weights based on the demographic group membership of the respondents in the final subset.
#' These weights will be appended to the end of the data frame(s) with the column name \code{new.weights}.
#' 
#' 
#' The output (see \code{\link{weightTAPSoutput}}) is of class weightTAPSoutput. This class implies the existence of certain slots that save useful information.
#' These slots are \code{df}, \code{attrit} and \code{stats}.
#' 
#' 
#' The slot \code{df} contains the dataframe(s) that represent the final subset data.
#' The final subset data keeps only the outcome and covariates of interest specified as well as a set of demographics and the new dynamic weights.
#' It also accounts for non-response in the outcome variable and attrition across waves through the waves specified.  
#' Respondents with missing values in the outcome variable for any of the waves desired are removed from the final dataframe.
#' 
#' 
#' The missing covariate data can be left as is, by specifying \code{method="none"}.
#' If imputation is desired, the argument \code{method} can be set to 'multi' for multiple imputation, or 'hotdeck' for hotdeck imputation.
#' If multiple imputation is done, the argument \code{m} should be set to the number of imputed dataframes to be created.
#' Depending on the imputation method selected, \code{df} can be a list of \code{m} elements or a list containing a single element. Each element of \code{df} stores a dataframe.
#' If \code{method="multi"} is specified, \code{df} contains a list of \code{m} dataframes. \n
#' 
#' To access the dataframes, use \code{getdf(objectname)}. Objectname corresponds to the object where the value of \code{weightTAPS()} was originally stored.
#' If hotdeck or no imputation was used, the final dataset is the first element of the \code{df} list, and can be accessed with \code{getdf(objectname)[[1]]}.
#' 
#' 
#' The slot \code{attrit} is a list of attrition rates from the first wave specified in the outcome argument.
#' Each quantity represents the percentage of people (by demographic group) that attritted TAPS through the waves specified.
#' It compares the initial composition of each demographic group (from the oldest wave specified) to the composition of the same demographic group in the final subset data delivered by \code{weightTAPS()}.
#' Large values, particularly large values relative to other values in the same sociodemographic category, indicate high rates of attrition.
#' It is important to highlight that high rates of attrition may cause problems in data analysis.
#' The slot \code{stats} lists each sociodemographic group's share of the overall population as represented in the final sample for each outcome.
#' 
#' The information contained in both the \code{attrit} and \code{stats} slots can be graphically illustrated using the \code{plot(objectname)} function.
#' Two different types of plots are displayed after running the plot function: a dot chart and a set of trend plots.
#' The dot chart shows the differences between the sociodemographic composition of the sample in the first wave specified and the final subset dataframe.
#' This information is disaggregated by the following sociodemographic groups: Age and Gender, Ethnicity, Education, Income, Region and Metropolitan status, and Internet use.
#' The trend plots presented illustrate the changing composition of the sample by demographic group across the waves specified.
#' The lines shown in each plot correspond to the different categories within each of the groups mentioned.
#' The lines show the percentage of the final subset data belonging to each category by wave. 
#' The plots aim to show the variation in the composition of the sociodemographic groups through the waves specified.
NULL
