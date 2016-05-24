#' Subset TAPS data and find weights
#'
#' Subset TAPS data by outcome and covariates, model the attrition rates, impute data for attrited individuals, and find weights for analysis
#'
#' @param interact A logical vector indicating if the function is to be run interactively. If TRUE arguments are not needed - default is TRUE
#' @param outcome A character vector of the names of outcome variables of interest. It is highly suggested that the outcome variables be entered starting with the earliest wave.
#' @param covars A character vector of the names of the covariate variables of interest
#' @param weight A logical argument specifying whether to use TAPS base weights or not - default is FALSE
#' @param refusedasNA A logical argument specifying whether to conisder the response 'Refused' as a missing value - default and suggested value is TRUE
#' @param method A character object indicating type of imputation to be used. hotdeck for hotdeck imputation, multi (default) for mulitple imputation, none for no imputation
#' @param na.delete A logical argument specifying whether to eliminate rows with NAs before calculating weights if method chosen is none - default is TRUE. Only set to FALSE if planning to use NA observations.
#' @param m A numeric argument specifying number of imputed data sets to produce if using multiple imputation - default is 5
#' @param pop.base A numeric object specifying which CPS data to use as a baseline. 1 is Dec. 2011, 2 is Dec. 2012, 3 is Dec. 2013. Default is 1.
#' @param trunc_at A numeric object specifying where to truncate the weights (what should the max weight be?) - default is 5
#' @param stringsAsFactors A logical vector indicating whether non-numeric variables should be factors rather than strings - default is TRUE
#'
#' @return An object of class weightTAPSoutput with the following slots:
#' \itemize{
#'  \item \code{df} List of dataframes (or single dataframe) of subsetted TAPS data with imputed data for missing values in covariates, with column for weights
#'  \item \code{attrit} A list of vectors of attrition rates for populations
#'  \item \code{stats} A list of population statistics from wave to wave
#'  }
#' @docType methods
#' @author David Carlson \email{carlson.david@@wustl.edu},  Michelle Torres \email{smtorres@@wustl},and Taeyong Park \email{t.park@@wustl.edu}
#' @examples
#' 
#' myOutcome <- c("APPRCONGS2","APPRCONGS6")
#' myCovars <- c("POLKNOW3S2","POLKNOW6S2") 
#' test<-weightTAPS(interact=FALSE,outcome=myOutcome,covars=myCovars,weight=FALSE,refusedasNA=TRUE,method="hotdeck",m=5,pop.base=1,trunc_at=5,stringsAsFactors=TRUE)
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{variablesTAPS}} \code{\link{subsetTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname weightTAPS
#' @aliases weightTAPS,ANY-method
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
#' If \code{method="multi"} is specified, \code{df} contains a list of \code{m} dataframes.
#' 
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
#' 
#' The information contained in both the \code{attrit} and \code{stats} slots can be graphically illustrated using the \code{plot(objectname)} function.
#' Two different types of plots are displayed after running the plot function: a dot chart and a set of trend plots.
#' The dot chart shows the differences between the sociodemographic composition of the sample in the first wave specified and the final subset dataframe.
#' This information is disaggregated by the following sociodemographic groups: Age and Gender, Ethnicity, Education, Income, Region and Metropolitan status, and Internet use.
#' The trend plots presented illustrate the changing composition of the sample by demographic group across the waves specified.
#' The lines shown in each plot correspond to the different categories within each of the groups mentioned.
#' The lines show the percentage of the final subset data belonging to each category by wave. 
#' The plots aim to show the variation in the composition of the sociodemographic groups through the waves specified.
#'  @export
setGeneric(name="weightTAPS",
           def=function(interact=TRUE,outcome=NULL,covars=NULL,weight=FALSE,refusedasNA=TRUE,method="multi",na.delete=TRUE,m=5,pop.base=1,trunc_at=5,stringsAsFactors=TRUE)
           {standardGeneric("weightTAPS")}
)

setMethod(f="weightTAPS",
          definition=function(interact=TRUE,outcome=NULL,covars=NULL,weight=FALSE,refusedasNA=TRUE,method="multi",na.delete=TRUE,m=5,pop.base=1,trunc_at=5,stringsAsFactors=TRUE){

            if(interact){
              outcome<-readline("Input the first outcome variable of interest:")
              i<-1
              while(outcome[i]!=""){
                outcome<-c(outcome,readline("Input the next outcome variable of interest \n (leave blank if done):"))
                i<-i+1
              }
              outcome<-outcome[-i]
              covars<-readline("Input the first covariate variable of interest \n (leave blank if none):")
              i<-1
              while(covars[i]!=""){
                covars<-c(covars,readline("Input the next covariate variable of interest \n (leave blank if done):"))
                i<-i+1
              }
              covars<-covars[-i]
              weight<-as.logical(readline("Enter TRUE if you would like to use the TAPS \n base weights, FALSE otherwise (FALSE is recommended):"))
              refusedasNA<-as.logical(readline("Enter TRUE if you would like to consider Refused \n responses as NA, FALSE otherwise (TRUE is recommended):"))
              if(is.na(covars[1])) covars<-NULL
              if(!is.null(covars)){
                method<-readline("Enter hotdeck if you would like to use \n hotdeck imputation, multi for multiple imputation, or none:")
                if(method=="multi") m<-as.numeric(readline("Enter number of datasets to be created using imputation:"))
              }else{
                method<-"none"
              }
              if(method == 'none') na.delete <- as.logical(readline('Enter TRUE if you would like to \n delete rows with NAs (TRUE recommended):'))
              pop.base<-as.numeric(readline("Enter 1 (suggested) to use CPS data \n from Dec. 2011, 2 for Dec. 2012, 3 for Dec. 2013:"))
              trunc_at<-as.numeric(readline('What value would you like to truncate the weights \n (5 is recommended):'))
              stringsAsFactors <- as.logical(readline('Enter TRUE if you would like strings \n as factors, FALSE otherwise:'))
              
            }
           
            data<-subsetTAPS(outcome,covars,weight,refusedasNA)
            if(method == 'none' & na.delete) data<-na.omit(data)
            if(!is.null(covars)){
              test<-apply(data[covars],2,function(x){all(is.na(x))})
              continue<-"continue"
              if(any(test)) continue<-readline(paste("Covariates selected ",covars[test]," have all missing values.\n To ignore, type continue, to terminate function type terminate: "))
              if(continue=="terminate") stop("Function terminated")
            }
            
            if(weight){
              data<-data[!is.na(data$baseweights),]
              data$baseweights <- data$baseweights*(nrow(data)/sum(data$baseweights)) 
            }
            names<-c(~agegend,~ethm,~educat,~regmetro,~incomcat,~ppnet)
            data[is.na(data$strata),"strata"]<-"NA"
            data$strata<-factor(data$strata)
            new.weights<-simpleWeight(data,pop.margins[pop.base][[1]],names,trunc_at)
            data<-cbind(data,new.weights)
            attrit.data<-attritTAPS(outcome)
            stats.data<-wavesTAPS(outcome)
            if(method=="none") df.data<-data
            if(method=="hotdeck") df.data<-hotdeckImp(data,outcome,covars)
            if(method=="multi") df.data<-multipleImp(data,outcome,covars,m)
            if(!stringsAsFactors){
              if(class(df.data)=='data.frame'){
                i <- sapply(df.data, is.factor)
                df.data[i] <- lapply(df.data[i], as.character)
              }
              else{
                for(i in 1:length(df.data)){
                  k <- sapply(df.data[[i]], is.factor)
                  df.data[[i]][k] <- lapply(df.data[[i]][k], as.character)
                }
              }
            }
            
            if(interact) print(new("weightTAPSoutput",df=df.data,attrit=attrit.data,stats=stats.data))
            return(new("weightTAPSoutput",df=df.data,attrit=attrit.data,stats=stats.data))
          }
)
