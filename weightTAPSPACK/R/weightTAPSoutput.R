#' An object outputted by weightTAPS function
#' 
#' Objects of class \code{weightTAPSoutput} are created by the \code{weightTAPS} function
#'
#' An object of the class `weightTAPSoutput' has the following slots:
#' \itemize{
#' \item \code{df} List of dataframes of subsetted TAPS data with imputed data for non-response in covariates, with column for weights - if method none chosen, of class data.frame
#' \item \code{attrit} A list of vectors of attrition rates by sociodemographic groups
#' \item \code{stats} A list of population statistics from wave to wave
#' }
#'
#' @author David Carlson \email{carlson.david@@wustl.edu},  Michelle Torres \email{smtorres@@wustl}, and Taeyong Park \email{t.park@@wustl.edu}
#' @docType methods
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{subsetTAPS}} \code{\link{wavesTAPS}} \code{\link{simpleWeight}} \code{\link{attritTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}  \code{\link{TAPScum}} \code{\link{TAPSimputeddemographics}}
#' @rdname weightTAPSoutput
#' @aliases weightTAPSoutput-class initialize,weightTAPS-method getdf,weightTAPSoutput-method getattrit,weightTAPSoutput-method getstats,weightTAPSoutput-method print,weightTAPSoutput-method show,weightTAPSoutput-method plot,weightTAPSoutput-method weightTAPSoutput
#' @details
#' The slot \code{df} contains the dataframe(s) that represent the final subset data.
#' The final subset data keeps only the outcome and covariates of interest specified. It also accounts for non-response in the outcome variable and attrition. 
#' Respondents with missing values in the outcome variable for any of the waves desired are removed from the final dataframe.
#' Missing values in sociodemographic variables are imputed for the respondents that remain in the sample, in order to compute their proper weight.
#' The missing covariate data can be left as is, by specifying \code{method="none"}.
#' If imputation is desired, the argument \code{method} can be set to 'multi' for multiple imputation, or 'hotdeck' for hotdeck imputation.
#' If multiple imputation is done, the argument \code{m} should be set to the number of imputed dataframes to be created.
#'
#'
#' The slot \code{attrit} is a list of attrition rates from the first wave specified in the outcome argument.
#' Each quantity represents the percentage of people (by demographic group) that attritted TAPS through the waves specified.
#' It compares the initial composition of each demographic group (from the oldest wave specified) to the composition of the same demographic group in the final subset data delivered by \code{weightTAPS()}.
#' Large values, particularly large values relative to other values in the same sociodemographic category, indicate high rates of attrition.
#'
#'    
#' The slot \code{stats} lists each sociodemographic group's share of the overall population as represented in the final sample for each outcome.
#' The information contained in both the \code{attrit} and \code{stats} slots can be graphically illustrated using the \code{plot(objectname)} function.
#' Two different types of plots are displayed after running the plot function: a dot chart and a set of trend plots.
#' The dot chart shows the differences between the sociodemographic composition of the sample in the first wave specified and the final subset dataframe.
#' This information is disaggregated by the following sociodemographic groups: Age and Gender, Ethnicity, Education, Income, Region and Metropolitan status, and Internet use.
#' The trend plots presented illustrate the changing composition of the sample by demographic group across the waves specified.
#' The lines shown in each plot correspond to the different categories within each of the groups mentioned.
#' The lines show the percentage of the final subset data belonging to each category by wave. 
#' The plots aim to show the variation in the composition of the sociodemographic groups through the waves specified.
#'
#'
#' \code{print(objectname)} will show a summary of the attrition rates.
#' \code{show(objectname)} will run the \code{print} function.
#' @export
setClass(Class="weightTAPSoutput", 
         representation = representation(
           df = "ANY",
           attrit = "list",
           stats = "list"
         ),
         prototype = prototype(
           df=list(),
           attrit=list(),
           stats=list()
         )
)

setMethod("initialize", "weightTAPSoutput", 
          function(.Object, ...){
            value=callNextMethod()
            return(value)
          }
) 

#' @rdname weightTAPSoutput
#' @export
setGeneric("getdf",
           function(object="weightTAPSoutput")  {
             standardGeneric("getdf")
           }
)

setMethod("getdf", "weightTAPSoutput",
          function(object){ 
            return(object@df)
          }
)

#' @rdname weightTAPSoutput
#' @export
setGeneric("getattrit",
           function(object="weightTAPSoutput")  {
             standardGeneric("getattrit")
           }
)

setMethod("getattrit", "weightTAPSoutput",
          function(object){ 
            return(object@attrit)
          }
)

#' @rdname weightTAPSoutput
#' @export
setGeneric("getstats",
           function(object="weightTAPSoutput")  {
             standardGeneric("getstats")
           }
)

setMethod("getstats", "weightTAPSoutput",
          function(object){ 
            return(object@stats)
          }
)

#' @export
setMethod(f="print", "weightTAPSoutput",
          definition=function(x){ 
            cat("The dataframe can be accessed with the function getdf or with name@df.\n")
            cat("The population statistics from each wave can be accessed with the function getstats or with name@stats.\n")
            cat("Use show or plot function to see changing population statistics.\n")
            cat("\nThe vectors of attrition rates\n")
            print(x@attrit)
          }
)

#' @export
setMethod(f="show", "weightTAPSoutput",
          definition=function(object){ 
            print(object)
          }
)

#' @export
setMethod(f="plot", "weightTAPSoutput",
          definition=function(x){
            statWaves<-x@stats
            holdpar <- par(no.readonly=TRUE)
            par(ask=TRUE, oma=c(4,1,1,1))
            dotchart(as.matrix(unlist(x@attrit))[,1],labels=names((unlist(x@attrit))),cex=0.5, 
                     main="Atrrition Rates", col=rainbow(31, start=0.7, end=0.1))            
            for (i in 1:length(statWaves)){
              par(ask=TRUE, oma=c(4,1,1,1), mar=c(5.1,4.1,4.1,2.1))
              Waves <- colnames(statWaves[[i]])
              nWaves<-ncol(statWaves[[i]])
              Groups<-rownames(statWaves[[i]])
              nGroups<-nrow(statWaves[[i]])
              namesSocio<-c("Age/Gender", "Ethnicity", "Education", "Region", "Income", "Internet users")
              xrange <- seq(1, nWaves,1)
              ymax <- round(max(statWaves[[i]]),0)+1
              vec.x<-rep(xrange, each=nGroups)
              vec.y<-as.numeric(statWaves[[i]])
              plot(vec.x, vec.y, type="n", xaxt="n",
                   ylab="% of sample", main=namesSocio[[i]], xlab="Outcome", ylim=c(0, ymax), xlim=c(1,nWaves))
              axis(1,at=xrange,labels=Waves)
              colors <- rainbow(nGroups) 
              linetype <- c(1:nGroups) 
              plotchar <- seq(18,18+nGroups,1)
              txt.x<-NULL
              for (j in 1:length(vec.x)){
                if(j %% 2==0){txt.x[j]<-vec.x[j]+0.02}else{txt.x[j]<-vec.x[j]-0.02}
              }
              for (k in 1:nGroups) { 
                lines(statWaves[[i]][k,], type="o", col=colors[k], lty=linetype[k], pch=plotchar[k])
              } 
              text(txt.x, vec.y, labels=vec.y, cex=0.6)
              par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
              plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              legend("bottom", legend=Groups, xpd=TRUE, horiz=TRUE, cex=0.6, col=colors,
                     pch=plotchar, lty=linetype, title="Groups")
            }
            par(holdpar)
          }
)
