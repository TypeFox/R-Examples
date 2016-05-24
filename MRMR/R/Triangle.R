#' is.Triangle
#' @description
#' Tests whether the object is a triangle
#' @return 
#' TRUE if the object is a triangle, FALSE if it is not
#' @export
#' @param object The object to be tested
is.Triangle = function(object)
{
  is(object, "Triangle")
}

checkTriangle = function(object)
{
  errors = character()
  if (length(errors) == 0) TRUE else errors
}

#' Triangle class
#' 
#' @description
#' Triangle is an S4 class used to store aggregated loss data. All triangles must have a defined set of OriginPeriods, a defined set of DevelopmentIntervals 
#' and data along those axes. A triangle may carry additional descriptive information such as line of business, geographic region and so on.
#' 
#' @details
#' One will rarely, if ever use the setClass method directly. The function \code{\link{newTriangle}} will generally be used to create a new Triangle object
#' 
#' @seealso \code{\link{newTriangle}}
#' 
#' @name Triangle-class
#' @rdname Triangle-class
#' @exportClass Triangle
#' 
setClass("Triangle"
         , representation(TriangleData = "data.frame"
                          , TriangleName = "character"
                          , OriginPeriodType = "character"
                          , DevelopmentInterval = "Period"
                          , StaticMeasures = "character"
                          , StochasticMeasures = "character"
                          , Groups = "character")
#          , sealed = TRUE
#         , validity = #some function
)

#' Create a Triangle object.
#' @param TriangleData A dataframe 
#' @param OriginPeriods The name of the column in the TriangleData which holds the origin period.
#' @param DevelopmentLags The column which holds the development lags.
#' @param OriginEnd If the OriginPeriods argument refers to the start date of an origin period, this column holds the end dates.
#' @param OriginLength If origin period is not an interval, this is used to construct the origin period.
#' @param StartDay If origin period is not an interval, this is used to construct the origin period.
#' @param StartMonth If origin period is not an interval, this is used to construct the origin period.
#' @param DevelopmentPeriod If DevelopmentLags is not a period object, this is used to contruct DevelopmentLags.
#' @param EvaluationDates A vector of dates corresponding to the data in TriangleData.
#' @param OriginPeriodType A character value describing the type of origin period.
#' @param TriangleName A character value used to refer to the Triangle object.
#' @param StaticMeasures A character vector which names the static measures in the Triangle object.
#' @param StochasticMeasures A character vector which names the stochastic measures in the Triangle object.
#' @param Groups A character vector which names the column which contains grouping data.
#' @param Cumulative Boolean indicating if the stochastic measures are cumulative or incremental.
#' @param Verbose Boolean indicating whether or not warnings should be displayed.
#' 
#' @export newTriangle
#' @include TriangleAdjustMeasures.R
#' @include CreateDevelopmentLags.R
#' @include CreateEvaluationDates.R
#' @include CreateOriginPeriods.R
#' 
# User-friendly constructor
newTriangle = function(TriangleData
                       , OriginPeriods = NULL
                       , DevelopmentLags = NULL
                       , OriginEnd = NULL
                       , OriginLength = years(1)
                       , StartDay = 1
                       , StartMonth = 1
                       , DevelopmentPeriod = months(1)
                       , EvaluationDates = NULL
                       , OriginPeriodType = "Accident Year"
                       , TriangleName = NULL
                       , StaticMeasures = NULL
                       , StochasticMeasures = NULL
                       , Groups = NULL
                       , Cumulative = TRUE
                       , Verbose = TRUE)
{
  arguments <- as.list(match.call())
  
  OriginPeriods = eval(arguments$OriginPeriods, TriangleData)
  if(is.null(OriginPeriods)) stop ("No origin period was specified.")
  
  if (!is.interval(OriginPeriods)) {
    OriginPeriods = CreateOriginPeriods(OriginPeriods, OriginEnd, OriginLength, StartDay, StartMonth, Verbose)
  }
  
  DevelopmentLags = eval(arguments$DevelopmentLags, TriangleData)
  if(is.null(DevelopmentLags)) stop ("No development lag information was provided.")
  
  if (!is.period(DevelopmentLags)){
    DevelopmentLags = CreateDevelopmentLags(DevelopmentLags, DevelopmentPeriod, EvaluationDates, OriginPeriods, Verbose)
  }
  
  CommonDevInterval = DevelopmentLags[order(DevelopmentLags)]
  CommonDevInterval = CommonDevInterval[1]
  if(Verbose){
    DevInteger = DevelopmentLags / CommonDevInterval  
  } else {
    DevInteger = suppressMessages(DevelopmentLags / CommonDevInterval)
  }
  
  if(is.null(EvaluationDates)) {
    EvaluationDates = CreateEvaluationDates(OriginPeriods, DevelopmentLags)
  } else {
    # Throw a warning if the evaluation dates which were passed in are not consistent with what they ought to be.
  }
  
  # It's possible that the user has fed data with overlap. 
  # Annual data with two start dates in the same year and annual development period, zB.
  # I might decide to check for this and throw a warning, but for now, I'll just blame the user.
  
  dfNewTriangleData = data.frame(OriginPeriod = OriginPeriods
                                 , DevelopmentLag = DevelopmentLags
                                 , EvaluationDate = EvaluationDates
                                 , DevInteger = DevInteger)
  
  dfNewTriangleData$OriginPeriodStart = int_start(dfNewTriangleData$OriginPeriod)
  dfNewTriangleData$OriginPeriodEnd = int_end(dfNewTriangleData$OriginPeriod)
  
  dfNewTriangleData$CalendarPeriodStart = dfNewTriangleData$EvaluationDate - CommonDevInterval + days(1)
  dfNewTriangleData$CalendarPeriodEnd = dfNewTriangleData$EvaluationDate 
  dfNewTriangleData$CalendarPeriod = with(dfNewTriangleData, new_interval(CalendarPeriodStart, CalendarPeriodEnd))
  
  if (is.null(Groups))
  {
    dfNewTriangleData$Group = "All"
    Groups = "Group"
  } else {
    dfNewTriangleData = cbind(dfNewTriangleData, TriangleData[Groups])
  }
  
  if (!is.null(StaticMeasures)) {
    dfNewTriangleData = cbind(dfNewTriangleData, TriangleData[StaticMeasures])
  } else {
    StaticMeasures = ""
  }
  
  if (is.null(StochasticMeasures)) stop ("You've not supplied any stochastic measures for this triangle. Idiot.")

  dfStochasticMeasures = TriangleData[StochasticMeasures]
  stochasticMeasureNames = FormMeasureNames(dfStochasticMeasures, Cumulative)
  names(dfStochasticMeasures) = stochasticMeasureNames

  dfNewTriangleData = cbind(dfNewTriangleData, dfStochasticMeasures) 
  
  if(Cumulative) {
    dfNewTriangleData = CreateIncrementals(dfNewTriangleData, stochasticMeasureNames, Groups)
    stochasticMeasureNames = c(stochasticMeasureNames, gsub("Cumulative", "Incremental", stochasticMeasureNames))
  } else {
    dfNewTriangleData = CreateCumulative(dfNewTriangleData, stochasticMeasureNames, Groups)
    stochasticMeasureNames = c(stochasticMeasureNames, gsub("Incremental","Cumulative", stochasticMeasureNames))
  }
  
  dfNewTriangleData = CreatePriors(dfNewTriangleData, stochasticMeasureNames, Groups)
  
  if (is.null(TriangleName)) TriangleName = ""
  
  row.names(dfNewTriangleData) = NULL
  
  tri = new("Triangle"
            , TriangleData = dfNewTriangleData
            , TriangleName = TriangleName
            , OriginPeriodType = OriginPeriodType
            , DevelopmentInterval = CommonDevInterval
            , StaticMeasures = StaticMeasures
            , StochasticMeasures = CleanMeasureNames(stochasticMeasureNames)
            , Groups = Groups)
  
  tri
}

TestDataFrame = function(){
  dfTest = data.frame(AccidentYear = c(2002, 2002, 2002, 2003, 2003, 2004)
             , Month = c(12, 24, 36, 12, 24, 12)
             , Paid = c(2318,  7932, 13822, 1743,  6240, 2221)
             , Reported = c(12811, 20370, 26656, 9651, 16995, 16995)
             , EP = c( 61183,  61183,  61183,  69175,  69175,  99322))
  
  dfTest
}