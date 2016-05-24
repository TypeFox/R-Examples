is.TriangleModel = function(object)
{
  is(object, "TriangleModel")
}

checkTriangleModel = function(object)
{
  errors = character()

  if (length(errors) == 0) TRUE else errors
}

setOldClass("htest")

#' TriangleModel class
#' 
#' @description
#' Triangle is an S4 class used to store a model fit to a Triangle object.
#' 
#' @details
#' Some stuff
#' 
#' @seealso \code{\link{Triangle-class}}
#' 
#' @name TriangleModel-class
#' @rdname TriangleModel-class
#' @exportClass TriangleModel
#' 
setClass("TriangleModel"
         , representation(ModelData = "data.frame"
                          , Response = "character"
                          , Predictor = "character"
                          , FitCategory = "character"
                          , Alpha = "numeric"
                          , Tail = "numeric"
                          , Fit = "lm"
                          , Formula = "formula"
                          , TailFunction = "function"
                          , Triangle = "Triangle"
                          , SW = "htest"
                          , BP = "htest")
         # , sealed = TRUE
         # , validity = #some function
)

TailFunction = function(x, Tail)
{
  y = ifelse(x >= Tail, "Tail", x)
  if (!"Tail" %in% y){
    y[y == max(y)] = "Tail"
  }
  y
}

#' Create a new TriangleModel object
#' @export newTriangleModel
#' 
#' @include Triangle.R
#' 
#' @param Triangle A Triangle object
#' @param Response Character vector indicating the response being measured
#' @param Predictor Character vector indicating the variable used to predict the response
#' @param FitCategory Character vector indicating the column used to categorize the predictor variable
#' @param Intercept Boolean indicating whether or not to include an intercept
#' @param Alpha Numeric indicating the parameter used to weight the predictors
#' @param Tail Integer indicating the maximum development lag for grouping
#' 
newTriangleModel = function(Triangle
                            , Response
                            , Predictor
                            , FitCategory
                            , Intercept = FALSE
                            , Alpha = 0
                            , Tail = NULL)
{
  dfTriangleData = Triangle@TriangleData
  df = dfTriangleData[,c("OriginPeriod", "DevelopmentLag", "EvaluationDate", "DevInteger"
                         , "OriginPeriodStart", "OriginPeriodEnd", "CalendarPeriodStart"
                         , "CalendarPeriodEnd", "CalendarPeriod")]
  df = cbind(df, dfTriangleData[Response])
  df = cbind(df, dfTriangleData[Predictor])
  
  df = df[!is.na(df[Predictor]), ]
  
  if (length(FitCategory) > 1){
    # do something
    stop("Not yet configured for multiple groups.")
  }
  
  df$FitCategory = df[,FitCategory]
  
  if (is.null(Tail)) Tail = max(df$DevInteger) - 1
  
  if (FitCategory == "DevInteger"){
    df$FitCategory = TailFunction(df$FitCategory, Tail)
  }
  
  df$FitCategory = as.factor(df$FitCategory)
  
  strFormula = paste0(Response, " ~ ", Predictor, ":FitCategory")
  
  if (Intercept){
    strFormula = paste0(strFormula, " + 1:FitCategory")
  } else {
    strFormula = paste0(strFormula, " + 0")
  }
  
  weights = 1 / df[,Predictor] ^ (Alpha/2)
  
  theFormula = as.formula(strFormula)
  
  Fit = lm(theFormula, data = df, weights = weights)
  
  df$Residual = residuals(Fit)
  df$Predicted = predict.lm(Fit)
  
  SW = shapiro.test(df$Residual)
  BP = bptest(Fit)
  
  TriangleModel = new("TriangleModel"
                      , ModelData = df
                      , Response = Response
                      , Predictor = Predictor
                      , FitCategory = FitCategory
                      , Alpha = 0
                      , Tail = Tail
                      , Fit = Fit
                      , Formula = theFormula
                      , Triangle = Triangle
                      , SW = SW
                      , BP = BP)
  
  TriangleModel
}