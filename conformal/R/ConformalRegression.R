## Isidro Cortes Ciriano. 2016-10-01
## Conformal Prediction for caret regression models


#############################################
### Conformal Prediction for Regression
#############################################

ConformalRegression <- setRefClass(
  "ConformalRegression",
  fields = list(
    PointPredictionModel = "ANY",
    ErrorModel = "ANY",
    confidence = "numeric",
    data.new = "ANY",
    alphas="numeric",
    errorPredictions = "numeric",
    pointPredictions = "numeric",
    Intervals = "numeric",
    plot = "ANY"
  ),
  methods = list(
    initialize = function(confi = 0.8)
    {
      "This method is called when you create an instance of the class."
      if (confi > 1 || confi < 0)
        stop("Confidence must be between 0 and 1")
      confidence <<- confi
      cat("Conformal Prediction Class for Regression Instantiated")
      cat("\n")
    },
    CalculateAlphas = function(model=NULL,error_model=NULL,
                               ConformityMeasure=StandardMeasure)
    {
      if(!is.function(StandardMeasure))
        stop("To calculate the alphas, a Conformity Measure need to be supplied.\n
             This functions must take three input vectors as arguments.")
      if(is.null(model) || is.null(error_model))
        stop("To calculate the alphas, a point prediction model and an error model 
             need to be suppplied")
      if(model$modelType != "Regression" || error_model$modelType != "Regression")
        stop("Both the point prediction and the error model need to be of type = 'Regression'")
      
      PointPredictionModel <<- model
      ErrorModel <<- error_model
      print("Calculating alphas..")
      cat('\n')
      predObsCV <- GetCVPreds(model)[,1:3]
      predObsCV <- predObsCV[order(predObsCV$rowIndex),]
      predObsCV_error <- GetCVPreds(error_model)[,1:3]
      predObsCV_error <- predObsCV_error[order(predObsCV_error$rowIndex),]
      
      out <- StandardMeasure(predObsCV$obs,predObsCV$pred,predObsCV_error$pred)
      alphas <<- out
    },
    GetConfidenceIntervals = function(new.data=NULL)
    {
      if (is.null(new.data)){
        stop("\nArgument 'data.new' cannot be empty.\nNew datapoints are required as input\n")
      }
      else{
        data.new <<- new.data
      }
      #require("caret") || stop("Pacakge 'caret' is required")
      
      print("Predicting (i) the value, and (ii) the error for the new data..")
      cat('\n')
      pred <- as.vector(predict(PointPredictionModel$finalModel, newdata = new.data))
      pointPredictions <<- pred
      pred_error <- as.vector(predict(ErrorModel$finalModel, newdata = new.data)) 
      errorPredictions <<- pred_error
      ### this formula can be tailored
      out <- ((alphas[length(alphas)*confidence]) * pred_error)
      Intervals <<- out

    },
    CorrelationPlot = function(obs=NULL,pred=pointPredictions,intervals=Intervals,margin = NULL, main = "", ylab = "Predicted", 
                               xlab = "Observed", PointSize =3, ColMargin = "blue", ErrorBarCol= "red",
                               ErrorBarSize = 0.5, ErrorBarWidth = 0.5, ErrorBarPosition= "identity", 
                               ErrorBarStat = "identity",TextSize = 15, 
                               TitleSize = 18, XAxisSize = 18, YAxisSize = 18, TitleAxesSize = 18, ErrorBarAlpha=0.8,
                               tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLabX = 0, 
                               PointColor = "black", PointAlpha = 1, PointShape = 15, MarginWidth = 1)
      {
     if(!(is.vector(obs)))
       stop("You must provide the observed values for the new data.")
     if (length(obs) != length(pred)) 
       stop("Both vectors have to be of equal length")
     if (sum(is.na(obs))>0 | sum(is.na(pred)) >0 | sum(is.na(intervals)) >0)
       stop("One of the input vectors (obs,pred, or intervals) contain NAs")
     if(is.null(margin))
       margin <- ( abs(max(obs,na.rm=T)) - abs(min(obs,na.rm=T)) ) / 4
     
    print("This method generates a correlation plot for the observed against the predicted values for a set of datapoints, e.g. a test set.")
     
     
     Data <- data.frame(Observed = obs, Predicted = pred,Intervals=intervals)
     p <- ggplot(Data, aes(x = Observed, y = Predicted)) +
       geom_abline(slope = 1,intercept = margin/2, colour = ColMargin, size = MarginWidth) + 
       geom_abline(slope = 1, intercept = -(margin/2), colour = ColMargin, 
                   size = MarginWidth) +
       geom_point(size = PointSize, colour = PointColor, 
                  shape = PointShape, alpha = PointAlpha) + 
       geom_errorbar(mapping = aes(x= Observed, ymax= Predicted + Intervals,
                                   ymin= Predicted-Intervals), 
                     position = ErrorBarPosition, 
                     stat = ErrorBarStat, 
                     alpha = ErrorBarAlpha,
                     width = ErrorBarWidth, 
                     color = ErrorBarCol, 
                     size = ErrorBarSize)  + theme_bw() + ggtitle(main) + 
       ylab(ylab) + xlab(xlab) + 
       ylim(c(min(c(obs, pred)-Data$Intervals),max(c(obs, pred)+Data$Intervals))) + 
       xlim(c(min(c(obs, pred)-Data$Intervals), max(c(obs,pred)+Data$Intervals))) + 
       theme(text = element_text(size = TextSize),
             axis.text.x = element_text(size = XAxisSize, angle = AngleLabX),
             axis.title.y = element_text(size = TitleAxesSize),
             axis.text.y = element_text(size = YAxisSize), 
             plot.title = element_text(size = TitleSize), 
             legend.key = element_blank(),
             plot.margin = unit(c(tmar, rmar, bmar, lmar), "cm")) 
     
     plot <<- p
    }
  )
)

