#' R6 class for Analyzing Fitbit  Data
#'
#' FitAnalyzer is an R6 class for analyzing Fitbit data. It is an opinionated implementation of a particular workflow for analysis. 
#' For people attempting to conduct their own analysis in a different fashion you should use the more generic functions implemented in FitUtil. \cr \cr
#' The workflow implemented for FitAnalyzer is the following: \cr
#' 1.	Create the FitAnalyzer with the goal variable for analysis. Eg: Calories or steps or distance. The goal variable is your personal goal that you are trying to analyze better. \cr
#' 2.	Call \code{findImportantVariables} to understand the most important variables unique to you that enable meeting your goal. \cr
#' 3.	Call \code{showMostImportantCharts} to get relevant charts that are unique to your data \cr
#' 4.	Call \code{predictGoal} to get a prediction on performance of the goal \cr \cr
#' You can conduct two types of analysis based on the type of dataset in consideration. \code{analysis.type} can be 'intra.day' or 'daily' analysis. 
#' 
#' @docType class
#' @format A \code{\link{R6Class}} generator object
#' @keywords data
#' 
#' @importFrom R6 R6Class
#' @importFrom dplyr arrange
#' @importFrom caret varImp
#' @importFrom gbm gbm predict.gbm gbm.perf relative.influence
#' @export FitAnalyzer
#' 
#' @section Methods:
#' \describe{
#'   \item{\code{getAnalysisFrame(folder, analysis.type)}}{This method uses \code{analysis.type} as an argument to return a data.frame that is clean and augmented with additional features like weekend.}
#'   \item{\code{findImportantVariables(tsDataFrame, seed = 12345)}}{Finds the most important variables that are enabling meeting the goals for the person, by creating a `glm` model and ranking the variables based on the coefficients of the model.}
#'   \item{\code{getFit()}}{Returns the `glm` fit object.}
#'   \item{\code{showMostImportantCharts(tsDataFrame)}}{Plots charts for the most relevant goals, with actual data and moving average using \code{geom_smooth()}.
#'   \cr \code{tsDataFrame}: a data.frame containing the fitibit activities.}
#'   \item{\code{predictGoal(x)}}{Gives a prediction on the goal performance, based on `glm` (daily) or `gbm` (intraday).}
#' }
#' 

##
## Begin niraj9@ code
##

FitAnalyzer <- R6::R6Class (
    "FitAnalyzer",
    
    public = list (
        
        initialize = function (goal = "calories") {
            private$goal <- goal
        },
        
        # Get Analysis frame
        getAnalysisFrame = function (folder = NA, analysis.type) {
            private$folder <- folder
            private$analysis.type <- analysis.type
            master <- NULL

            if (analysis.type == "intra.day") {
              master <-
                    createIntraFrame(folder)
              master <-
                    augmentIntraData(master)
            } else {
                master <- 
                    createTsMasterFrame(folder)
                master <- 
                    markValidRows(master)
                master <-
                    master[master$valid == TRUE, ]
                master <- augmentData(master)
            }
            return (master)
        },
        
        # Find important variables

        findImportantVariables = function (tsDataFrame, seed = 12345) {
            set.seed(seed)
            if (!is.null(private$fit)){
                return (private$imp.vars)
            }

            ifelse(private$analysis.type == "intra.day",
                   private$createIntraFit(tsDataFrame),
                   private$createDailyFrameFit(tsDataFrame)
                   )

            return (private$imp.vars)
        },
        
        # Get fit
        getFit = function () {
            return (private$fit)
        },
        
##
## End niraj9@ code
##
        
                
##
## Begin lassence@ code
##
        
        # Plot most important charts 
        showMostImportantCharts = function (tsDataFrame) {
            
            # Intraday plot
            if (private$analysis.type == "intra.day") {
                # Get important variables 
                intra.vars <- names(sort(private$imp.vars, decreasing = TRUE))
                intra.vars <- intra.vars[grep('intra.', intra.vars)]
                # Plot chart for 4 most important variables
                buildChartIntra(data = tsDataFrame,
                                y.axes = intra.vars[1:4])
                
            # Day time series plot
            } else {
                buildChartDay(
                    data = tsDataFrame,
                    y.axes = unlist(private$imp.vars$name)[1:4])
            }
        },
        
##
## End lassence@ code
##


##
## Begin niraj9@ code
##

        # Predict goals
        predictGoal = function (x) {
            response <- NULL
            response <- 
                  ifelse (private$analysis.type == "intra.day",
                          gbm::predict.gbm(private$fit, newdata = x,
                                           n.trees = private$gbm.best.iter),
                          predict.glm(private$fit, 
                                      newdata = as.data.frame(x), 
                                      type = "response"))

            return (response)
        }
        
    ),
    
    # Private variables
    private = list (
        
        folder = NULL,
        goal = NULL,
        imp.vars = NULL,
        analysis.type  = NULL,
        fit = NULL,
        gbm.best.iter = NULL,
        
        createDailyFrameFit = function (master) {
            y <-
                createGoalVariableVector(master, goal = private$goal)
            x <-
                createDependentVariableFrame(master, goal = private$goal)
            glm.fit <-
                glm(y ~ ., data = x, family = "gaussian")
            imp <- caret::varImp(glm.fit, scale = TRUE)
            imp$name <- rownames(imp)
            imp <- dplyr::arrange(imp, -Overall)
            private$fit <- glm.fit
            private$imp.vars <- imp
        },
        
        createIntraFit = function (master, cv.folds) {
            master$date <- NULL
            
            gbm.txt <- paste("gbm::gbm(formula = " ,
                             private$goal ,
                             "~ .,data = master,
                              distribution = 'gaussian',
                              n.trees = 500,
                              shrinkage = .05,
                              interaction.depth = 5,
                              bag.fraction = .5,
                              train.fraction = .8,
                              verbose = FALSE)", sep = "")
            gbm.fit <- eval(parse(text = gbm.txt))
            private$fit <- gbm.fit
            private$gbm.best.iter <-
                gbm::gbm.perf(gbm.fit, method = "test", plot.it = FALSE)
            private$imp.vars <-
                gbm::relative.influence(gbm.fit, n.trees = 500, scale = TRUE)
            private$imp.vars <- sort(private$imp.vars, decreasing = TRUE)
        }
    )
)

##
## End niraj9@ code
##


