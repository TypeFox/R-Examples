## webgui.R
##   - Web-based User Interface to RGP  
##
## RGP - a GP system for R
## 2010-14 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

RGP_PORT <- 6011

RGP_COLORS <- list(
  RED = "#D70026FF",
  YELLOW = "#E57600FF",
  BLUE = "#056D8FFF",
  GREEN = "#3BC500FF",
  GRAY = "#A0A0A0FF",
  DARK_GRAY = "#404040FF"
)

RGP_RUN_STATES <- list(PAUSED = 1, RUNNING = 2, RESET = 3, STOP = 4)
RGP_WORKER_MESSAGES <- list(PROGRESS = 1, NEWBEST = 2, STATISTICS = 3, RESULT = 4, RESET = 5, ALERT = 6)
RGP_HISTORY_LENGTH <- 1000

bootstrapButton <- function (inputId, label, class = "", icon = "", style = "") {
  tags$button(id = inputId, type = "button", style = style, class = paste("btn action-button", class), tags$i(class = icon), label)
}

infoPanel <- function(content, title = "Info") {
  div(class = "alert alert-info alert-block",
    tags$button(class = "close", "data-dismiss" = "alert", HTML("&times;")),
    tags$h4(tags$i(class = "fa fa-info-circle"), title),
    p(content))
}

dataPanel <- tabPanel("Data", value = "dataPanel",
  div(class = "row-fluid",
    sidebarPanel(
      tags$legend("CSV File Upload"),
      fileInput("csvFile", "CSV File",
                accept = c("text/csv", "text/comma-separated-values", "text/plain", ".csv")),
      checkboxInput("header", "Header", TRUE),
      radioButtons("sep", "Separator",
                   c(Comma = ",",
                     Semicolon = ";",
                     Tab = "\t"),
                   "Comma"),
      radioButtons("quote", "Quote",
                   c(None="",
                     "Double Quote" = "\"",
                     "Single Quote" = "'"),
                     "Double Quote"),
      tags$legend(style = "padding-top: 24px", "Data Partitioning"),
      sliderInput("trainingDataShare", "Training Data Share", 
                  min = 0.1, max = 1, value = .5, step = .01),
      selectInput("trainingDataPosition", "Training Data Position",
                  choices = c("Random", "First Rows", "Middle Rows", "Last Rows")),
      numericInput("dataPositionRandomSeed", "Random Seed", 
                   min = 0, value = 1, step = 1)),
    mainPanel(
      infoPanel("Use the 'Data' panel to import and pre-process your data set. To get started, upload a data file in CSV format by using the controls to the left. You can control the paritioning into training and validation sets with the 'Data Partitioning' controls.", title = "Data Panel"),
      tabsetPanel(
        tabPanel("Table", dataTableOutput("dataTable")),
        tabPanel("Plot", plotOutput("dataPlot"))),
      uiOutput("dataPanelHelpUi"))))
 
objectivePanel <- tabPanel("Objective", value = "objectivePanel",
  div(class = "row-fluid",
    sidebarPanel(
      tags$legend("Search Objective"),
      selectInput("dependentVariable", "Dependent Variable",
                  choices = c("")),
      textInput("formula", "Formula"),
      textInput("buildingBlocks", "Building Blocks",
                value = 'c("+", "-", "*", "/", "sin", "cos", "exp", "log", "sqrt")'),
      selectInput("errorMeasure", "Error Measure", 
                  choices = c("SMSE", "SSSE", "RMSE", "SSE", "MAE")),
      checkboxInput("enableComplexityCriterion", "Enable Complexity Criterion",
                    value = TRUE)),
    mainPanel(
      infoPanel("Use the 'Objective' panel to define the model objective. Choose the dependent variable, i.e. the variable to predict, via the combo-box to the left.", title = "Objective Panel"),
      plotOutput("dependentVariablePlot"),
      selectInput("dependentVariablePlotAbscissa", "Abscissa",
                  choices = c("(Row Number)")))))

runPanel <- tabPanel("Run", value = "runPanel",
  div(class = "row-fluid",
    sidebarPanel(
      tags$legend("Run Control"),
      bootstrapButton("startRunButton", "Start Run", icon = "fa fa-play", class = "btn-primary btn-block"),
      bootstrapButton("pauseRunButton", "Pause Run", icon = "fa fa-pause", class = "btn-block"),
      bootstrapButton("resetRunButton", "Reset Run", icon = "fa fa-eject", class = "btn-danger btn-block", style = "margin-top: 18px"),
      tags$legend(style = "padding-top: 24px", "Run Parameters"),
      numericInput("mu", "Mu (Population Size)", 
                   min = 2, max = 1000, value = 100, step = 1),
      numericInput("lambda", "Lambda (Number of Children / Generation)", 
                   min = 2, max = 100, value = 50, step = 1),
      numericInput("nu", "Nu (Number of New Individuals / Generation)", 
                   min = 2, max = 100, value = 50, step = 1),
      sliderInput("crossoverProbability", "Crossover Probability", 
                  min = 0, max = 1, value = .5, step = .01),
      sliderInput("subtreeMutationProbability", "Subtree Mutation Probability Weight", 
                  min = 0, max = 1, value = 1, step = .01),
      sliderInput("functionMutationProbability", "Function Mutation Probability Weight", 
                  min = 0, max = 1, value = 0, step = .01),
      sliderInput("constantMutationProbability", "Constant Mutation Probability Weight", 
                  min = 0, max = 1, value = 0, step = .01),
      checkboxInput("enableAgeCriterion", "Enable Age Criterion",
                    value = TRUE),
      sliderInput("parentSelectionProbability", "Parent Selection Probability", 
                  min = 0, max = 1, value = 1, step = .01),
      selectInput("selectionFunction", "Selection Function", 
                  choices = c("Crowding Distance", "Hypervolume")),
      sliderInput("fitnessSubSamplingShare", "Fitness Subsampling Share", 
                  min = 0, max = 1, value = 1, step = .01),
      numericInput("randomSeed", "Random Seed", 
                   min = 0, value = 1, step = 1)),
    mainPanel(
      infoPanel("Start the model search run by pressing the 'Start Run' button to the left. You can monitor the run's progress visually with the tools below. Once solutions of satisfactory quality start to appear, press 'Pause Run' and change to the 'Results' panel to analyse the results in more detail. You can continue a paused run by pressing 'Start Run' again or start over by pressing 'Reset Run'.", title = "Run Panel"),
      tabsetPanel(
        tabPanel("Progress", plotOutput("progressPlot", height = 1000)), 
        tabPanel("Pareto Front", plotOutput("paretoPlot", height = 768)), 
        tabPanel("Best Solution", plotOutput("bestSolutionPlot"),
                                  tableOutput("bestSolutionTable")),
        tabPanel("Statistics", tableOutput("runStatisticsTable"))))))
 
resultsPanel <- tabPanel("Results", value = "resultsPanel",
  tags$head(tags$style(type = "text/css", "tfoot { display: none; }")), # hack to disable column filtering in result table
  tags$head(tags$script(src = "scripts/jquery.sparkline.min.js")),
  tags$head(tags$script(type = "text/javascript", HTML("$(function() { $.extend($.fn.dataTable.defaults, { 'fnDrawCallback': function(oSettings) { $('.solutionSparkline').sparkline('html', { type: 'line', tagValuesAttribute: 'trueyvalues', disableHiddenCheck: true, height: '40px', width: '200px', lineColor: '#D70026', fillColor: false, disableInteraction: false, spotColor: false, minSpotColor: false, maxSpotColor: false, composite: false, chartRangeMin: 0, chartRangeMax: 1 }); $('.solutionSparkline').sparkline('html', { type: 'line', tagValuesAttribute: 'indyvalues', disableHiddenCheck: true, lineColor: '#056D8F', fillColor: false, disableInteraction: false, spotColor: false, minSpotColor: false, maxSpotColor: false, composite: true, chartRangeMin: 0, chartRangeMax: 1 }); } }); });"))),
  div(class = "row-fluid",
    infoPanel("The 'Results' panel shows a table of the results of a model search run in paused state. Also available is a variable presence plot, that shows how often each variable occurs in the set of result models. To show results, start a run in the 'Run' panel, then press 'Pause Run'.", title = "Results Panel"),
    tabsetPanel(
      tabPanel("Result Pareto Front", 
               dataTableOutput("resultParetoFrontTable")),
      tabPanel("Variable Importance",
               plotOutput("resultVariableImportancePlot")))))

##' HTML for rgpui 
##'
##' HTML for rgpui, represented as a DOM.
##'
##' @examples
##' cat(as.character(rgpuiHtml))
##'
##' @export
##' @import rgp 
##' @import emoa 
##' @import shiny
rgpuiHtml <- bootstrapPage(
  tags$head(tags$link(rel = "stylesheet", href = "css/font-awesome.min.css")),
  div(class = "container-fluid",
    div(class = "row-fluid",
      headerPanel(list(img(src = "images/logo_rgp.png"), span("RGP", tags$small("Symbolic Regression UI"))),
                  windowTitle = "RGP")),
    div(class = "row-fluid",
      uiOutput("alertUi")),
    div(class = "row-fluid",
      tabsetPanel(
        dataPanel,
        objectivePanel,
        runPanel,
        resultsPanel,
        id = "rgpPanels")),
    tags$footer(style = "padding-top: 8px; padding-bottom: 24px", HTML("&copy; 2010-13 Oliver Flasch, "), a("rsymbolic.org", href = "http://rsymbolic.org"))))
 
workerProcessMain <- function() {
  serverConnection <- socketConnection(port = RGP_PORT, server = TRUE, open = "rwb", blocking = TRUE)
  stopProcess <- FALSE
  command <- list(op = RGP_RUN_STATES$PAUSED)
  population <- NULL
  runStatistics <- NULL
  alertList <- list()

  while (!stopProcess) {
    tryCatch({
      Sys.sleep(0.5)
      if (socketSelect(list(serverConnection), timeout = 0)) {
        command <- unserialize(serverConnection)
        message(paste("job received command: ", command)) # TODO
      }
      if (RGP_RUN_STATES$PAUSED == command$op) {
        # do nothing
      } else if (RGP_RUN_STATES$RUNNING == command$op) {
        runResult <- workerProcessRun(serverConnection, population, runStatistics, command$params)
        command <- runResult$command
        population <- runResult$population
        runStatistics <- runResult$runStatistics
      } else if (RGP_RUN_STATES$RESET == command$op) {
        # reset worker process state...
        population <- NULL
        runStatistics <- NULL
        alertList <- list()
        serialize(list(msg = RGP_WORKER_MESSAGES$RESET,
                       params = list()),
                  serverConnection)
          command <- list(op = RGP_RUN_STATES$PAUSED) # change to PAUSED state
        } else if (RGP_RUN_STATES$STOP == command$op) {
        stopProcess <- TRUE
      } else {
        warning("RGP background job: unknown command: ", command)
        stopProcess <- TRUE
      }
    }, error = function(e) {
      message("workerProcessMain: Catched error '", e, "'")
      alertList <<- c(list(list(time = Sys.time(), type = "Error", content = e)), alertList)
      serialize(list(msg = RGP_WORKER_MESSAGES$ALERT,
                     params = list(alertList = alertList)),
                serverConnection)
      command <<- list(op = RGP_RUN_STATES$PAUSED) # change to PAUSED state
    })
  }

  close(serverConnection)
}

rescaleIndividual <- function(ind, srDataFrame, independentVariables, dependentVariable) {
  indX <- srDataFrame[, independentVariables]
  indY <- if (is.data.frame(indX)) apply(indX, 1, function(x) do.call(ind, as.list(x))) else ind(indX)
  trueY <- srDataFrame[, dependentVariable]
  indY <- if (length(indY) == 1) rep(indY, length(trueY)) else indY
  b = cov(trueY, indY) / var(indY)
  a = mean(trueY) - b * mean(indY)
  rescaledInd <- function(...) a + b * ind(...)
  return (rescaledInd)
}

errorMeasureFromName <- function(errorMeasureName) {
  switch(errorMeasureName,
         "SMSE" = smse,
         "SSSE" = ssse,
         "RMSE" = rmse,
         "SSE" = sse,
         "MAE" = mae,
         stop("rgpWebUi: unkown error measure name: ", errorMeasureName))
}

workerProcessRun <- function(serverConnection, population, runStatistics, params) {
  # check params for problems...
  if (is.null(params$dataFrame)) stop("rgpWebUi: No valid input data.")

  command <- list(op = RGP_RUN_STATES$RUNNING)

  serverState <- params$serverState
  set.seed(params$randomSeed) # TODO do not set seed when continuing a paused run, instead load it from the run state

  srFormula <- as.formula(params$formulaText)
  srDataFrame <- params$dataFrame$trainingData

  tryCatch({
    funSet <- do.call(functionSet, as.list(eval(parse(text = params$buildingBlocks))))
  }, error = function(e) { stop("rgpWebUi: Invalid building block set: '" , params$buildingBlocks, "'.") })
  inVarSet <- do.call(inputVariableSet, as.list(params$independentVariables))
  constSet <- numericConstantSet

  mutationFunction <- if (params$subtreeMutationProbability == 1 && params$functionMutationProbability == 0 && params$constantMutationProbability == 0) {
    function(ind) {
      subtreeMutantBody <- mutateSubtreeFast(body(ind), funSet, inVarSet, -10.0, 10.0, insertprob = 0.5, deleteprob = 0.5, subtreeprob = 1.0, constprob = 0.5, maxsubtreedepth = 8)
      makeClosure(subtreeMutantBody, inVarSet$all, envir = funSet$envir)
    }
  } else if (params$subtreeMutationProbability == 0 && params$functionMutationProbability == 1 && params$constantMutationProbability == 0) {
    function(ind) {
      functionMutantBody <- mutateFuncFast(body(ind), funSet, mutatefuncprob = 0.1)
      makeClosure(functionMutantBody, inVarSet$all, envir = funSet$envir)
    }
  } else if (params$subtreeMutationProbability == 0 && params$functionMutationProbability == 0 && params$constantMutationProbability == 1) {
    function(ind) {
      constantMutantBody <- mutateNumericConstFast(body(ind), mutateconstprob = 0.1, mu = 0.0, sigma = 1.0)
      makeClosure(constantMutantBody, inVarSet$all, envir = funSet$envir)
    }
  } else {
    function(ind) {
      mutantBody <- body(ind)
      weightSum <- params$subtreeMutationProbability + params$functionMutationProbability + params$constantMutationProbability
      rouletteWheelPosition <- runif(1, min = 0, max = weightSum)
      if (0 == weightSum) {
        return (ind)
      } else if (rouletteWheelPosition < params$subtreeMutationProbability) {
        mutantBody <- mutateSubtreeFast(mutantBody, funSet, inVarSet, -10.0, 10.0, insertprob = 0.5, deleteprob = 0.5, subtreeprob = 1.0, constprob = 0.5, maxsubtreedepth = 8)
      } else if (rouletteWheelPosition < params$subtreeMutationProbability + params$functionMutationProbability) {
        mutantBody <- mutateFuncFast(mutantBody, funSet, mutatefuncprob = 0.1)
      } else if (rouletteWheelPosition <= weightSum) {
        mutantBody <- mutateNumericConstFast(mutantBody, mutateconstprob = 0.1, mu = 0.0, sigma = 1.0)
      }
      makeClosure(mutantBody, inVarSet$all, envir = funSet$envir)
    }
  }

  errorMeasure  <- errorMeasureFromName(params$errorMeasure)

  ndsSelectionFunction <- switch(params$selectionFunction,
                                 "Crowding Distance" = nds_cd_selection,
                                 "Hypervolume" = nds_hv_selection,
                                 stop("rgpWebUi: unkown NDS selection function name: ", params$selectionFunction))

  searchHeuristic <- makeAgeFitnessComplexityParetoGpSearchHeuristic(lambda = params$lambda,
                                                                     crossoverProbability = params$crossoverProbability,
                                                                     newIndividualsPerGeneration = params$nu,
                                                                     enableComplexityCriterion = params$enableComplexityCriterion,
                                                                     enableAgeCriterion = params$enableAgeCriterion,
                                                                     ndsParentSelectionProbability = params$parentSelectionProbability,
                                                                     ndsSelectionFunction = ndsSelectionFunction)

  commandStopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) {
    command$op == RGP_RUN_STATES$PAUSED
  }

  runStatistics <- if (is.null(runStatistics)) {
    list(startTime = Sys.time(),
         fitnessHistory = c(),
         complexityHistory = c(),
         ageHistory = c(),
         stepNumber = 0,
         evaluationNumber = 0,
         timeElapsed = 0,
         bestFitness = Inf,
         dominatedHypervolumeHistory = c(),
         generationsHistory = c(),
         lastBestFitness = Inf)
  } else {
    runStatistics
  }

  currentRunStepNumber <- 0
  currentRunEvaluationNumber <-0
  currentRunTimeElapsed <- 0
  currentBestFitness <- Inf

  progressMonitor <- function(pop, objectiveVectors, fitnessFunction,
                              stepNumber, evaluationNumber, bestFitness, timeElapsed, indicesToRemove) {
    currentRunStepNumber <<- stepNumber
    currentRunEvaluationNumber <<- evaluationNumber
    currentRunTimeElapsed <<- timeElapsed
    currentBestFitness <<- currentBestFitness

    # offset step, time and evaluation counters with counters of previous run...
    stepNumber <- stepNumber + runStatistics$stepNumber
    evaluationNumber <- evaluationNumber + runStatistics$evaluationNumber
    timeElapsed <- timeElapsed + runStatistics$timeElapsed
    bestFitness <- min(bestFitness, runStatistics$bestFitness)

    # TODO maybe do not do this in every step
    if (socketSelect(list(serverConnection), timeout = 0)) {
      command <<- unserialize(serverConnection)
      message(paste("job received command: ", command)) # TODO
    }

    if (bestFitness < runStatistics$lastBestFitness) {
      alarm() # beep when a new best individual is found TODO transmit alarm to web client
      runStatistics$lastBestFitness <<- bestFitness
      bestIndividual <- pop[order(objectiveVectors$fitnessValues)][[1]] 
      rescaledBestIndividual <- if (params$errorMeasure == "SMSE" || params$errorMeasure == "SSSE") {
        rescaleIndividual(bestIndividual, srDataFrame, params$independentVariables, params$dependentVariable)
      } else {
        bestIndividual
      }
      message("NEW best solution (not rescaled):")
      message(sprintf(" %s", deparse(bestIndividual)))
      serialize(list(msg = RGP_WORKER_MESSAGES$NEWBEST,
                     params = list(bestIndividual = bestIndividual,
                                   rescaledBestIndividual = rescaledBestIndividual,
                                   stepNumber = stepNumber,
                                   evaluationNumber = evaluationNumber,
                                   bestFitness = bestFitness,
                                   timeElapsed = timeElapsed)),
                serverConnection)
    }

    if (stepNumber %% 10 == 0) { # every 10th generation...
      message(sprintf("evolution step %i, fitness evaluations: %i, best fitness: %f, time elapsed: %f",
                      stepNumber, evaluationNumber, bestFitness, timeElapsed))
      points <- rbind(objectiveVectors$fitnessValues, objectiveVectors$complexityValues, objectiveVectors$ageValues)
      finitePoints <- points[, !apply(is.infinite(points), 2, any)]
      bestFitnessIndex <- which.min(objectiveVectors$fitnessValues)
      historyIndexToReplace <- floor(runif(1, min = 1, max = RGP_HISTORY_LENGTH))
      runStatistics$fitnessHistory <<- c(if (length(runStatistics$fitnessHistory) <= RGP_HISTORY_LENGTH) runStatistics$fitnessHistory else runStatistics$fitnessHistory[-historyIndexToReplace], log(objectiveVectors$fitnessValues[bestFitnessIndex])) # cut-off at RGP_HISTORY_LENGTH number of entries
      runStatistics$complexityHistory <<- c(if (length(runStatistics$complexityHistory) <= RGP_HISTORY_LENGTH) runStatistics$complexityHistory else runStatistics$complexityHistory[-historyIndexToReplace], objectiveVectors$complexityValues[bestFitnessIndex])
      runStatistics$ageHistory <<- c(if (length(runStatistics$ageHistory) <= RGP_HISTORY_LENGTH) runStatistics$ageHistory else runStatistics$ageHistory[-historyIndexToReplace], objectiveVectors$ageValues[bestFitnessIndex])
      runStatistics$dominatedHypervolumeHistory <<- c(if (length(runStatistics$dominatedHypervolumeHistory) <= RGP_HISTORY_LENGTH) runStatistics$dominatedHypervolumeHistory else runStatistics$dominatedHypervolumeHistory[-historyIndexToReplace], dominated_hypervolume(finitePoints))
      runStatistics$generationsHistory <<- c(if (length(runStatistics$generationsHistory) <= RGP_HISTORY_LENGTH) runStatistics$generationsHistory else runStatistics$generationsHistory[-historyIndexToReplace], stepNumber)
      serialize(list(msg = RGP_WORKER_MESSAGES$PROGRESS,
                     params = list(stepNumber = stepNumber,
                                   evaluationNumber = evaluationNumber,
                                   timeElapsed = timeElapsed,
                                   generations = runStatistics$generationsHistory,
                                   fitnessHistory = runStatistics$fitnessHistory,
                                   complexityHistory = runStatistics$complexityHistory,
                                   ageHistory = runStatistics$ageHistory,
                                   dominatedHypervolumeHistory = runStatistics$dominatedHypervolumeHistory,
                                   poolFitnessValues = objectiveVectors$poolFitnessValues,
                                   poolComplexityValues = objectiveVectors$poolComplexityValues,
                                   poolAgeValues = objectiveVectors$poolAgeValues,
                                   poolIndicesToRemove = indicesToRemove)),
                serverConnection)
    }
  }

  # initialize population if NULL, otherwise reuse existing population...
  population <- if (is.null(population)) {
    message("workerProcessRun: INITIALIZING population")
    fastMakePopulation(params$mu, funSet, inVarSet, 8, -10.0, 10.0)
  } else {
    message("workerProcessRun: CONTINUING with existing population")
    population
  }

  message("workerProcessRun: STARTING GP run")
  sr <- suppressWarnings(symbolicRegression(srFormula,
                                            data = srDataFrame,
                                            functionSet = funSet,
                                            errorMeasure = errorMeasure,
                                            stopCondition = commandStopCondition,
                                            population = population,
                                            populationSize = params$mu,
                                            individualSizeLimit = 128, # individuals with more than 128 nodes (inner and leafs) get fitness Inf
                                            searchHeuristic = searchHeuristic,
                                            mutationFunction = mutationFunction,
                                            envir = environment(),
                                            verbose = FALSE,
                                            progressMonitor = progressMonitor))
  message("workerProcessRun: GP run done")

  # update runStatistics counters...
  runStatistics$stepNumber <- runStatistics$stepNumber + currentRunStepNumber 
  runStatistics$evaluationNumber <- runStatistics$evaluationNumber + currentRunEvaluationNumber
  runStatistics$timeElapsed <- runStatistics$timeElapsed + currentRunTimeElapsed 
  runStatistics$bestFitness <- currentBestFitness
 
  # send results to server
  serialize(list(msg = RGP_WORKER_MESSAGES$RESULT,
                 params = list(result = sr,
                               serverState = serverState,
                               data = srDataFrame,
                               validationData = params$dataFrame$validationData,
                               independentVariables = params$independentVariables,
                               dependentVariable = params$dependentVariable,
                               errorMeasure = params$errorMeasure)),
            serverConnection)

  return (list(command = command, result = sr, population = population, runStatistics = runStatistics))
}

server <- function(input, output, session) {
  serverState <- list(command = RGP_RUN_STATES$PAUSED, seed = NULL) 
  workerProcess <- mcparallel(workerProcessMain())
  Sys.sleep(1) # wait for background job to initialize 
  workerProcessConnection <- socketConnection(port = RGP_PORT, open = "rwb", blocking = TRUE) 

  independentVariables <- c()

  dataFrame <- reactive({
    dataFile <- input$csvFile

    if (is.null(dataFile))
      return (NULL)
    
    dataFrame <- read.csv(dataFile$datapath, header = input$header, sep = input$sep, quote = input$quote, colClasses = "numeric")

    updateSelectInput(session, "dependentVariable",
                      choices = colnames(dataFrame), selected = tail(colnames(dataFrame), 1))
    updateSelectInput(session, "dependentVariablePlotAbscissa",
                      choices = c("(Row Number)", colnames(dataFrame)), selected = "(Row Number)") 

    trainingDataIndices <- if ("Random" == input$trainingDataPosition) {
      set.seed(input$dataPositionRandomSeed)
      sample(1:nrow(dataFrame), size = floor(input$trainingDataShare * nrow(dataFrame)), replace = FALSE)
    } else if ("First Rows" == input$trainingDataPosition) {
      1:floor(input$trainingDataShare * nrow(dataFrame))      
    } else if ("Middle Rows" == input$trainingDataPosition) {
      skipLength <- floor((1 - input$trainingDataShare) * nrow(dataFrame) / 2)
      skipLength:(nrow(dataFrame) - skipLength) 
    } else if ("Last Rows" == input$trainingDataPosition) {
      floor((1 - input$trainingDataShare) * nrow(dataFrame)):nrow(dataFrame)
    } else {
      stop("rgpWebUi: Unknown training data position type: '", input$trainingDataPosition, "'.")
    }

    return (list(data = dataFrame,
                 trainingData = dataFrame[trainingDataIndices,],
                 validationData = dataFrame[-trainingDataIndices,],
                 trainingDataIndices = trainingDataIndices))
  })

  # generate symbolic regression formula from input data column headers
  observe({
    allVariables <- colnames(dataFrame()$data)
    dependentVariable <- input$dependentVariable 
    independentVariables <<- allVariables[allVariables != dependentVariable]
    formulaText <- paste(dependentVariable, "~", paste(independentVariables, collapse = " + "))
    updateTextInput(session, "formula", value = formulaText)
  })

  # extract independentVariables from symbolic regression formula
  observe({
    tryCatch({
      f <- as.formula(input$formula)
      independentVariables <<- attr(terms(f), "term.labels")
    }, error = function(e) NULL)
  })

  output$dataPanelHelpUi <- renderUI({
    if (is.null(dataFrame())) {
      div(class = "hero-unit",
          h1("Welcome to RGP"),
          p("The RGP Symbolic Regression UI makes it easy to discover mathematical models for your data."),
          p(style = "font-size: 14px; font-weight: normal", "Getting starting with Symbolic Regression in RGP in easy:",
          tags$ol(tags$li("Upload a comma-separated (CSV) data file via the controls to the left."),
                  tags$li("Set the model objective in the 'Objective' tab."),
                  tags$li("Start the Symbolic Regression run in the 'Run' tab."),
                  tags$li("Monitor the run's progress in the 'Progress' tab."),
                  tags$li("When solutions of sufficient quality emerge, pause the run in the 'Run' tab."),
                  tags$li("Inspect the solutions in the 'Results' tab.")),
          "Visit ", tags$a(href = "http://rsymbolic.org", "rsymbolic.org"), " for detailed information about the RGP system."))
    } else { "" }
  })
  
  output$dataPlot <- renderPlot({
    if (is.null(dataFrame())) {
      NULL
    } else {
      trainingDataMarkerColumn <- numeric(nrow(dataFrame()$data))
      trainingDataMarkerColumn <- rep(1, nrow(dataFrame()$data))
      trainingDataMarkerColumn[dataFrame()$trainingDataIndices] <- 16 
      plot(dataFrame()$data, col = RGP_COLORS$RED, pch = trainingDataMarkerColumn)
      legend("bottomright", c("Training", "Validation"),
             bty = "n",
             pch = c(16, 1),
             col = c(RGP_COLORS$RED, RGP_COLORS$RED))
    }
  })

  output$dataTable <- renderDataTable({
    if (is.null(dataFrame())) {
      NULL
    } else {
      trainingDataMarkerColumn <- rep("Validation", nrow(dataFrame()$data))
      trainingDataMarkerColumn[dataFrame()$trainingDataIndices] <- "Training"
      cbind(dataFrame()$data, list("Role" = trainingDataMarkerColumn))
    }
  })
  
  output$dependentVariablePlot <- renderPlot({
    if (is.null(dataFrame())) {
      NULL
    } else {
      trainingDataMarkerColumn <- numeric(nrow(dataFrame()$data))
      trainingDataMarkerColumn <- rep(1, nrow(dataFrame()$data))
      trainingDataMarkerColumn[dataFrame()$trainingDataIndices] <- 16 
      if ("(Row Number)" == input$dependentVariablePlotAbscissa) {
        plot(dataFrame()$data[, input$dependentVariable], col = RGP_COLORS$RED, pch = trainingDataMarkerColumn,
             xlab = "Row Number", ylab = input$dependentVariable,
             main = "Dependent Variable Plot")
        lines(dataFrame()$data[, input$dependentVariable], col = RGP_COLORS$GRAY)
      } else {
        plot(x = dataFrame()$data[, input$dependentVariablePlotAbscissa],
             y = dataFrame()$data[, input$dependentVariable],
             col = RGP_COLORS$RED, pch = trainingDataMarkerColumn,
             xlab = input$dependentVariablePlotAbscissa, ylab = input$dependentVariable,
             main = "Dependent Variable Plot")
        lines(x = dataFrame()$data[, input$dependentVariablePlotAbscissa],
              y = dataFrame()$data[, input$dependentVariable],
              col = RGP_COLORS$GRAY)
        legend("bottomright", c("Training", "Validation"),
               bty = "n",
               pch = c(16, 1),
               col = c(RGP_COLORS$RED, RGP_COLORS$RED))
      }
    }
  })
  
  observe({ if (input$startRunButton > 0) {
    serverState$command <<- RGP_RUN_STATES$RUNNING 
    serialize(list(op = serverState$command, params = list(serverState = serverState,
                                                           buildingBlocks = isolate(input$buildingBlocks),
                                                           independentVariables = independentVariables,
                                                           dependentVariable = isolate(input$dependentVariable),
                                                           mu = isolate(input$mu),
                                                           lambda = isolate(input$lambda),
                                                           nu = isolate(input$nu),
                                                           crossoverProbability = isolate(input$crossoverProbability),
                                                           subtreeMutationProbability = isolate(input$subtreeMutationProbability),
                                                           functionMutationProbability = isolate(input$functionMutationProbability),
                                                           constantMutationProbability = isolate(input$constantMutationProbability),
                                                           enableAgeCriterion = isolate(input$enableAgeCriterion),
                                                           enableComplexityCriterion = isolate(input$enableComplexityCriterion),
                                                           parentSelectionProbability = isolate(input$parentSelectionProbability),
                                                           selectionFunction = isolate(input$selectionFunction),
                                                           fitnessSubSamplingShare = isolate(input$fitnessSubSamplingShare),
                                                           errorMeasure = isolate(input$errorMeasure),
                                                           randomSeed = isolate(input$randomSeed),
                                                           formulaText = isolate(input$formula),
                                                           dataFrame = isolate(dataFrame()))), workerProcessConnection)
  }})

  observe({ if (input$pauseRunButton > 0) {
    serverState$command <<- RGP_RUN_STATES$PAUSED
    serialize(list(op = serverState$command), workerProcessConnection) 
  }})

  observe({ if (input$resetRunButton > 0) {
    serverState$command <<- RGP_RUN_STATES$RESET
    # reset local state
    serverState$seed <- NULL
    # send reset command to worker process
    serialize(list(op = serverState$command), workerProcessConnection) 
    # change state to PAUSED
    serverState$command <<- RGP_RUN_STATES$PAUSED
  }})

  workerProcessMessage <- reactive({
    invalidateLater(100, session) # each 100 milliseconds
    return (if (socketSelect(list(workerProcessConnection), timeout = 0)) {
      unserialize(workerProcessConnection)
    } else {
      NULL
    })
  })

  lastWorkerProcessMessages <- reactiveValues(progress = NULL, newBest = NULL, statistics = NULL, result = NULL, alert = NULL)

  observe({
    msg <- workerProcessMessage() # make sure that unserialize(workerProcessConnection) keeps called regularly
    if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$PROGRESS) {
      lastWorkerProcessMessages$progress <- msg
    } else if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$NEWBEST) {
      lastWorkerProcessMessages$newBest <- msg
    } else if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$STATISTICS) {
      lastWorkerProcessMessages$statistics <- msg
    } else if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$RESULT) {
      lastWorkerProcessMessages$result <- msg
    } else if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$RESET) {
      lastWorkerProcessMessages$progress <- NULL
      lastWorkerProcessMessages$newBest <- NULL
      lastWorkerProcessMessages$statistics <- NULL
      lastWorkerProcessMessages$result <- NULL
      lastWorkerProcessMessages$alert <- NULL
    } else if (!is.null(msg) && msg$msg == RGP_WORKER_MESSAGES$ALERT) {
      lastWorkerProcessMessages$alert <- msg
    } else if (!is.null(msg)) {
      stop("rgpWebUi: unknown worker process message:", msg)
    } else {
      # ignore NULL messages
    }
  })
  
  output$alertUi <- renderUI({
    if (!is.null(lastWorkerProcessMessages$alert)) {
      alertList <- lastWorkerProcessMessages$alert$params$alertList
      do.call(tagList, Map(function(alert) {
        div(class = "alert alert-error alert-block",
          tags$button(class = "close", "data-dismiss" = "alert", HTML("&times;")),
          tags$h4(tags$i(class = "fa fa-warning"), alert$type, paste("(", alert$time, ")", sep = "")),
          HTML(alert$content))
      }, alertList))
    }
  })

  output$progressPlot <- renderPlot({
    if (!is.null(lastWorkerProcessMessages$progress)) {
      params <- lastWorkerProcessMessages$progress$params 
      oldPar <- par(no.readonly = TRUE)
      layout(matrix(1:4, 4, 1, byrow = TRUE))
      plot(params$generations, params$fitnessHistory, type = "l",
           main = "Fittest Individual Fitness", xlab = "Generation", ylab = "log Fitness")
      plot(params$generations, params$complexityHistory, type = "l", col = "red",
           main = "Fittest Individual Complexity", xlab = "Generation", ylab = "Complexity (Visitation Length)")
      plot(params$generations, params$ageHistory, type = "l", col = "green",
           main = "Fittest Individual Age", xlab = "Generation", ylab = "Age (Generations)")
      plot(params$generations, params$dominatedHypervolumeHistory, type = "l", col = "gray",
           main = "Dominated Hypervolume", xlab = "Generation", ylab = "Hypervolume")
      par(oldPar)
    }
  })

  output$paretoPlot <- renderPlot({
    if (!is.null(lastWorkerProcessMessages$progress)) {
      params <- lastWorkerProcessMessages$progress$params 
      plotParetoFront(params$poolFitnessValues, params$poolComplexityValues, params$poolAgeValues,
                      params$poolIndicesToRemove,
                      main = sprintf("Selection Pool Fitness Pareto Plot (%d Individuals)", length(params$poolFitnessValues)),
                      xlab = "Fitness (Prediction Error)", ylab = "Complexity (Visitation Length)")
    }
  })

  output$bestSolutionPlot <- renderPlot({
    if (!is.null(lastWorkerProcessMessages$newBest)) {
      trainingDataMarkerColumn <- numeric(nrow(dataFrame()$data))
      trainingDataMarkerColumn <- rep(1, nrow(dataFrame()$data))
      trainingDataMarkerColumn[dataFrame()$trainingDataIndices] <- 16
      params <- lastWorkerProcessMessages$newBest$params 
      ind <- params$rescaledBestIndividual
      indX <- dataFrame()$data[, independentVariables]
      indY <- if (is.data.frame(indX)) apply(indX, 1, function(x) do.call(ind, as.list(x))) else ind(indX)
      plot(dataFrame()$data[, input$dependentVariable], col = RGP_COLORS$RED, pch = trainingDataMarkerColumn,
           xlab = "Row Number", ylab = input$dependentVariable,
           main = "Best Solution Plot")
      lines(dataFrame()$data[, input$dependentVariable], col = RGP_COLORS$GRAY)
      points(indY, col = RGP_COLORS$BLUE, pch = trainingDataMarkerColumn - 1)
      lines(indY, col = RGP_COLORS$DARK_GRAY)
      legend("bottomright", c("Data", "Best Solution"),
             bty = "n",
             pch = c(16, 15),
             col = c(RGP_COLORS$RED, RGP_COLORS$BLUE))
    }
  })

  output$bestSolutionTable <- renderTable({
    if (!is.null(lastWorkerProcessMessages$newBest)) {
      params <- lastWorkerProcessMessages$newBest$params 
      data.frame(list(Attribute = c("Formula", "Error", "Generation", "Fitness Evaluation Number", "Time Elapsed"),
                      Value = c(do.call(paste, c(as.list(deparse(params$bestIndividual)), sep = "")),
                                params$bestFitness,
                                params$stepNumber,
                                params$evaluationNumber,
                                formatSeconds(params$timeElapsed))))
    }
  }, include.rownames = FALSE)
  
  output$runStatisticsTable <- renderTable({
    if (!is.null(lastWorkerProcessMessages$progress)) {
      params <- lastWorkerProcessMessages$progress$params 
      data.frame(list(Attribute = c("Generation", "Fitness Evalutation Number", "Time Elapsed", "Fitness Evaluations / Second"),
                      Value = c(params$stepNumber,
                                params$evaluationNumber,
                                formatSeconds(params$timeElapsed),
                                sprintf("%.2f", params$evaluationNumber / params$timeElapsed))))
    }
  }, include.rownames = FALSE)

  output$resultParetoFrontTable <- renderDataTable({
    if (!is.null(lastWorkerProcessMessages$result)) {
      params <- lastWorkerProcessMessages$result$params
      srResult <- params$result
      population <- srResult$population
      fitnessValues <- srResult$searchHeuristicResults$fitnessValues
      complexityValues <- srResult$searchHeuristicResults$complexityValues
      objectiveValues <- rbind(fitnessValues, complexityValues)
      ndsRanks <- nds_rank(objectiveValues)
      paretoFrontMask <- ndsRanks == 1
      uniqueMask <- !duplicated(population)
      mask <- paretoFrontMask & uniqueMask 
      deparseInd <- function(ind) do.call(paste, c(as.list(deparse(ind)), sep = ""))
      indToYstring <- function(ind) {
        rescaledInd <- if (params$errorMeasure == "SMSE" || params$errorMeasure == "SSSE") {
          rescaleIndividual(ind, dataFrame()$data, params$independentVariables, params$dependentVariable)
        } else {
          ind 
        }
        indX <- dataFrame()$data[, independentVariables]
        indY <- if (is.data.frame(indX)) apply(indX, 1, function(x) do.call(rescaledInd, as.list(x))) else rescaledInd(indX)
        trueY <- dataFrame()$data[, params$dependentVariable]
        totalMin <- min(c(indY, trueY)) # normalize plot ranges
        totalRange <- max(c(indY, trueY)) - totalMin
        indYNormalized <- (indY - totalMin) / totalRange
        trueYNormalized <- (trueY - totalMin) / totalRange
        indYString <- do.call(paste, c(as.list(indYNormalized), sep = ","))
        trueYString <- do.call(paste, c(as.list(trueYNormalized), sep = ","))
        HTML(paste("<span class='solutionSparkline' indyvalues='", indYString, "' trueyvalues='", trueYString, "'></span>", sep = ""))
      }
      paretoFrontFormulas <- as.character(Map(deparseInd, population[mask]))
      paretoFrontFitnessValues <- fitnessValues[mask]
      errorMeasure <- errorMeasureFromName(params$errorMeasure)
      trueY <- params$validationData[, params$dependentVariable]
      paretoFrontIndYs <- Map(function(ind) {
        rescaledInd <- if (params$errorMeasure == "SMSE" || params$errorMeasure == "SSSE") {
          rescaleIndividual(ind, params$validationData, params$independentVariables, params$dependentVariable)
        } else {
          ind 
        }
        indX <- params$validationData[, independentVariables]
        indY <- if (is.data.frame(indX)) apply(indX, 1, function(x) do.call(rescaledInd, as.list(x))) else rescaledInd(indX)
        return (if (length(indY) == 1) rep(indY, length(trueY)) else indY)
      }, population[mask])
      paretoFrontValidationFitnessValues <- as.numeric(Map(function(indY) errorMeasure(trueY, indY), paretoFrontIndYs))
      paretoFrontValidationRSquaredValues <- as.numeric(Map(function(indY) rsquared(trueY, indY), paretoFrontIndYs))
      paretoFrontComplexityValues <- complexityValues[mask]
      paretoFrontPlots <- as.character(Map(indToYstring, population[mask]))

      result <- data.frame(list(Formula = paretoFrontFormulas,
                                TrainingError = paretoFrontFitnessValues,
                                ValidationError = paretoFrontValidationFitnessValues,
                                ValidationRSquared = paretoFrontValidationRSquaredValues,
                                Complexity = paretoFrontComplexityValues,
                                Plot = paretoFrontPlots)) 
      colnames(result) <- c("Formula", "Error (Training)", "Error (Validation)", "R^2 (Validation)", "Complexity", "Plot")

      return (result)
    }
  }, options = list(pageLength = 25, searching = 0, info = 0, columns = list(list(orderable = 0), NULL, NULL, NULL, NULL, list(orderable = 0)))) # TODO optimize data table options

  output$resultVariableImportancePlot <- renderPlot({
    if (!is.null(lastWorkerProcessMessages$result)) {
      params <- lastWorkerProcessMessages$result$params
      srResult <- params$result
      population <- srResult$population
      fitnessValues <- srResult$searchHeuristicResults$fitnessValues
      complexityValues <- srResult$searchHeuristicResults$complexityValues
      objectiveValues <- rbind(fitnessValues, complexityValues)
      ndsRanks <- nds_rank(objectiveValues)
      paretoFrontMask <- ndsRanks == 1
      paretoFront <- population[paretoFrontMask]
      variablePresenceParetoFront <- colSums(populationVariablePresenceMap(paretoFront))
      variablePresencePopulation <- colSums(populationVariablePresenceMap(population)) - variablePresenceParetoFront
      variablePresence <- rbind(variablePresenceParetoFront, variablePresencePopulation)
      barplot(variablePresence, main = "Independent Variable Presence",
              xlab = "Variable", ylab = "Count", col = c(RGP_COLORS$RED, RGP_COLORS$BLUE), border = NA,
              legend = c("Pareto Front Models", "All Models of Result Population"), args.legend = list(bty = "n"))
    }
  })
}

##' Web-based user interface to RGP 
##'
##' Start the web-based user interface to RGP. Currently, this user interface only
##' support Symbolic Regression.
##'
##' @param port The TCP port the web-based user interface to listen on, defaults to
##'   port 1447.
##'
##' @examples
##' \dontrun{
##' symbolicRegressionUi()
##' }
##'
##' @export
##' @import rgp 
##' @import emoa 
##' @import shiny
##' @import parallel
symbolicRegressionUi <- function(port = 1447) {
  addResourcePath("css", system.file("css", package = "rgpui"))
  addResourcePath("images", system.file("images", package = "rgpui"))
  addResourcePath("fonts", system.file("fonts", package = "rgpui"))
  addResourcePath("scripts", system.file("scripts", package = "rgpui"))
  runApp(list(ui = rgpuiHtml, server = server), port = port)
}

.onAttach <- function(libname, pkgname) {
  # show startup message
  packageStartupMessage("*** RGPUI version ", (sessionInfo())$otherPkg$rgpui$Version, " initialized successfully.\n",
                        "    Type 'symbolicRegressionUi()' to start the web-based symbolic regression UI,\n",
                        "    then navigate to http://localhost:1447.")
}

