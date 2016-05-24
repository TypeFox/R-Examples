#' USGS verification example: base parameteristion
#' @description USGS verification example
#' @references p.23, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
baseProps_ex1 <- function() {
  
  baseProperties(time = 1.5, ### days
                 basinWidth = 33.63, ### feet
                 basinLength = 33.63, ### feet
                 infiltrationRate = 1.333,
                 horizConductivity = 4, ### feet/day
                 iniHead = 10, ### feet
                 specificYield = 0.085,
                 numberTimeSteps = 150)
}

#' USGS verification example: model parameterisation (one distance)
#' @param x distance from the center of the recharge basin in the x direction (L)
#' (Default: 0 feet)
#' @param y distance from the center of the recharge basin in the y direction (L)
#' (Default: 0 feet)
#' @param baseProps basic model properties as retrieved by baseProps_ex1()
#' @param dbg If True additional messages on debug messages
#' @references p.23, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
#' @seealso \code{\link{baseProps_ex1}} for the USGS example parameterisation

conf_ex1 <- function(x = 0,
                      y = 0,
                      baseProps = baseProps_ex1(),
                      dbg=FALSE) {
  hantush(x = x,
          y = y,
          baseProps = baseProps,
          dbg = dbg
  )
}

#' USGS verification example: model parameterisation (multiple distances)
#' @param x vector of distances from the center of the recharge basin in the x
#' direction (L) (Default: 0 feet)
#' @param config name of configuration function (Default: conf_ex1)
#' @param dbg If True additional messages on debug messages
#' @param ... further arguments passed to function "config"
#' @references p.23, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
#' @examples
#' res <- confDistances_ex1()
#' if(FALSE){
#' #### Head for each time step (defined with parameter "numberTimeSteps)
#' xyplot(head ~ x | as.factor(sprintf("%f days", timeSteps)),
#'        data=res$timeSteps,
#'        type="b",
#'        las=1,
#'        as.table=TRUE)
#' #### Head at end of simulation
#' plot(head ~ x,
#'      data=res$simTime,
#'      type="b",
#'      las=1)
#'}
#' #### Water level increase at end of simulation & compare to alternative models
#' modelComparison <- compareModelResults(conf = res)
#' plotModelComparison(modelComparison = modelComparison)

confDistances_ex1 <- function(x = c(0, 0.3, 3.3, 6.6, 10, 20, 25, 30, 40, 50, 75, 100, 150, 200),
                               config= conf_ex1,
                               dbg = FALSE) {
  
  confFunction <- function(...) { config(...)}
  hantushDistances(x = x,
                   config = confFunction,
                   dbg = dbg)
  
}

#' USGS verification example: get model comparison table
#' @return data.frame with water level increase of different model approaches
#' @references Table 5, p.25, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
#' @examples
#' modelComparison <- getModelComparisonTable()

getModelComparisonTable <- function() {
  
  csvPath <- system.file("extdata",
                         "modelComparison.csv",
                         package = "kwb.hantush")
  modelComparison <- read.csv(file = csvPath, header = TRUE, dec = ".")
  
  return(modelComparison)
}


#' USGS verification example: compare R results for water level increase with other
#' models and calculate statitical goodness of fit values (e.g. RMSE, PBIAS, NSE)
#' @param conf list as retrieved by confDistances_ex1()
#' @return data.frame with R results and other model including all goodness of fit
#' criteria calculated with gof() of package hydrogof
#' @references Table 5, p.25, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
#' @seealso \code{\link{confDistances_ex1}} for the USGS example parameterisation with distances
compareModelResults <- function(conf = confDistances_ex1() ) {
  
  table <- getModelComparisonTable()
  
  res <- conf$simTime[,c("x", "head")]
  ### Water level increase
  res$head <- res$head - conf$baseProps$iniHead
  names(res)[names(res) == "head"] <- "R"
  
  modelComparison <-  merge(res, table)
  
  
  ##### Reformat
  newNames <- names(modelComparison)[3:ncol(modelComparison)]
  modelComparison <- reshape(data = modelComparison,times = newNames,
                             timevar = "model",
                             varying = list(newNames),
                             direction = "long")
  
  names(modelComparison)[names(modelComparison) == "Tailor"] <- "modelVal"
  
  gofRes <- data.frame()
  for (myModel in unique(modelComparison$model))
  {
    tmp <- modelComparison[modelComparison$model == myModel, ]
    tmp_gofRes <- t(hydroGOF::gof(sim = tmp$R,
                                  obs = tmp$modelVal,
                                  na.rm = TRUE,
                                  digits = 5
    )
    )
    colnames(tmp_gofRes) <- gsub(" %", "", colnames(tmp_gofRes))
    tmp_gofRes <- cbind(data.frame(model = myModel),as.data.frame(tmp_gofRes))
    
    gofRes <- rbind(gofRes, tmp_gofRes)
  }
  
  gofRes$ModelLabel <- sprintf("%s (RMSE: %1.3f feet, PBIAS: %2.1f %%, NSE: %1.3f)",
                               gofRes$model,
                               gofRes$RMSE,
                               gofRes$PBIAS,
                               gofRes$NSE)
  
  modelComparison <- merge(modelComparison, gofRes)
  
  return(modelComparison)
  
}

#' USGS verification example: plot model comparison results
#' @param modelComparison data.frame as retrieved by compareModelResults(),
#' Default: (compareModelResults())
#' @param title to be used as title above plot (Default: "")
#' @param ... further arguments passed to function xyplot()
#' @return model comparison goplots
#' @seealso \code{\link{compareModelResults}} for comparison with USGS benchmark models
#' @examples
#' ### Plot model comparison with title "Model comparison" and one plot for
#' ### each page
#' plotModelComparison(title = "Model comparison",
#'                     layout = c(1,1))

plotModelComparison <- function(modelComparison = compareModelResults(),
                                 title = "",
                                 ...) {
  
  print(lattice::xyplot(R + modelVal ~ x  | ModelLabel,
                        type = "b", pch = 16,
                        xlab = "x distance from the center of the recharge basin (feet) ",
                        ylab = "Water level increase (feet)",
                        data = modelComparison,
                        par.settings = lattice::simpleTheme(col = c("blue","red"),
                                                            pch = c(16,16), lty = c(1,1)),
                        auto.key = list(columns = 2, 
                                        text = c("R","Model"), 
                                        lines = TRUE, 
                                        points = FALSE),
                        main = title,
                        ...
  )
  )
}
