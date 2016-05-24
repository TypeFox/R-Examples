#' Load data and create an exreport experiment
#'
#' This function loads a data.frame, checks its properties and formats an
#' exreport experiment object. The columns of an experiments must contain at 
#' least two categorical columns to be identified as the method and problem 
#' variables and a thrid numerical column to be identified as an output variable. 
#' Additional columns can be added as parameters or additional outputs.
#'
#' @export
#' @param data A data.frame object satisfying the experiment format
#' @param methods The name of the variable which contains the methods, by default
#' is searches for a column named "method".
#' @param problems The name of the variable which contains the problems, by 
#' default is searches for a column named "problem".
#' @param parameters A list of the columns names to be identified as parameters. 
#' By default the remaining categorical columns 
#' are identified as parameters, so this list is useful only to identify numeric 
#' columns.
#' @param respectOrder A logical parameter which indicates if the order of the
#' elements of the method and problem columns must be respected by appearance or
#' ordered alphabeticaly. It affects to the look of data representations.
#' @param name A string which will identify the experiment in the report.
#' @param tol Tolerance factor to identify repeated experiments for duplicated 
#' rows.
#' @return A new exreport experiment object.
#' 
#' @seealso expCreateFromTable
#' 
#' @examples
#' # Creates experiment specifying column names and the numerical variables that
#' # are parameters
#' 
#' expCreate(wekaExperiment, 
#' methods="method", 
#' problems="problem", 
#' parameters="fold", 
#' name="Test Experiment")
#' 
expCreate <- function(data, methods="method", problems="problem", parameters = c(), respectOrder = FALSE, name, tol = 1E-9) {
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.data.frame(data))
    stop(.callErrorMessage("wrongParameterError", "data", "data.frame"))
  if (!is.character(parameters) & !is.null(parameters))
    stop(.callErrorMessage("wrongParameterError", "parameter", "character array"))
  
  # Check that method and column variables are present
  if (!(methods %in% colnames(data)))
    stop(.callErrorMessage("variableNotPresentError", methods))
  if (!(problems %in% colnames(data)))
    stop(.callErrorMessage("variableNotPresentError", problems))
  
  # Check if specified the parameters exists
  for (var in parameters) {    
    if (!(var %in% colnames(data)))
      stop(.callErrorMessage("variableNotPresentError", var))
  }
  
  # Convert columns to factors
  if (respectOrder){
    data[,methods]  <- factor(data[,methods], levels = unique(data[,methods]))
    data[,problems] <- factor(data[,problems], levels = unique(data[,problems]))
  }
  else{
    data[,methods]  <- as.factor(data[,methods])
    data[,problems] <- as.factor(data[,problems])
  }
  
  # Store the column names for methods and problems for further reference
  method <- methods
  problem <- problems
  
  # Search and process parameters columns, remaining categorical columns 
  # (not problem or methods) and those specified by the input
  params <- character()
  for (column in colnames(data)) {
    if (column %in% c(method, problem))
      next
    if (column %in% parameters | 
          is.factor(data[,column]) | 
          is.character(data[column])) {
      params <- c(params, column)
      # Convert it to a factor
      if (respectOrder)
        data[,column] <- factor(data[,column], levels = unique(data[,column]))
      else
        data[,column] <- as.factor(data[,column])
      next
    }    
  }
  
  # Check outputs columns. The remaining numerical columns are outputs
  outs <- character()
  for (column in colnames(data)) {
    if (column %in% c(method, problem))
      next
    if (column %in% params)
      next
    outs <- c(outs, column)    
  }
  
  # The problem must have at least one output
  if (length(outs) == 0)
    stop(.callErrorMessage("noOutputsError", var))
  
  # Create the experiments
  historic <- list(paste("Experiment ",
                         name ,
                         " loaded from a data set.",
                         sep=""))
  
  e  <- .experiment(data = data, 
                    method = method, 
                    problem = problem, 
                    params = params,
                    outs = outs, 
                    name = name, 
                    historic = historic)
  
  # Check if rows are duplicated, but do not revome them.
  dup <- expGetDuplicated(e,tol=tol)
  if(nrow(dup$data)>0)
    warning("you can get the duplicates through the function expGetDuplicated",call. = FALSE)  
  
  e
}




#' Create an exreport experiment from a tabular representation
#'
#' Create an exreport experiment object from a tabular representation. 
#' The input data must be a table having methods as rows and problems as columns.
#' The values in such table correspond to a particular output.
#' The resulting experiment can be characterized with static parameters.
#'
#' @export
#' @param data Input tabular data satisfying the previous constraints.
#' @param output String indicating the name of the output that the table values 
#' represent.
#' @param name A string which will identify the experiment in the report.
#' @param parameters A list of strings containing the names and values for the 
#' static configuration 
#' of the algorithm. The name of each element of the list will correspond with 
#' the name of a parameter
#' and the element with the value asigned.
#' @param respectOrder A logical parameter which indicates if the order of the
#' elements of the method and problem columns must be respected by appearance or
#' ordered alphabeticaly. It affects to the look of data representations.
#' @return A new exreport experiment object.
#' 
#' @seealso expCreate
#' 
#' @examples
#' # We generate a data frame where the methods are rows and the problems columns 
#' # from the wekaExperiment problem. (This is only an example, normally you
#' # would prefer to load a proper experiment and process it.)
#' 
#' library(reshape2)
#' df <- dcast(wekaExperiment[wekaExperiment$featureSelection=="no",], 
#' method ~ problem, 
#' value.var="accuracy", 
#' fun.aggregate = mean)
#' 
#' # We can create it and parametrice accordingly:
#' expCreateFromTable(df, output="accuracy", name="weka")
#' 
#' # Optionally we can set a fixed value for parameters, and ordered by appearance:
#' expCreateFromTable(df, output="accuracy", name="weka", 
#' parameters=list(featureSelection = "no"), respectOrder=TRUE)
#' 
expCreateFromTable <- function(data, output, name, parameters=list(), respectOrder = FALSE) {
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.data.frame(data))
    stop(.callErrorMessage("wrongParameterError", "data", "data.frame"))
  if (!is.character(output))
    stop(.callErrorMessage("wrongParameterError", "output", "character"))
  if (!is.list(parameters))
    stop(.callErrorMessage("wrongParameterError", "parameters", "list"))
  # Check the format of the parameter list
  if (length(parameters)!=0 & is.null(names(parameters)))
    stop(.callErrorMessage("noNamesError"))
  # Check that the output and of the parameters name is not problem or method or parameter: 
  if (output%in%names(parameters) || output=="method" || output=="problem")
    stop(.callErrorMessage("repeatedOutputError"))
  if (length(unique(names(parameters))) != length(names(parameters)) || "method"%in%names(parameters) || "problem"%in%names(parameters))
    stop(.callErrorMessage("repeatedParamsError"))
  
  # Melt the tabular data
  convertedData <- reshape2::melt(data, id.vars=1)
  
  # Rename columns, order of metl is fixed:
  names(convertedData) <- c("method", "problem", output)
  
  # Create the experiment:
  exp <- expCreate(convertedData, respectOrder = respectOrder, name = name)
  
  # Add fixed parameters:
  if (length(parameters)!=0)
    exp <- expExtend(exp, parameters)
  
  exp
}