.errorStrings <- list(
  wrongParameterError       = "Error: %s must be a variable of type %s",
  variableNotPresentError   = "Error: Specified %s variable not present in data",
  noOutputsError            = "Error: At least one output variable must be present in the input data",
  noNamesError              = "Error: Please provide valid names for the variables input list",
  noValueNamesError         = "Error: Please provide valid names for the values input array",
  parametersDifferError     = "Error: The parameter of the given experiments differ",
  commonOutputsError        = "Error: There is no common output variables in the given experiments",
  failedIntegrityError      = "Error: The output of one or more duplicated rows differs, please check the experiment integrity",
  requireInstantiationError = "Error: There are parameters present in the experiment, please perform instantiation before",
  parameterInRangeError     = "Error: The parameter %s require values in %s",
  wilcoxonNMethodsError     = "Error: Wilcoxon Test require only two methods to be applied",
  friedman2MethodsError     = "Error: Only two methods available for Friedman Test, please perform Wilcoxon Test instead",
  reportableError           = "Error: Only Reportable objects can be added to the report",
  invalidFolderError        = "Error: Please provide a valid destination folder",
  repeatedParamsError       = "Error: Please provide an unique name for each parameter in the input list, different method and problem",
  repeatedOutputError       = "Error: Please provide an unique name for each output variable, different from parameters, method and problem",
  invalidOutputName         = "Error: There is at least one common name in %s and %s. Please, use different names",
  directoryAlreadyExists    = "Error: The directory %s already exists"
)

.warningStrings <- list(
  differentMethodName       = "Warning: The method name is different for %s and %s. The name \"%s\" will be used",
  differentProblemName      = "Warning: The problem name is different for %s and %s. The name \"%s\" will be used",
  outputsRemoved            = "Warning: Only the outputs '%s' will remain in the experiment",
  unmatchedInstancesRemoved = "There are unmatched instances between experiments. The combined experiment will omit them",
  outputsRenamed            = "There are common output names (%s). They will be renamed",
  varsInTemplateNotExist    = "The variables [%s] in the text have not been provided",
  friedmanTestNotRejected   = "Friedman test was not rejected with p-value %.4e >= %.4f alpha, computed p-values will not be significant."
)

.callErrorMessage <- function(error, ...)
{
  sprintf(.errorStrings[[error]], ...)
}

.callWarningMessage <- function(warn, ...)
{
  sprintf(.warningStrings[[warn]], ...)
}

.mapId2Name <- list(
  friedman          = "Friedman", #used in statistical_test, testFriedman function
  imanDaven         = "Iman Davenport", #used in statistical_test, testFriedman function
  min               = "minimize", #used in statistical_test, testFriedman function
  max               = "maximize", #used in statistical_test, testFriedman function
  pvalue            = ".phMetricPvalues", #used in web_reports, testPosthocTable function
  wtl               = ".phMetricWTL",  #used in web_reports, testPosthocTable function
  rank              = ".phMetricRanks"  #used in web_reports, testPosthocTable function
)