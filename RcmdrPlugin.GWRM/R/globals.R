#Global variables need it.

#' @importFrom utils globalVariables
if (getRversion() >= '2.15.1')
  globalVariables(c('top','kvalue','buttonsFrame','lhsVariable', 'rhsVariable', 'subsetVariable', 'onHelp', 'xBox', 'outerOperatorsFrame', 
                    'formulaFrame','subsetFrame','lhsEntry','setBusyCursor','setIdleCursor'))