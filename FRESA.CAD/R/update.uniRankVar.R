#' @method update uniRankVar
update.uniRankVar <-
function(object,...)
{
	parameters <- list(...);
	if (is.null(parameters$variableList)) { variableList = object$variableList;} else {variableList = parameters$variableList;}
	if (is.null(parameters$formula)) { formula = object$formula;} else {formula = parameters$formula;}
	if (is.null(parameters$Outcome)) { Outcome = object$Outcome;} else {Outcome = parameters$Outcome;}
	if (is.null(parameters$data)) {data = object$data;} else {data = parameters$data;}
	if (is.null(parameters$categorizationType)) { categorizationType = object$categorizationType;} else {categorizationType = parameters$categorizationType;}
	if (is.null(parameters$type)) { type = object$type;} else {type = parameters$type;}
	if (is.null(parameters$rankingTest)) { rankingTest = object$rankingTest;} else {rankingTest = parameters$rankingTest;}
	if (is.null(parameters$cateGroups)) { cateGroups = object$cateGroups;} else {cateGroups = parameters$cateGroups;}
	if (is.null(parameters$raw.dataFrame)) { raw.dataFrame = object$raw.dataFrame;} else {raw.dataFrame = parameters$raw.dataFrame;}
	if (is.null(parameters$description)) { description = object$description;} else {description = parameters$description;}
	if (is.null(parameters$uniType)) { uniType = object$uniType;} else {uniType = parameters$uniType;}
	if (is.null(parameters$FullAnalysis)) { FullAnalysis = TRUE} else {FullAnalysis = parameters$FullAnalysis;}
	
	
	updatedRank <- uniRankVar(variableList,formula,Outcome,data,categorizationType,type,rankingTest,cateGroups,raw.dataFrame,description,uniType,FullAnalysis) 

	return (updatedRank);
	
}
