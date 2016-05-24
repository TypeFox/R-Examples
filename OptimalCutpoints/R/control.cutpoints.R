control.cutpoints <-
function(
	costs.ratio = 1,
  	CFP = 1,
  	CFN = 1,   
  	valueSp = 0.85,
  	valueSe = 0.85,  	
  	maxSp = TRUE,
  	generalized.Youden = FALSE,
  	costs.benefits.Youden = FALSE,
  	costs.benefits.Efficiency = FALSE,
  	weighted.Kappa = FALSE,
  	standard.deviation.accuracy = FALSE,
  	valueNPV = 0.85,
  	valuePPV = 0.85,
  	maxNPV = TRUE,
  	valueDLR.Positive = 2,
  	valueDLR.Negative = 0.5,
  	adjusted.pvalue = c("PADJMS","PALT5", "PALT10"),
  	ci.SeSp = c("Exact","Quadratic","Wald","AgrestiCoull","RubinSchenker"),
  	ci.PV = c("Exact","Quadratic","Wald","AgrestiCoull","RubinSchenker","Transformed","NotTransformed","GartNam"),
  	ci.DLR = c("Transformed","NotTransformed","GartNam"))
  	list(costs.ratio = costs.ratio, CFP = CFP , CFN = CFN, valueSp = valueSp, valueSe = valueSe, maxSp = maxSp, generalized.Youden = generalized.Youden, costs.benefits.Youden = costs.benefits.Youden, costs.benefits.Efficiency = costs.benefits.Efficiency, weighted.Kappa = weighted.Kappa, standard.deviation.accuracy = standard.deviation.accuracy, valueNPV = valueNPV, valuePPV = valuePPV, maxNPV = maxNPV, valueDLR.Positive = valueDLR.Positive, valueDLR.Negative = valueDLR.Negative, adjusted.pvalue = match.arg(adjusted.pvalue), ci.SeSp = match.arg(ci.SeSp), ci.PV = match.arg(ci.PV), ci.DLR = match.arg(ci.DLR))
