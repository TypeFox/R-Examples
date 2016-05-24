HasGroups = function(aTriangle){
  
  length(aTriangle@Groups) >= 2
}

OriginDevMatch = function(objTriangle){
  maxRows = length(unique(objTriangle@OriginPeriod))
  maxRows = maxRows * length(unique(objTriangle@DevelopmentLag))
  
  maxRows <= nrow(objTriangle@TriangleData)
}

TriangleOriginPeriods = function(objTriangle){
  op = unique(objTriangle@TriangleData$OriginPeriodStart)
  op = op[order(op)]
  op
}

TriangleCalendarPeriods = function(objTriangle){
  cp = unique(objTriangle@TriangleData$CalendarPeriodStart)
  cp = cp[order(cp)]
  cp
}

is.StochasticMeasure = function(objTriangle, MeasureName){
  MeasureName = CleanMeasureNames(MeasureName)
  MeasureName %in% objTriangle@StochasticMeasures
}

is.StaticMeasure = function(objTriangle, MeasureName){
  MeasureName = CleanMeasureNames(MeasureName)
  MeasureName %in% objTriangle@StaticMeasures
}

# TriangleGroups = function(aTriangle){
#   groups = unique(aTriangle@TriangleData$Group)
#   groups
# }