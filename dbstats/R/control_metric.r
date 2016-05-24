
 ########################
 #### control_metric ####
 ########################

 ## description: Internal function. Check if the metric is one we have defined.
 ##              If apply the default metric.
 ##
 ##        Inputs:  metric (character. 'metric' parameter, input by the user.)
 ##        Outputs: metric (returned the metric, if it's a defined metric. If not,
 ##                         apply the default metric 'euclidean') 
 ##


control_metric<-function(metric){
   m <- c("euclidean", "manhattan", "gower")
   
   # if character parameter: "metric" contains only part of the word, finish
   # filling the word, in accordance with the metrics m.
   metricaux <- m[pmatch(metric,m)]
   if (is.na(metricaux)){
     warning(gettextf("the metric %s is not defined. Will apply the default metric 'euclidean'",metric))
     metric<-"euclidean"
   }else
     metric<-metricaux
  
  return(metric)
}