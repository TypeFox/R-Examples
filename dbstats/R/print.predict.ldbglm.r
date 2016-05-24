
  #######################################
  #### print.predict.ldbglm function ####
  #######################################

print.predict.ldbglm<-function(x,...){
  cat("\n prediction for new data: \n")
  print(x$fit)
  cat("\n")
}