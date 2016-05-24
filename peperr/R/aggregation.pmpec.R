`aggregation.pmpec` <-
function(full.data, response, x, model, cplx=NULL, times=NULL, type=c("apparent", "noinf"), 
   fullsample.attr=NULL, ...){ 
   require(survival)   

   xnames <- names(x) 
   time <- response[,"time"]
   status <- response[,"status"]
   
   uncens <- which(status == 1)
   time.uncens <- time[uncens]
   if (!is.null(fullsample.attr)){
      times <- fullsample.attr
   }
   if (is.null(times)){
      if(length(unique(time.uncens))<100){
         eval.times <- c(0,sort(time.uncens))
      } else {
         quant <- quantile(time.uncens, probs=0.95)
         eval.times <-  sort(time.uncens[time.uncens<=quant])
         if (length(eval.times)>199){
            space <- floor(length(eval.times)/100)
            eval.times <- eval.times[(1:100)*space]
         }
      }
   } else {
      eval.times <- sort(times)
   }
   eval.times <- c(0, unique(eval.times[eval.times>0]))
# # #    if(is.list(cplx)){
# # #    if (type=="apparent"){
# # #       error <- pmpec(object=model, response=response, x=x, times=eval.times,
# # #          model.args=list(complexity=cplx$stepno), type="PErr", 
# # #          external.time=full.data$time, external.status=full.data$status, ...)
# # #       }
# # # 
# # #    if (type=="noinf"){
# # #       error <- pmpec(object=model, response=response, x=x, times=eval.times,
# # #          model.args=list(complexity=cplx$stepno), type="NoInf", 
# # #          external.time=full.data$time, external.status=full.data$status, ...)
# # #       }
# # # } else {
if (type=="apparent"){
      error <- pmpec(object=model, response=response, x=x, times=eval.times,
         model.args=list(complexity=cplx), type="PErr", 
         external.time=full.data$time, external.status=full.data$status, ...)
      }

   if (type=="noinf"){
      error <- pmpec(object=model, response=response, x=x, times=eval.times,
         model.args=list(complexity=cplx), type="NoInf", 
         external.time=full.data$time, external.status=full.data$status, ...)
      }

# # # }

   attr(error, "addattr") <- eval.times
   error
}

