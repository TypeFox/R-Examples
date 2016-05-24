`compare` <-
function(obj1, obj2){
  if(!("msc" %in% class(obj1)) | !("msc" %in% class(obj2))) stop("both arguments must be msc objects")
  if(length(obj1$msc.List) != length(obj2$msc.List) | obj1$stepSize != obj2$stepSize) stop("objects not comparable")
  obj1$Sum.Stats = obj1$Sum.Stats - obj2$Sum.Stats
  obj1$Sum.Stats[,1:length(obj1$msc.List)] = obj2$Sum.Stats[,1:length(obj2$msc.List)]
  print(summary(obj1$Sum.Stats[,-(1:length(obj1$msc.List))]))
  plot(obj1)
}

