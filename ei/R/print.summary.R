#prints summary function
#@x summary object


print.summary <- function(x, ...){

sum.object <- x
dec <- sum.object$precision
sum.object = sum.object[1:length(sum.object)-1]
top <- list()
num = 1
for (key in names(sum.object)){
  val <- sum.object[[key]]
  if ((is.character(val) || is.numeric(val)) && length(val) < 2){
     top[[num]] <- sum.object[[key]]
     names(top[[num]]) = key
     num=num+1
   }
}
if("Resamp" %in% names(sum.object)){
message(cat(names(top[[1]]), "=", top[[1]], ",",names(top[[2]]),
 "=", top[[2]],",",names(top[[3]]), "=", top[[3]],",",
 names(top[[4]]), "=", top[[4]],",",names(top[[5]]), "=", top[[5]]))
}

if(!("Resamp" %in% names(sum.object))){

message(cat(names(top[[1]]), "=", top[[1]], ",",names(top[[2]]),
 "=", top[[2]],",",names(top[[3]]), "=", top[[3]],",",
 names(top[[4]]), "=", top[[4]]))

}
for (key in names(sum.object)){
  val <- sum.object[[key]]
  if (!((is.character(val) || is.numeric(val)) && length(val) < 2)){
     message()
     message(key)
     print(floor(val * 10^dec)/(10^dec))
   }
}
}
