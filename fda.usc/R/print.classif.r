print.classif<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\n-Call:\n", deparse(x$C), "\n", sep = "")
#   cat("\n-Probability of correct classification by group (prob.classification):\n")
#   print(x$prob.classification)
if (x$C[[1]]=='classif.knn'){
    cat("\n-Optimal number of neighbors: knn.opt=",x$knn.opt,
    "\nwith highest probability of correct classification max.prob=",
    x$max.prob,"\n")}
else {if (x$C[[1]]=='classif.kernel'| x$C[[1]]=='classif.np') {
   cat("\n-Optimal bandwidth: h.opt=",x$h.opt,"with highest probability of
   correct classification: max.prob=",x$max.prob,"\n")
   }
else {if  (x$C[[1]]=='classif.gsam' | x$C[[1]]=='classif.gsam2boost'|
x$C[[1]]=='classif.glm' |x$C[[1]]== 'classif.glm2boost' |
x$C[[1]]=='classif.gkam' | x$C[[1]]=='classif.gkam2boost'|
 x$C[[1]]=='classif.tree' |   x$C[[1]]=='classif.tree2boost'){
   cat("\n-Probability of correct classification: ",round(x$max.prob,4),"\n")
     }
     if (x$C[[1]]=='classif.DD'){
    cat("\n-Probability of correct classification: ",round(1-x$misclassification,4),"\n")}
     }
}
cat("\n")
invisible(x)
}

