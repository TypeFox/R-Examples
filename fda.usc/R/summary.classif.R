summary.classif<-function (object, ...)
{
 cat("     - SUMMARY - \n")
 cat("\n-Probability of correct classification by group (prob.classification):\n")
 print(object$prob.classification)
 cat("\n-Confusion matrix between the theoretical groups (by rows)
  and estimated groups (by column) \n")
 print(table(object$group,object$group.est))
 if (object$C[[1]]=='classif.np'){
    if (object$ty=="S.NW") {
   cat("\n-Vector of probability of correct classification
    by banwidth (h):\n")
    print(round(1-object$gcv,4))
#   cat("\n-Functional measure of closeness (optimal distance, h.opt):\n")
#   print(round(object$h.opt,4))

cat("\n-Optimal bandwidth: h.opt=",object$h.opt,"with highest probability of
correct classification: max.prob=",object$max.prob,"\n")
  }  }
 if (object$C[[1]]=='classif.np'){
   if (object$ty=="S.KNN") {
   cat("\n-Vector of probability of correct classification
   by number of neighbors (knn):\n")
    print(round(1-object$gcv,4))
    cat("\n-Optimal number of neighbors: knn.opt=",object$h.opt,
    "\nwith highest probability of correct classification max.prob= ",
    object$max.prob,"\n")
    }}
    if (object$C[[1]]=='classif.gkam' | object$C[[1]]=='classif.gkam2boost' |
     object$C[[1]]=='classif.gsam2boost' | object$C[[1]]=='classif.gsam'
     | object$C[[1]]=='classif.glm'| object$C[[1]]=='classif.glm2boost'){
   cat("\n-Probability of correct classification: ",round(object$max.prob,4),"\n")
    }
    if (object$C[[1]]=='classif.tree'|object$C[[1]]=='classif.tree2boost'){
   cat("\n-Probability of correct classification: ",round(object$max.prob,4),"\n")}
   if (object$C[[1]]=='classif.DD'){
   cat("\n-Probability of correct classification: ",round(1-object$misclassification,4),"\n")}
cat("\n")
output<-object
}

