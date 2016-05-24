
print.matched.pscore <- function(x,
                                 ...){

  object <- x
  
  cat("\n Matched by: ", object$matched.by, "\n")
  
  cat("\n Matching parameter:\n")

  print(matrix(c(round(object$match.parameters$caliper,3),
                 object$match.parameters$ratio,
                 object$match.parameters$who.treated),
               nrow=3,ncol=1,
               dimnames=
               list(c("Caliper size:",
                      "Ratio:",
                      "Who is treated?:"), c(""))))
  

  cat("\n Matching information:\n")

  print(matrix(c(object$match.parameters$givenTmatchingC,
                 object$match.parameters$bestmatch.first),
               nrow=2,ncol=1,
               dimnames=
               list(c("Untreated to treated?:",
                      "Best match?:"), c(""))))
  

  cat("\n Matching data:\n")

  print(matrix(c(dim(object$data[object$treat==object$match.parameters$who.treated,])[1],
                 dim(object$data[object$treat==object$match.parameters$who.treated &
                                 object$match.index>0,])[1],
                 dim(object$data[object$treat!=object$match.parameters$who.treated,])[1],
                 dim(object$data[object$treat!=object$match.parameters$who.treated &
                                 object$match.index>0,])[1],
                 length(object$match.index[object$match.index>0]),
                 length(object$match.index[object$match.index==0]),
                 length(unique(object$match.index))-1,
                 sum(as.numeric(table(object$match.index[object$match.index>0])) !=
                     (object$match.parameters$ratio +1))),
               nrow=8,ncol=1,
               dimnames=
               list(c("Number of treated obs.:",
                      "Number of matched treated obs.:",
                      "Number of untreated obs.:",
                      "Number of matched untreated obs.:",
                      "Number of total matched obs.:",
                      "Number of not matched obs.:",
                      "Number of matching sets:",
                      "Number of incomplete matching sets:"), c(""))))

}


