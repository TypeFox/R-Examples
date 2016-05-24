
print.matched.data.frame <- function(x,
                                     ...){

  object <- x
  
  cat("\n Matched by: ", object$matched.by, "\n")
  
  cat("\n Matching parameter:\n")

  print(matrix(c(round(object$match.parameter$caliper,3),
                 object$match.parameter$ratio,
                 object$match.parameter$who.treated),
               nrow=3,ncol=1,
               dimnames=
               list(c("Caliper size:",
                      "Ratio:",
                      "Who is treated?:"), c(""))))
  

  cat("\n Matching information:\n")

  print(matrix(c(object$match.parameter$givenTmatchingC,
                 object$match.parameter$bestmatch.first),
               nrow=2,ncol=1,
               dimnames=
               list(c("Untreated to treated?:",
                      "Best match?:"), c(""))))

  df.treat <- as.vector(object$data[names(object$data)==object$name.treat])
  
  cat("\n Matching data:\n")
  
  print(matrix(c(length(object$match.index[df.treat==object$match.parameters$who.treated]),
                 length(object$match.index[object$match.index>0 &
                                           df.treat==object$match.parameters$who.treated]),
                 length(object$match.index[df.treat!=object$match.parameters$who.treated]),
                 length(object$match.index[object$match.index>0 &
                                           df.treat!=object$match.parameters$who.treated]),
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
                      "Number of incomplete matching sets:"), c("")))
        )
          
}


