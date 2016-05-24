
print.matched.data.frames <- function(x,
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
  
  vec.match <-
    c(object$match.index[[1]], object$match.index[[2]])
  
  
  cat("\n Matching data:\n")

  print(matrix(c(dim(object$data[[1]])[1],
                    length(object$match.index[[1]][object$match.index[[1]]>0]),

                    dim(object$data[[2]])[1],
                    length(object$match.index[[2]][object$match.index[[2]]>0]),

                    length(object$match.index[[1]][object$match.index[[1]]>0])+
                    length(object$match.index[[2]][object$match.index[[2]]>0]),
                    
                    length(object$match.index[[1]][object$match.index[[1]]==0])+
                    length(object$match.index[[2]][object$match.index[[2]]==0]),
                    
                    length(unique(vec.match))-1,
                    
                    sum(as.numeric(table(vec.match[vec.match>0])) !=
                        (object$match.parameters$ratio +1))),
                  nrow=8,ncol=1,
                  dimnames=
                  list(c("Number of treated obs.:",
                         "Number of matched treated obs.:",
                         "Number of untreated obs.:",
                         "Number of matched untreated obs.:",
                         "Number of matched obs.:",
                         "Number of not matched obs.:",
                         "Number of matching sets:",
                         "Number of incomplete matching sets:"), c(""))))
  
}


