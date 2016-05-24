
summary.matched.data.frames <- function(object,
                                        ...){

  mat1 <- matrix(c(round(object$match.parameters$caliper,3),
                   object$match.parameters$ratio,
                   object$match.parameters$who.treated),
                 nrow=3,ncol=1,
                 dimnames=
                 list(c("Caliper size:",
                        "Ratio:",
                        "Who is treated?:"), c("")))

  
  mat2 <- matrix(c(object$match.parameters$givenTmatchingC,
                   object$match.parameters$bestmatch.first),
                 nrow=2,ncol=1,
                 dimnames=
                 list(c("Untreated to treated?:",
                        "Best match?:"), c("")))

  
  vec.match <-
    c(object$match.index[[1]], object$match.index[[2]])

   mat3 <- matrix(c(dim(object$data[[1]])[1],
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
                         "Number of incomplete matching sets:"), c("")))

  mat <- list(match.var = object$matched.by,
              parameter = mat1,
              info      = mat2,
              number    = mat3)
  
  class(mat) <- "summary.matched.data.frames"
  
  mat
  
}


