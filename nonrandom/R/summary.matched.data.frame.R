
summary.matched.data.frame <- function(object,
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

  df.treat <- as.vector(object$data[names(object$data)==object$name.treat])

  mat3 <- matrix(c(length(object$match.index[df.treat==object$match.parameters$who.treated]),
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

  mat <- list(match.var = object$matched.by,
              parameter = mat1,
              info      = mat2,
              number    = mat3)

  class(mat) <- "summary.matched.data.frame"

  mat
  
}


