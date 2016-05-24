post <-
function(object, level="upper"){
 if (level=="upper" && !is.null(object$Misc$mform2)){
      data <- object$data
      mform2 <- object$Misc$mform2
      group <- object$data[mform2][row.names(unique(data[mform2])),]
      post.prob <- object$post.prob[row.names(unique(data[mform2])),] ; dimnames(post.prob)[[1]]<- group
      post.int  <- object$post.int[row.names(unique(data[mform2]))]; names(post.int)<- group
  } else {
      post.prob<- object$post.prob
      post.int <- object$post.int
  }
  classif <- apply(post.prob, 1, which.max)
  names(classif) <-  dimnames(post.prob)[[1]]
  
  if (is.null(object$Misc$mform2)){ level<- "none"}
  
  return(list(prob= post.prob, int=post.int, classif = classif, level=level))
}