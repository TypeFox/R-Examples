#' @export
meanStdG = function(G, means){
  G/(means%*%t(means))
}