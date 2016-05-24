`subscales` <-
function(items, scales, scale.names=NA, score.items=FALSE, check.reliability=FALSE, key=NA){
n.scales <- ncol(scales)

if(score.items) {
  save.names <- colnames(items)
  items <- as.data.frame(score(items,key,output.scored=TRUE)$scored)
  colnames(items) <- save.names
}

sets <- apply(scales,2,function(XXX) items[,XXX==1])
suppressWarnings(if(! is.na(scale.names)) names(sets) <- scale.names  else names(sets)<-paste("Q.",c(seq(1:n.scales)),sep=""))

if(check.reliability){
  for(i in 1:n.scales){
   sets[[i]] <- suppressWarnings(score(sets[[i]], output.scored=TRUE,rel=TRUE))
  }}
  else {  for(i in 1:n.scales){ sets[[i]] <- suppressWarnings(score(sets[[i]], output.scored=TRUE)) }}
sets
}

