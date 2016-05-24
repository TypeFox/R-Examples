setMethod("group", signature(object="flexmix"), function(object) {
  group <- object@group
  if (!length(group)) group <- group(object@model[[1]])
  group
})

setMethod("group", signature(object="FLXM"), function(object) {
  factor(seq_len(nrow(object@x)))
})

setMethod("group", signature(object="FLXMRglmfix"), function(object) {
  factor(seq_len(nrow(object@x)/sum(object@nestedformula@k)))
})
