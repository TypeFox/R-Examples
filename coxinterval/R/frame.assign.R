### expand 'assign' attribute to indicate column index in the model frame
frame.assign <- function(frame, terms = attr(frame, "terms"),
                         matrix = model.matrix(terms, frame))
{
  labels <- attr(terms, "term.labels")
  ## account for intercept term in model frame
  assign <- c(0, match(labels[attr(matrix, "assign")], names(frame)))
  assign <- data.frame(cbind(assign, 1:ncol(matrix)))
  names(assign) <- c("frame", "matrix")
  assign
}
