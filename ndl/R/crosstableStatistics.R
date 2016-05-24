crosstableStatistics <- function(ctable)
{ N <- sum(ctable);
  # observed as row margins, predicted as column margins
  # according to Menard (1995: 24-32)
  sum.row <- apply(ctable,1,sum);
  sum.col <- apply(ctable,2,sum);
  correct.with.model <- sum(diag(ctable));
  errors.with.model <- N - correct.with.model;
  errors.without.model.prediction <- N - max(sum.row);
  errors.without.model.classification <- sum(sum.row*((N-sum.row)/N));
  lambda.p <- 1-(errors.with.model/errors.without.model.prediction);
  d.lambda.p <- (errors.without.model.prediction/N-errors.with.model/N)/sqrt((errors.without.model.prediction/N)*(1-errors.without.model.prediction/N)/N);
  p.lambda.p <- 1-pnorm(d.lambda.p);
  tau.p <-  1-(errors.with.model/errors.without.model.classification);
  d.tau.p <- (errors.without.model.classification/N-errors.with.model/N)/sqrt((errors.without.model.classification/N)*(1-errors.without.model.classification/N)/N);
  p.tau.p <- 1-pnorm(d.tau.p);
  accuracy <- sum(diag(ctable))/N;
  recall.predicted <- diag(ctable)/sum.row;
  precision.predicted <- diag(ctable)/sum.col;
  statistics <- list(accuracy = accuracy, recall.predicted=recall.predicted, precision.predicted=precision.predicted, lambda.prediction = lambda.p, tau.classification = tau.p, d.lambda.prediction = d.lambda.p, d.tau.classification = d.tau.p, p.lambda.prediction = p.lambda.p, p.tau.classification = p.tau.p);
  return(statistics);
} 
