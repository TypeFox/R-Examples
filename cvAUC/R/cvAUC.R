cvAUC <- function(predictions, labels, label.ordering = NULL, folds = NULL) {
  
  # Pre-process the input
  clean <- .process_input(predictions = predictions, labels = labels, 
                          label.ordering = label.ordering, folds = NULL,
                          ids = NULL, confidence = NULL)
  
  pred <- ROCR::prediction(clean$predictions, clean$labels)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  fold_auc <- as.numeric(ROCR::performance(pred, measure = "auc", x.measure = "cutoff")@y.values)
  cvauc <- mean(fold_auc)
  
  return(list(perf = perf, fold.AUC = fold_auc, cvAUC = cvauc))
}

