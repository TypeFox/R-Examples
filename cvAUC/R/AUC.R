AUC <- function(predictions, labels, label.ordering = NULL) { 
  
  # Pre-process the input
  clean <- .process_input(predictions = predictions, labels = labels, 
                          label.ordering = label.ordering, folds = NULL,
                          ids = NULL, confidence = NULL)
  
  pred <- ROCR::prediction(clean$predictions, clean$labels)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  auc <- as.numeric(ROCR::performance(pred, measure = "auc", x.measure = "cutoff")@y.values)
  
  return(auc)
}
