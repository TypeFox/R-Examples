# AUTHOR: Timothy P. Jurka
# DESCRIPTION: Demo presented at the 2011 CAP Catania Conference.

library(RTextTools)

data(NYTimes)

matrix <- create_matrix(cbind(NYTimes["Title"],NYTimes["Subject"]), language="english")

container <- create_container(matrix,NYTimes$Topic.Code,trainSize=1:2500, testSize=2501:3100, virgin=FALSE)


# TRAINING
svm_model <- train_model(container,"SVM")
maxent_model <- train_model(container,"MAXENT")


# PREDICTION
svm_results <- classify_model(container,svm_model)
maxent_results <- classify_model(container,maxent_model)


# ANALYTICS
analytics <- create_analytics(container,cbind(svm_results,maxent_results))

analytics@label_summary
analytics@algorithm_summary
analytics@document_summary
analytics@ensemble_summary


# WRITE RESULTS TO CSV
write.csv(analytics@label_summary,"SampleData_LabelSummary.csv")
write.csv(analytics@algorithm_summary,"SampleData_AlgorithmSummary.csv")
write.csv(analytics@document_summary,"SampleData_DocumentSummary.csv")
write.csv(analytics@ensemble_summary,"SampleData_EnsembleSummary.csv")