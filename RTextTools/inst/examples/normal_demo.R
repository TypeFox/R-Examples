# AUTHOR: Timothy P. Jurka
# DESCRIPTION: This file demonstrates the standard operation of RTextTools.

# LOAD THE RTextTools LIBARY
library(RTextTools)


# READ THE CSV DATA
data(NYTimes)


# [OPTIONAL] SUBSET YOUR DATA TO GET A RANDOM SAMPLE
NYTimes <- NYTimes[sample(1:3100,size=3100,replace=FALSE),]


# CREATE A TERM-DOCUMENT MATRIX THAT REPRESENTS WORD FREQUENCIES IN EACH DOCUMENT
# WE WILL TRAIN ON THE Title and Subject COLUMNS
matrix <- create_matrix(cbind(NYTimes["Title"],NYTimes["Subject"]), language="english", removeNumbers=TRUE, stemWords=TRUE, weighting=weightTfIdf)


# CREATE A container THAT IS SPLIT INTO A TRAINING SET AND A TESTING SET
# WE WILL BE USING Topic.Code AS THE CODE COLUMN. WE DEFINE A 2000 
# ARTICLE TRAINING SET AND A 1000 ARTICLE TESTING SET.
container <- create_container(matrix,NYTimes$Topic.Code,trainSize=1:3000, testSize=3001:3100,virgin=FALSE)


# THERE ARE TWO METHODS OF TRAINING AND CLASSIFYING DATA.
# ONE WAY IS TO DO THEM AS A BATCH (SEVERAL ALGORITHMS AT ONCE)
models <- train_models(container, algorithms=c("GLMNET","MAXENT","SVM"))
results <- classify_models(container, models)


# ANOTHER WAY IS TO DO THEM ONE BY ONE.
glmnet_model <- train_model(container,"GLMNET")
maxent_model <- train_model(container,"MAXENT")
svm_model <- train_model(container,"SVM")

glmnet_results <- classify_model(container,glmnet_model)
maxent_results <- classify_model(container,maxent_model)
svm_results <- classify_model(container,svm_model)

# USE print_algorithms() TO SEE ALL AVAILABLE ALGORITHMS.
print_algorithms()


# VIEW THE RESULTS BY CREATING ANALYTICS
# IF YOU USED OPTION 1, YOU CAN GENERATE ANALYTICS USING
analytics <- create_analytics(container, results)

# IF YOU USED OPTION 2, YOU CAN GENERATE ANALYTICS USING:
analytics <- create_analytics(container,cbind(svm_results,maxent_results))

# RESULTS WILL BE REPORTED BACK IN THE analytics VARIABLE.
# analytics@algorithm_summary: SUMMARY OF PRECISION, RECALL, F-SCORES, AND ACCURACY SORTED BY TOPIC CODE FOR EACH ALGORITHM
# analytics@label_summary: SUMMARY OF LABEL (e.g. TOPIC) ACCURACY
# analytics@document_summary: RAW SUMMARY OF ALL DATA AND SCORING
# analytics@ensemble_summary: SUMMARY OF ENSEMBLE PRECISION/COVERAGE. USES THE n VARIABLE PASSED INTO create_analytics()

head(analytics@algorithm_summary)
head(analytics@label_summary)
head(analytics@document_summary)
head(analytics@ensemble_summary)

# WRITE OUT THE DATA TO A CSV
write.csv(analytics@algorithm_summary,"SampleData_AlgorithmSummary.csv")
write.csv(analytics@label_summary,"SampleData_LabelSummary.csv")
write.csv(analytics@document_summary,"SampleData_DocumentSummary.csv")
write.csv(analytics@ensemble_summary,"SampleData_EnsembleSummary.csv")