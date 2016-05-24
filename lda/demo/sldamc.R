library(lda)
set.seed(8675309)

## Use the newsgroup data set.
data(newsgroup.train.documents)
data(newsgroup.test.documents)
data(newsgroup.train.labels)
data(newsgroup.test.labels)
data(newsgroup.vocab)
data(newsgroup.label.map)

chosenClasses<-1:3
classes<-which(newsgroup.train.labels %in% chosenClasses)

num.topics <- 30
num.classes<-length(chosenClasses)

## Initialize the params
params <- sample(c(1,-1),num.topics*(num.classes-1),replace=TRUE)
integerAnnotations<-as.integer(factor(newsgroup.train.labels[classes]))-1L #needs to be zero based and integer

result <- slda.em(documents=newsgroup.train.documents[classes],
                  K=num.topics,
                  vocab=newsgroup.vocab,
                  num.e.iterations=10,
                  num.m.iterations=4,
                  alpha=1.0, eta=0.1,
                  integerAnnotations,
                  params,
                  variance=0.25,
                  lambda=1.0,
                  logistic=TRUE,
                  method="sLDA")

classesInTest<-which(newsgroup.test.labels %in% chosenClasses)
predictions <- slda.predict(newsgroup.test.documents[classesInTest],
                            result$topics, 
                            result$model,
                            alpha = 1.0,
                            eta=0.1)

predProbs<-cbind(0,predictions)
predInd<-t(apply(predProbs,1,which.max))
cat("Individual classification accuracy with the first", num.classes,"classes of newsgroup data set with 10 e, 4 m steps,",num.topics, "topics\n", sum(predInd==newsgroup.test.labels[classesInTest])/length(classesInTest),"\n")
cat("Confusion matrix:\n")
print(table(predInd,newsgroup.test.labels[classesInTest]))

