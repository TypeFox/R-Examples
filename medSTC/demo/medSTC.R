library(medSTC)
set.seed(1)

## Use the newsgroup data set.
data(newsgroup.train.documents)
data(newsgroup.test.documents)
data(newsgroup.train.labels)
data(newsgroup.test.labels)

num.topics <- 50
num.classes<-3

cat("Note that the demo takes several minutes.\n")
model<-medSTC(documents=newsgroup.train.documents, 
				mlabels=newsgroup.train.labels, 
				ntopics=num.topics, 
				initial_c=0.5, 
				lambda=0.1, 
				rho=0.01,  
				delta_ell=360, 
				supervised=TRUE, 
				primal_svm=1, 
				var_max_iter=25, 
				convergence=1e-4, 
				em_max_iter=40, 
				em_convergence=1e-4, 
				svm_alg_type=2) 


predictions <- predict(model,newsgroup.test.documents)
predInd<-predictions$assignments
cat("Individual classification accuracy with the first", num.classes,"classes of newsgroup data set with 25 e, 40 m steps,",num.topics, "topics\n", sum(predInd==newsgroup.test.labels)/length(newsgroup.test.labels),"\n")
cat("Confusion matrix:\n")
print(table(predInd,newsgroup.test.labels))

getTopWordsForTopics<-
function(model, vocab, topk=5){
	data(stopwords)
	 z<-match(stopwords,newsgroup.vocab,nomatch=0)
	 z<-z[z>0]
	x<-apply(model$state$LogProbabilityOfWordsForTopics,2,function(x) { y<-order(x,decreasing=TRUE); y<-y[!y%in%z]; vocab[y[1:topk]]   })
	colnames(x)<-paste("topic",1:ncol(model$state$LogProbabilityOfWordsForTopics))	
	x
}
data(newsgroup.vocab)
getTopWordsForTopics(model, newsgroup.vocab)