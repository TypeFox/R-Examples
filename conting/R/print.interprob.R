print.interprob <-
function(x,digits = max(3, getOption("digits") - 3),...){

nam<-x$term
df<-data.frame(post_prob=x$prob)
row.names(df)<-nam

cat("Posterior probabilities of log-linear parameters:\n")
print(df,digits=digits)}
