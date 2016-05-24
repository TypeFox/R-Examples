`score.transform` <-
function(scores, mu.new=0, sd.new=1, normalize=FALSE){

percentile <- trunc(rank(scores))/length(scores)

if(normalize) scores.new <- qnorm(percentile) 
else{
  mu.old <- mean(scores)
  sd.old <- sd(scores)
  scores.new <- (scores - mu.old)/sd.old
}

scores.new <- (scores.new*sd.new) + mu.new
out <- list(new.scores = scores.new, p.scores = percentile)
out
}

