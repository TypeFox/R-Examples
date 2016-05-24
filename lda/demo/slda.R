set.seed(8675309)

## Use the political blogs data set.
data(poliblog.documents)
data(poliblog.vocab)
data(poliblog.ratings)

num.topics <- 10

## Initialize the params
params <- sample(c(-1, 1), num.topics, replace=TRUE)

result <- slda.em(documents=poliblog.documents,
                  K=num.topics,
                  vocab=poliblog.vocab,
                  num.e.iterations=10,
                  num.m.iterations=4,
                  alpha=1.0, eta=0.1,
                  poliblog.ratings / 100,
                  params,
                  variance=0.25,
                  lambda=1.0,
                  logistic=FALSE,
                  method="sLDA")

## Make a pretty picture.
require("ggplot2")
Topics <- apply(top.topic.words(result$topics, 5, by.score=TRUE),
                2, paste, collapse=" ")
coefs <- data.frame(coef(summary(result$model)))
theme_set(theme_bw())
coefs <- cbind(coefs, Topics=factor(Topics, Topics[order(coefs$Estimate)]))
coefs <- coefs[order(coefs$Estimate),]
qplot(Topics, Estimate, colour=Estimate, size=abs(t.value), data=coefs) +
  geom_errorbar(width=0.5, aes(ymin=Estimate-Std..Error,
                  ymax=Estimate+Std..Error)) + coord_flip()

predictions <- slda.predict(poliblog.documents,
                            result$topics, 
                            result$model,
                            alpha = 1.0,
                            eta=0.1)

qplot(predictions,
      fill=factor(poliblog.ratings),
      xlab = "predicted rating",
      ylab = "density",
      alpha=I(0.5),
      geom="density") +
  geom_vline(aes(xintercept=0)) +
  theme(legend.position = "none")

predicted.docsums <- slda.predict.docsums(poliblog.documents,
                                          result$topics, 
                                          alpha = 1.0,
                                          eta=0.1)
predicted.proportions <- t(predicted.docsums) / colSums(predicted.docsums)

qplot(`Topic 1`, `Topic 2`, 
      data = structure(data.frame(predicted.proportions), 
                       names = paste("Topic", 1:10)), 
      size = `Topic 3`)
