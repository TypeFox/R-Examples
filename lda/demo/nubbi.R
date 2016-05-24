require("ggplot2")

data(nyt)

model <- 
  nubbi.collapsed.gibbs.sampler(nyt$contexts, nyt$pair.contexts,
                                nyt$pairs,
                                K.individual = 12,
                                K.pair = 6,
                                vocab = nyt$vocab, num.iterations=25,
                                alpha = 0.1, eta = 0.1, xi = c(0.1, 0.1, 1))


## Randomly select the top 25 entities
entity.counts <- table(nyt$pairs)
top.entities <- sort(entity.counts, decreasing=TRUE)[1:25]
top.entities <- as.numeric(names(top.entities))

pairs.to.keep <- (nyt$pairs[,1] %in% top.entities) &
                 (nyt$pairs[,2] %in% top.entities)

strength <- apply(model$document_source_sums[pairs.to.keep,],
                  1, function(x) x[3] / sum(x))


pair.topic <- sapply(model$document_sums[(length(nyt$contexts) + 1):
                                         (length(model$document_sums))],
                     which.max)

qplot(x=nyt$names[nyt$pairs[pairs.to.keep,1] + 1],
      y=nyt$names[nyt$pairs[pairs.to.keep,2] + 1],
      size=strength,
      colour=as.factor(pair.topic[pairs.to.keep]),
      xlab = "",
      ylab = "") +
  scale_color_discrete("Relationship topic") +
  scale_size("Relationship strength") +
  theme(axis.text.x = element_text(angle=90, hjust=1))

