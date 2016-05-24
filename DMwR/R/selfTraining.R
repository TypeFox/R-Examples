
# =====================================================
# Function that can be used to to apply self training on
# a given classifier.
# =====================================================
# Luis Torgo, Feb 2010
# =====================================================
# Example run:
# library(DMwR)
# library(e1071)
# data(iris)
# idx <- sample(150,100)
# tr <- iris[idx,]
# ts <- iris[-idx,]
# nb <- naiveBayes(Species ~ .,tr)
# table(predict(nb,ts),ts$Species)
# trSS <- tr
# nas <- sample(100,50)
# trSS[nas,'Species'] <- NA
# func <- function(m,d) {
#    p <- predict(m,d,type='raw')
#    data.frame(cl=colnames(p)[apply(p,1,which.max)],p=apply(p,1,max))
# }
# nbSSbase <- naiveBayes(Species ~ .,trSS[-nas,])
# table(predict(nbSSbase,ts),ts$Species)
# nbSS <- SelfTrain(Species ~ .,trSS,learner('naiveBayes',list()),'func')
# table(predict(nbSS,ts),ts$Species)
#
SelfTrain <- function(form,data,
                      learner,
                      predFunc,
                      thrConf=0.9,
                      maxIts=10,percFull=1,
                      verbose=F)
  {
    N <- NROW(data)
    it <- 0
    sup <- which(!is.na(data[,as.character(form[[2]])]))
    repeat {
      it <- it+1
      model <- runLearner(learner,form,data[sup,])
      probPreds <- do.call(predFunc,list(model,data[-sup,]))
      new <- which(probPreds[,2] > thrConf)
      if (verbose) cat('IT.',it,'\t nr. added exs. =',length(new),'\n')
      if (length(new)) {
        data[(1:N)[-sup][new],as.character(form[[2]])] <- probPreds[new,1]
        sup <- c(sup,(1:N)[-sup][new])
      } else break
      if (it == maxIts || length(sup)/N >= percFull) break
    }
    return(model)
  }
