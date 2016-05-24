### R code from vignette source 'consensus_auc.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: consensus_auc.Rnw:21-22
###################################################
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 2: consensus_auc.Rnw:25-33
###################################################
library(MASS)
data(fgl)
set.seed(100)
ind.train <- sample(nrow(fgl)
                   , size = floor(nrow(fgl)*.7)
                   , replace = FALSE
                   )
ind.eval <- setdiff(seq(1:nrow(fgl)), ind.train)


###################################################
### code chunk number 3: consensus_auc.Rnw:36-49
###################################################
library(rpart)
set.seed(123)
fgl.rpart <- rpart(type ~ .
                 , data = fgl[ind.train, ]
                 , parms = list(split = "information")
                 )

newcp <- max(fgl.rpart$cptable[,"CP"] *
	     as.numeric(fgl.rpart$cptable[,"xerror"] <
			sum(fgl.rpart$cptable[dim(fgl.rpart$cptable)[1]
                                              , c("xerror","xstd")]))
	     ) + 1e-13
fgl.rpart.pruned <- prune(fgl.rpart, cp = newcp)


###################################################
### code chunk number 4: consensus_auc.Rnw:52-57
###################################################
library(nnet)
fgl.multinom <- multinom(type ~ .
                            , data = fgl[ind.train, ]
			    , trace = FALSE
                            )


###################################################
### code chunk number 5: consensus_auc.Rnw:60-69
###################################################
library(mda)
confusion(predict(fgl.rpart.pruned
		    , newdata = fgl[ind.eval, ]
		    , type = "class")
          ,fgl[ind.eval, ]$type)
confusion(predict(fgl.multinom
		    , newdata = fgl[ind.eval, ]
		    , type = "class") 
	  , fgl[ind.eval, ]$type)


###################################################
### code chunk number 6: consensus_auc.Rnw:72-84
###################################################
library(HandTill2001)
auc(multcap(response = fgl[ind.eval, ]$type
	     , predicted = predict(fgl.rpart.pruned
				   , newdata = fgl[ind.eval, ])
	     ))

auc(multcap(response = fgl[ind.eval, ]$type
	     , predicted = predict(fgl.multinom
				   , newdata = fgl[ind.eval, ]
				   , type = "probs"
				   )
	     ))


###################################################
### code chunk number 7: consensus_auc.Rnw:89-95
###################################################
set.seed(100)
ind.inner.train <- sample(ind.train
                         , size = floor(length(ind.train)*0.7)
                         , replace = FALSE
                         )
ind.inner.eval <- setdiff(ind.train, ind.inner.train)


###################################################
### code chunk number 8: consensus_auc.Rnw:98-112
###################################################
wa.fgl.multinom <- multinom(fgl.multinom
                            , data = fgl[ind.inner.train, ]
			    , trace = FALSE
                            )
wa.fgl.rpart <- rpart(type ~ .
                     , data = fgl[ind.inner.train, ]
                     , parms = list(split = "information")
                     )
newcp <- max(wa.fgl.rpart$cptable[,"CP"] *
             as.numeric(wa.fgl.rpart$cptable[,"xerror"] <
                        sum(wa.fgl.rpart$cptable[dim(wa.fgl.rpart$cptable)[1]
                                                , c("xerror","xstd")]))
             ) + 1e-13
wa.fgl.rpart.pruned <- prune(wa.fgl.rpart, cp = newcp)


###################################################
### code chunk number 9: consensus_auc.Rnw:115-123
###################################################
li <- list()
li$rpart$auc <- auc(multcap(response = fgl[ind.inner.eval, ]$type
                             , predicted = predict(wa.fgl.rpart.pruned
                                 , newdata = fgl[ind.inner.eval, ])))
li$mllm$auc <- auc(multcap(response = fgl[ind.inner.eval, ]$type
                            , predicted = predict(wa.fgl.multinom
                                , newdata = fgl[ind.inner.eval, ]
                                , type = "probs")))


###################################################
### code chunk number 10: consensus_auc.Rnw:126-131
###################################################
li$rpart$predictions <- predict(fgl.rpart.pruned
                                , newdata = fgl[ind.eval, ])
li$mllm$predictions <- predict(fgl.multinom
                               , newdata = fgl[ind.eval, ]
                               , type = "probs")


###################################################
### code chunk number 11: consensus_auc.Rnw:134-137
###################################################
predicted <- Reduce('+', lapply(li, function(x)
				x$auc * x$predictions)
) / Reduce('+', sapply(li,"[", "auc"))


###################################################
### code chunk number 12: consensus_auc.Rnw:140-143
###################################################
auc(multcap(response = fgl[ind.eval, ]$type
             , predicted = predicted)
     )


###################################################
### code chunk number 13: consensus_auc.Rnw:146-154
###################################################
classes.predicted <-
  factor(x =
         apply(predicted
               , 1
               , function(x)
               dimnames(predicted)[[2]][which.max(x)])
	 , levels = levels(fgl[ind.eval, ]$type)	 
         )


