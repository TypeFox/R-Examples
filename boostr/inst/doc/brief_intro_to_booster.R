
## ----, cache=FALSE-------------------------------------------------------
library(mlbench)
data(Glass)
set.seed(1234)
boostedSVM1 <- 
boostr::boostWithArcX4(x = list(train = e1071::svm),
                       B = 3,
                       data = Glass,
                       .procArgs = list(
                         .trainArgs=list(
                           formula=formula(Type~.),
                           cost=100)))

boostedSVM1


## ----, cache=FALSE-------------------------------------------------------
set.seed(1234)
boostedSVM2 <-
boostr::boost(x = list(train=e1071::svm),
              B = 3,
              reweighter = boostr::arcx4Reweighter,
              aggregator = boostr::arcx4Aggregator,
              data = Glass,
              .procArgs = list(
                .trainArgs=list(
                  formula=formula(Type~.),
                  cost=100)),
              .boostBackendArgs = list(
                .reweighterArgs=list(m=0)))

boostedSVM2

identical(boostr::reweighterOutput(boostedSVM1),
          boostr::reweighterOutput(boostedSVM2))


