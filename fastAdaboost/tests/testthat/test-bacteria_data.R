library(fastAdaboost)
library(MASS)
  
test_that("Bacteria dataset works",{
    data(bacteria)
    boost_obj <- adaboost(y~.,bacteria , 10)
    pred <- predict(boost_obj,bacteria)
    print(paste("Adaboost Error:",pred$error))
    print( table(pred$class,bacteria$y) )
})
#  
test_that("bacteria dataset works with a selection of variables",{
    data(bacteria)
    boost_obj <- adaboost(y~ap+hilo+week,bacteria , 10)
    pred <- predict(boost_obj, bacteria)
    print(paste("Adaboost Error:",pred$error))
    print( table(pred$class,bacteria$y) )
})

test_that("bacteria dataset works with unlabeled data",{
    data(bacteria)
    boost_obj <- adaboost(y~ap+hilo+week,bacteria , 10)
    bacteria_2 <-bacteria
    unlabeled_cols <-names(bacteria_2)[!(names(bacteria_2) %in% "y") ]
    bacteria_2 <- bacteria_2[,unlabeled_cols]
    pred <- predict(boost_obj, bacteria_2)
    print(paste("Adaboost Error:",pred$error))
    expect_true(is.na(pred$error))
})
  
