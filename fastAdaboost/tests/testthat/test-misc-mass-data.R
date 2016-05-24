library(fastAdaboost)
library(MASS)

#here we have a combination of factors and numerical variables
data(Cars93)
set.seed(999)


test_that("Cars93 data works, for adaboost M1",{
  ada_obj <- adaboost(Origin~Length+Wheelbase+Width+AirBags+DriveTrain, Cars93, 10)
  pred <- predict(ada_obj, Cars93)
  print(paste("Adaboost Error on Cars93:", pred$error))
  print( table(pred$class, Cars93$Origin) )
})

test_that("Cars93 data works, for real adaboost",{  
  ada_obj <- real_adaboost(Origin~Length+Wheelbase+Width+AirBags+DriveTrain, Cars93, 10)
  pred <- predict(ada_obj, Cars93)
  #print(pred$class)
  print(paste("Real Adaboost Error on Cars93:", pred$error))
  print( table(pred$class, Cars93$Origin) )
  
})

