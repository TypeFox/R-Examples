context("Main functions: hybridEnsemble, predict and CVhybridEnsemble")



test_that("Test output hybridEnsemble and predict", {
skip_on_cran()
data(Credit)
hE <-hybridEnsemble(x=Credit[1:100,names(Credit) != 'Response'],
                    y=Credit$Response[1:100],
                    verbose=FALSE,
                    combine=c("rbga","DEopt","GenSA","malschains","psoptim","soma","tabu","NNloglik","GINNLS","LHNNLS"),
                    RF.ntree=50,
                    AB.iter=50,
                    NN.size=5,
                    NN.decay=0,
                    SV.gamma = 2^-15,
                    SV.cost = 2^-5,
                    SV.degree=2,
                    SV.kernel='radial',
                    tabu.iters=20,
                    tabu.listSize=c(5))
  
predictions <- predict(hE, newdata=Credit[1:100,names(Credit) != 'Response'])

expect_equal(class(hE),"hybridEnsemble")
expect_output(str(hE),"List of 40")
expect_output(str(predictions),"List of 14")

})



test_that("Test output CVhybridEnsemble", {
skip_on_cran()
data(Credit)

CVhE <- CVhybridEnsemble(x=Credit[1:200,names(Credit) != 'Response'],
                    y=Credit$Response[1:200],
                    verbose=FALSE,
                    filter=0.05,
                    KF.rp=1,
                    RF.ntree=50,
                    AB.iter=50,
                    NN.size=5,
                    NN.decay=0,
                    SV.gamma = 2^-15,
                    SV.cost = 2^-5,
                    SV.degree=2,
                    SV.kernel='radial')
expect_output(str(CVhE),"List of 4")
})


