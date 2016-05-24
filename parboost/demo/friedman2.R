# # Demo not run
# data(friedman2)
# cl <- makeCluster(2)
# parboost_model <- parboost(cluster_object = cl, data = friedman2,
#                            nsplits = 4, seed = 17, formula = y ~ .,
#                            baselearner="bbs", postprocessing = "glm",
#                            control = boost_control(mstop=500))
# stopCluster(cl)
# print(parboost_model)
# summary(parboost_model)
# head(predict(parboost_model))
