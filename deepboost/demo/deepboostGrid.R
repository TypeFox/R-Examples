library(deepboost)

data("sonar")
formula <- R ~ .
best_params <-
  deepboost.gridSearch(formula, sonar)

boost <- deepboost(formula, sonar,
                   num_iter = best_params[2][[1]],
                   beta = best_params[3][[1]],
                   lambda = best_params[4][[1]],
                   loss_type = best_params[5][[1]]
)

print(boost)

preds <- predict(boost, sonar)
