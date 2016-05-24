
library(caret)
library(ada)
library(deepboost)

# read datasets
data("adult")
data("australian")
data("banana")
data("bupa")
data("coli2000")
data("haberman")
data("heart")
data("magic")
data("pima")
data("sonar")

# create lists of datasets and formulas
datasets <- list(adult=adult, aust=australian, banana=banana, bupa=bupa, coli=coli2000,
                 haber=haberman, heart=heart, magic=magic, pima=pima, sonar=sonar)
formulas <- list(X..50K ~ X39 + X77516 + X13 + X2174 +  X0 + X40,
                 X0.3 ~ .,
                 X.1.0 ~ .,
                 X1 ~ .,
                 X0.45 ~ .,
                 negative ~ .,
                 X2.2 ~ .,
                 g ~ .,
                 tested_positive ~ .,
                 R ~ .)

results <- data.frame(dataset = numeric(0), ensemble_size = numeric(0), ada_acc = numeric(0), ada_sd = numeric(0),
                      ada_time = numeric(0), deep_acc = numeric(0), deep_sd = numeric(0), deep_time = numeric(0),
                      t_test = numeric(0))
# for each number of iterations
for(num_iter in c(5,10,20,50)){
  # for each data set
  for(i in c(2,4,6,7,9,10)){
    ds <- datasets[[i]]
    levels(ds[,length(ds)]) <- c(1,-1)
    formula <- formulas[[i]]
    ada_acc <- rep(0,5)
    deep_acc <- rep(0,5)
    ada_t <- 0
    deep_t <- 0
    # 5 different 10folds
    for(j in 1:5){
      flds <- createFolds(1:nrow(ds), k = 10)
      for(k in 1:10){
        l <- (k%%10)+1
        eval_train <- ds[-flds[[l]],]
        eval_test <- ds[flds[[l]],]
        train <- ds[-flds[[k]],]
        test <- ds[flds[[k]],]

        beta_vals = c(2^-0, 2^-1, 2^-2, 2^-3, 2^-4, 2^-5, 2^-6)
        lambda_vals = c(0.0001, 0.005, 0.01, 0.05, 0.1, 0.5)
        dpbGrid <-  expand.grid(beta = beta_vals,
                                lambda = lambda_vals)

        # train ADABOOST
        best_acc = 0
        best_nu = 0
        for(nu in beta_vals){
          eval_model <- ada(formula, eval_train, iter = num_iter, nu=nu)
          acc <-  sum(predict(eval_model, eval_test) == eval_test[,length(eval_test)]) / nrow(eval_test)
          if(acc > best_acc){
            best_acc <- acc
            best_nu <- nu
          }
        }

        t <- Sys.time()
        ab_model <- ada(formula, train, iter = num_iter, nu=best_nu)
        ada_acc[j] <- ada_acc[j] + sum(predict(ab_model, test) == test[,length(test)]) / nrow(test)
        ada_t <- ada_t + round(difftime(Sys.time(), t, units = "secs"), 2)


        # train DEEPBOOST
        best_acc = 0
        best_lambda = 0
        best_beta = 0
        for(grow in 1:nrow(dpbGrid)){
          beta <- dpbGrid[grow,"beta"]
          lambda <- dpbGrid[grow,"lambda"]
          eval_model <- deepboost.formula(formula, eval_train, num_iter = num_iter, beta = beta, lambda = lambda, verbose = F)
          acc <-  sum(predict(eval_model, eval_test) == eval_test[,length(eval_test)]) / nrow(eval_test)
          if(acc > best_acc){
            best_acc <- acc
            best_lambda <- lambda
            best_beta <- beta
          }
        }

        t <- Sys.time()
        db_model <- deepboost.formula(formula, train, num_iter = num_iter, beta = best_beta, lambda = best_lambda, verbose = F)
        deep_acc[j] <- deep_acc[j] + sum(predict(db_model, test) == test[,length(test)]) / nrow(test)
        deep_t <- deep_t + round(difftime(Sys.time(), t, units = "secs"), 2)
      }
      ada_acc[j] <- ada_acc[j]/10.0
      deep_acc[j] <- deep_acc[j]/10.0
    }
    # caluculate results
    ada_acc_mean <- round(mean(ada_acc), 4)
    #ada_auc_mean <- mean(ada_auc)
    deep_acc_mean <- round(mean(deep_acc), 4)
    #deep_auc_mean <- mean(deep_auc)
    ada_acc_sd <- round(sd(ada_acc), 6)
    #ada_auc_sd <- sd(ada_auc)
    deep_acc_sd <- round(sd(deep_acc), 6)
    #deep_auc_sd <- sd(deep_auc)
    acc_t_test <- t.test(ada_acc, deep_acc, paired=TRUE)$p.value < 0.05
    #auc_t_test <- t.test(ada_auc, deep_auc, paired=TRUE)$p.value < 0.05

    # print to file
    fname <- paste('./', names(datasets)[i], num_iter, ".res", sep='')
    res <- data.frame(dataset = names(datasets)[i], ensemble_size = num_iter, ada_acc = ada_acc_mean,
                      ada_sd = ada_acc_sd, ada_time = ada_t, deep_acc = deep_acc_mean,
                      deep_sd = deep_acc_sd, deep_time = deep_t,  t_test = acc_t_test)
    write.csv(res, fname, row.names = FALSE)
    print(paste(ada_t+deep_t, 'seconds for dataset:', names(datasets)[i], ',ensemble size:', num_iter))
    results <- rbind(results, res)
  }
}
write.csv(results, './results.txt', row.names = FALSE)
