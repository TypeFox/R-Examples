inspect.lm.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000, 1000000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    y <- rnorm(sample.size)
    unpurged.model <- lm(y ~ x)
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.lm', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.lm', ctr, '.rds', sep=''))
  }
}

inspect.glm.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000, 1000000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    y <- as.factor(runif(sample.size) > 0.5)
    z <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- glm(y ~ x + z, family=binomial)
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.glm', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.glm', ctr, '.rds', sep=''))
  }
}

inspect.merMod.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000, 1000000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    y <- rnorm(sample.size)
    z <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- lmer(y ~ x + (1|z))
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.merMod', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.merMod', ctr, '.rds', sep=''))
  }
}

inspect.glmerMod.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000, 1000000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    y <- as.factor(runif(sample.size) > 0.5)
    z <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- glmer(y ~ x + (1|z), family=binomial())
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.glmerMod', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.glmerMod', ctr, '.rds', sep=''))
  }
}

inspect.rpart.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000, 1000000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    w <- rnorm(sample.size)
    y <- x + w + rnorm(sample.size)
    unpurged.model <- rpart(y ~ x + w)
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.rpart', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.rpart', ctr, '.rds', sep=''))
  }
}

inspect.randomForest.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    w <- rnorm(sample.size)
    y <- x + w + rnorm(sample.size)
    unpurged.model <- randomForest(y ~ x + w, ntree=10, nodesize=100)
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.randomForest', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.randomForest', ctr, '.rds', sep=''))
  }
}

inspect.coxph.file.sizes <- function() {
  ctr <- 0
  for (sample.size in c(1000, 10000, 100000)) {
    ctr <- ctr + 1
    x <- rnorm(sample.size)
    y.time <- abs(rnorm(sample.size))
    y.status <- ifelse(runif(sample.size) > 0.5, 0, 1)
    unpurged.model <- survival::coxph(survival::Surv(y.time, y.status) ~ x)
    purged.model <- purge(unpurged.model)
    saveRDS(unpurged.model, paste('unpurged.coxph', ctr, '.rds', sep=''))
    saveRDS(purged.model, paste('purged.coxph', ctr, '.rds', sep=''))
  }
}

run.all.inspect.file.sizes <- function() {
  inspect.lm.file.sizes()
  inspect.glm.file.sizes()
  inspect.merMod.file.sizes()
  inspect.glmerMod.file.sizes()
  inspect.rpart.file.sizes()
  inspect.randomForest.file.sizes()
  inspect.coxph.file.sizes()
}
