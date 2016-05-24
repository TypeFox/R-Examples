vals <- 0:4
probs <- c(1,4,6,4,1) / 16      # providing probabilities directly
sum(vals * probs)
sum(0:4 * dbinom(0:4,4,0.5))     # using the fact that X is binomial
