# list() v. unlist()
(ms <- list(a = 1:3, b = c("maple", "pine")))
unlist(ms)

# Comparison: constructor vs. coercion
identical(list(3), as.list(3))
identical(list(1:3), as.list(1:3))
  
# Missing value
test <- c(7, 5, 3, 9, 1 / 0)
test[2] <- NA; test[3] <- NaN; test[4] <- Inf
test
is.na(test); is.nan(test); is.finite(test); is.infinite(test)
test1 <- c(4.5, 6); test1 == 6; is.na(test1)
test2 <- c(NA, NA); test2 == NA; is.na(test2)