
# Q1

id <- seq(1, 6)  # id <- 1:6

age <- c(30, 32, 28, 39, 20, 25)

edu <- rep(0, 6)

class <- factor(rep(c("poor", "middle"), each=3))

# Q2

IndianMothers <- data.frame(id, age, edu, class)

# Q3

ageSummary <- round(c(min=min(age), max=max(age),
                      median=median(age),
                      mean=mean(age), sd=sd(age)),
                    1)
