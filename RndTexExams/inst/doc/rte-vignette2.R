## ---- eval=TRUE, fig.width= 5--------------------------------------------
library(RndTexExams, quietly = TRUE, warn.conflicts = TRUE)

n.sim.questions <- 15

first.names <- c('John', 'Marcelo','Ricardo', 'Tarcizio')
last.names <- c('Smith', 'Johnson', 'P.')

name.grid <- expand.grid(first.names,last.names)

# Create names with all combinations of first and last name
my.names <- paste(name.grid[,1], name.grid[,2])

set.seed(15)

# A simulated grade of students (TRUE = correct, FALSE = wrong)
correction.mat <- matrix(sample(c(TRUE,FALSE),
                                size = length(my.names)*n.sim.questions,
                                replace = TRUE),nrow = length(my.names))

# Simulate some cheating
idx.cheater.1 <- 5 # std 5 and 10 will have similar correct answers
idx.cheater.2 <- 10 
proportion.to.cheat <- 0.5  # proportion of same correct answers
q.to.cheat <- floor(proportion.to.cheat*n.sim.questions)
correction.mat[idx.cheater.1, ] <-  c(rep(TRUE,q.to.cheat),
                                      rep(FALSE,n.sim.questions-q.to.cheat))
correction.mat[idx.cheater.2, ] <- correction.mat[idx.cheater.1, ]

# bind names and correction matrix
df.grade <- cbind(data.frame(exam.names = my.names),correction.mat)

# test for cheating
test.cheating.out <- rte.test.cheating(df.grade)


