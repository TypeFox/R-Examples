# Create a factor and change its levels
state.region; table(state.region)
dat <- c("bold", "full", "half", "full")
man2 <- man <- factor(x = dat); levels(man); nlevels(man)
man; attributes(man)
levels(man2) <- c("bold", "health", "half")
levels(man2)[3] <- "quasi"

# Converting a factor with letters into a character string
man4 <- man3 <- unclass(man)  # Correct to an integer vector
class(man4) <- "factor"; identical(man, man4)  # back to a factor
man5 <- as.character(man)  # Correct to a character string
man6 <- as.numeric(man)
man7 <- as.numeric(as.character(man))  # Warning: Coercion creates NAs

# Convert a factor with digits into a numeric vector
num <- factor(x = c(2, 3, 10, 15, 50, 10))
num2 <- as.character(num)
num3 <- as.numeric(num)
num4 <- as.numeric(as.character(num))             # Correct and the best
levels(num)[c(1,2,3,4,5,3,4,5)]                   # Learn subscript
num4 <- as.numeric(levels(num)[as.numeric(num)])  # Correct but lengthy
num4 <- as.numeric(levels(num)[num])              # Correct and 2nd best
num; num2; num3; num4

# Add new elements to a factor vector
bob1 <- num2; bob2 <- num4; bob3 <- num
bob1[7:8] <- c("60", "10")  # add two elements to a character vector 
bob2[7:8] <- c(60, 10)      # add two elements to a numerical vector 
bob3[7] <- 60               # warning: invalid factor level, NA generated
levels(bob3) <- c(levels(num), "60")
bob3[7:8] <- c(60, 10)      # Add two elements to a factor vector
bob1; bob2; bob3

# How to create an ordered factor?
h.F <- factor(x = dat, levels = c("bold", "full", "half"), ordered = FALSE)
h.T <- factor(x = dat, levels = c("bold", "half", "full"), ordered = TRUE)
h.T <- ordered(x = dat, levels = c("bold", "half", "full"))
is.factor(h.F); is.ordered(h.F)
h.T

# How to order an existing unordered factor?
h.FT <- factor(x = h.F, levels = levels(h.F)[c(2, 3, 1)], ordered = TRUE)
n.FT <- factor(x = num, levels = levels(num)[c(2,3,1,4,5)], ordered = TRUE)
h.ss <- relevel(x = h.F, ref = "full")  # change one only
h.FT; n.FT; h.ss