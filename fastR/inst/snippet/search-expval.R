# expected value of binary search
vals <- 1:8
probs <- 2^(vals-1) / 255
sum(vals*probs)

# expected value of linear search
vals <- 1:255
probs <- rep(1/255,255)
sum(vals*probs)
