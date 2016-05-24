vals <- 0:50
probs.x <- dnbinom(vals,3,0.5)
probs.y <- dnbinom(vals,3,0.2)
var.x <- sum(vals^2 * probs.x) - sum(vals * probs.x)^2; var.x
var.y <- sum(vals^2 * probs.y) - sum(vals * probs.y)^2; var.y

# better approximation using more terms:
vals <- 0:500
probs.x <- dnbinom(vals,3,0.5)
probs.y <- dnbinom(vals,3,0.2)
var.x <- sum(vals^2 * probs.x) - sum(vals * probs.x)^2; var.x
var.y <- sum(vals^2 * probs.y) - sum(vals * probs.y)^2; var.y
