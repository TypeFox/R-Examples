"boot.boot" <-
function (x, i, Epsilon=0.5)
{
    reps <- 50
    x <- x[i]
    mean <- mean(x)
    boots <- boot(x, meanInBoot, reps)$t
    C1 <- quantile(boots, 0.05)
    C2 <- quantile(boots, 0.95)

    result <- (C1 > -Epsilon & C2 < +Epsilon)
}

