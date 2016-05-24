HItest <-
function (class, MLE, thresholds = c(2, 1)) 
{
    LLC <- apply(class[, 1:6], 1, max)
    LLM <- MLE[, 3]
    dAIC <- (2 * 1 - 2 * LLC) - (2 * 2 - 2 * LLM)
    c1 <- class$LLD > thresholds[1]
    c2 <- (LLM - LLC) < thresholds[2]
    data.frame(S = MLE[, 1], H = MLE[, 2], Best.class = class$Best, 
        LL.class = LLC, LLD.class = class$LLD, LL.max = MLE[, 
            3], dAIC, c1, c2)
}
