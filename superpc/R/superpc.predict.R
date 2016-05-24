superpc.predict <- function (object, data, newdata, threshold, n.components = 3, 
    prediction.type = c("continuous", "discrete", "nonzero"), 
    n.class = 2) 
{
    this.call <- match.call()
    prediction.type <- match.arg(prediction.type)
    if (n.class > 3) {
        stop("Maximum number of survival classes is 3")
    }
    which.features <- (abs(object$feature.scores) >= threshold)
    x.sml <- data$x[which.features, ]
    n.pc <- n.components
    x.sml.svd <- mysvd(x.sml, n.components = n.components)
    if (prediction.type == "nonzero") {
        if (!is.null(data$featurenames)) {
            out <- data$featurenames[which.features]
        }
        else {
            out <- (1:nrow(data$x))[which.features]
        }
    }
    if (prediction.type == "continuous" | prediction.type == 
        "discrete") {
        xtemp = newdata$x[which.features, ]
        xtemp = t(scale(t(xtemp), center = x.sml.svd$feature.means, 
            scale = F))
        scal = apply(scale(abs(x.sml.svd$u), center = F, scale = x.sml.svd$d), 
            2, sum)
        cur.v <- scale(t(xtemp) %*% x.sml.svd$u, center = FALSE, 
            scale = scal * x.sml.svd$d)
        xtemp0 = data$x[which.features, ]
        xtemp0 = t(scale(t(xtemp0), center = x.sml.svd$feature.means, 
            scale = F))
        cur.v0 <- scale(t(xtemp0) %*% x.sml.svd$u, center = FALSE, 
            scale = scal * x.sml.svd$d)
    }
    result <- superpc.fit.to.outcome(object, data, cur.v0, print = FALSE)$results
    if (object$type == "survival") {
        coef = result$coef
    }
    if (object$type == "regression") {
        coef = result$coef[-1]
    }
    if (prediction.type == "continuous") {
        out <- scale(cur.v, center = FALSE, scale = sign(coef))
        v.pred.1df = apply(scale(out, center = FALSE, scale = 1/abs(coef)), 
            1, sum)
    }
    else if (prediction.type == "discrete") {
        out0 <- scale(cur.v0, center = FALSE, scale = sign(coef))
        v.pred0.1df = apply(scale(out0, center = FALSE, scale = 1/abs(coef)), 
            1, sum)
        out <- scale(cur.v, center = FALSE, scale = sign(coef))
        v.pred.1df = apply(scale(out, center = FALSE, scale = 1/abs(coef)), 
            1, sum)
        for (j in 1:ncol(out)) {
           # br = quantile(cur.v0[, j], (0:n.class)/n.class)
             br = quantile(out0[, j], (0:n.class)/n.class) ## yp
           # out[, j] <- cut(out[, j], breaks = br, n.class, labels = FALSE)
             out[,j] = ifelse(out[,j] <= br[2], 1, 2) ## yp
          #  out[is.na(out[, j]), j] <- 1
        }
        br = quantile(v.pred0.1df, (0:n.class)/n.class)
       # v.pred.1df <- cut(v.pred.1df, breaks = br, labels = FALSE)
       # v.pred.1df[is.na(v.pred.1df)] <- 1
        v.pred.1df = ifelse(v.pred.1df <= br[2], 1, 2)  ## yp
    }
    if (is.matrix(out)) {
        dimnames(out) = list(NULL, rep(prediction.type, ncol(out)))
    }
    junk <- list(v.pred = out, u = x.sml.svd$u, d = x.sml.svd$d, 
        which.features = which.features, v.pred.1df = v.pred.1df, 
        n.components = n.pc, coef = result$coef, call = this.call, 
        prediction.type = prediction.type)
    return(junk)
}

