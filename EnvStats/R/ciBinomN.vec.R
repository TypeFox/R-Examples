ciBinomN.vec <-
function (half.width, p.hat.or.p1.hat, p2.hat, conf.level, sample.type, 
    ratio, ci.method, correct, n.or.n1.min, n.or.n1.max, tol.half.width, 
    tol.p.hat, tol, maxiter) 
{
    if (sample.type == "two.sample") {
        n2 <- round(n.or.n1.min * ratio)
        dum.list <- ciBinomHalfWidth(n.or.n1 = n.or.n1.min, p.hat.or.p1.hat = p.hat.or.p1.hat, 
            n2 = n2, p2.hat = p2.hat, conf.level = conf.level, 
            sample.type = "two.sample", ci.method = ci.method, 
            correct = correct, warn = FALSE)
        if ((dum.list$half.width - half.width) <= tol.half.width & 
            abs(p.hat.or.p1.hat - dum.list$p1.hat) <= tol.p.hat & 
            abs(p2.hat - dum.list$p2.hat) <= tol.p.hat) {
            ret.val <- c(n1 = dum.list$n1, p1.hat = dum.list$p1.hat, 
                n2 = dum.list$n2, p2.hat = dum.list$p2.hat, half.width = dum.list$half.width)
        }
        else {
            fcn.for.root.two <- function(n.weird, p1.hat.weird, 
                p2.hat.weird, ratio.weird, conf.weird, ci.method.weird, 
                correct.weird) {
                half.width - ciBinomHalfWidth(n.or.n1 = round(n.weird), 
                  p.hat.or.p1.hat = p1.hat.weird, n2 = round(round(n.weird) * 
                    ratio), p2.hat = p2.hat.weird, conf.level = conf.weird, 
                  sample.type = "two.sample", ci.method = ci.method.weird, 
                  correct = correct.weird, warn = FALSE)$half.width
            }
            f.lower <- fcn.for.root.two(n.weird = n.or.n1.min, 
                p1.hat.weird = p.hat.or.p1.hat, p2.hat.weird = p2.hat, 
                ratio.weird = ratio, conf.weird = conf.level, 
                ci.method.weird = ci.method, correct.weird = correct)
            f.upper <- fcn.for.root.two(n.weird = n.or.n1.max, 
                p1.hat.weird = p.hat.or.p1.hat, p2.hat.weird = p2.hat, 
                ratio.weird = ratio, conf.weird = conf.level, 
                ci.method.weird = ci.method, correct.weird = correct)
            old.n <- uniroot(fcn.for.root.two, p1.hat.weird = p.hat.or.p1.hat, 
                p2.hat.weird = p2.hat, ratio.weird = ratio, conf.weird = conf.level, 
                ci.method.weird = ci.method, correct.weird = correct, 
                lower = n.or.n1.min, upper = n.or.n1.max, f.lower = f.lower, 
                f.upper = f.upper, tol = tol, maxiter = maxiter)$root
            old.n <- ceiling(old.n)
            new.n.vec <- n.or.n1.min:old.n
            dum.list <- ciBinomHalfWidth(n.or.n1 = new.n.vec, 
                p.hat.or.p1.hat = p.hat.or.p1.hat, n2 = round(ratio * 
                  new.n.vec), p2.hat = p2.hat, conf.level = conf.level, 
                sample.type = "two.sample", ci.method = ci.method, 
                correct = correct, warn = FALSE)
            index <- (dum.list$half.width - half.width) <= tol.half.width & 
                abs(p.hat.or.p1.hat - dum.list$p1.hat) <= tol.p.hat & 
                abs(p2.hat - dum.list$p2.hat) <= tol.p.hat
            if (!any(index)) {
                while (!any(index) & max(new.n.vec) < n.or.n1.max) {
                  old.n <- max(new.n.vec) + 1
                  new.n.vec <- old.n:min(old.n + 100, n.or.n1.max)
                  dum.list <- ciBinomHalfWidth(n.or.n1 = new.n.vec, 
                    p.hat.or.p1.hat = p.hat.or.p1.hat, n2 = round(ratio * 
                      new.n.vec), p2.hat = p2.hat, conf.level = conf.level, 
                    sample.type = "two.sample", ci.method = ci.method, 
                    correct = correct, warn = FALSE)
                  index <- (dum.list$half.width - half.width) <= 
                    tol.half.width & abs(p.hat.or.p1.hat - dum.list$p1.hat) <= 
                    tol.p.hat & abs(p2.hat - dum.list$p2.hat) <= 
                    tol.p.hat
                }
            }
            index <- ((1:length(new.n.vec))[index])[1]
            ret.val <- c(n1 = dum.list$n1[index], p1.hat = dum.list$p1.hat[index], 
                n2 = dum.list$n2[index], p2.hat = dum.list$p2.hat[index], 
                half.width = dum.list$half.width[index])
        }
    }
    else {
        dum.list <- ciBinomHalfWidth(n.or.n1 = n.or.n1.min, p.hat.or.p1.hat = p.hat.or.p1.hat, 
            conf.level = conf.level, sample.type = "one.sample", 
            ci.method = ci.method, correct = correct, warn = FALSE)
        if ((dum.list$half.width - half.width) <= tol.half.width & 
            abs(p.hat.or.p1.hat - dum.list$p.hat) <= tol.p.hat) {
            ret.val <- c(n = n.or.n1.min, p.hat = dum.list$p.hat, 
                half.width = dum.list$half.width)
        }
        else {
            fcn.for.root.one <- function(n.weird, p.hat.weird, 
                conf.weird, ci.method.weird, correct.weird) {
                half.width - ciBinomHalfWidth(n.or.n1 = round(n.weird), 
                  p.hat.or.p1.hat = p.hat.weird, conf.level = conf.weird, 
                  sample.type = "one.sample", ci.method = ci.method.weird, 
                  correct = correct.weird, warn = FALSE)$half.width
            }
            f.lower <- fcn.for.root.one(n.weird = n.or.n1.min, 
                p.hat.weird = p.hat.or.p1.hat, conf.weird = conf.level, 
                ci.method.weird = ci.method, correct.weird = correct)
            f.upper <- fcn.for.root.one(n.weird = n.or.n1.max, 
                p.hat.weird = p.hat.or.p1.hat, conf.weird = conf.level, 
                ci.method.weird = ci.method, correct.weird = correct)
            old.n <- uniroot(fcn.for.root.one, p.hat.weird = p.hat.or.p1.hat, 
                conf.weird = conf.level, ci.method.weird = ci.method, 
                correct.weird = correct, lower = n.or.n1.min, 
                upper = n.or.n1.max, f.lower = f.lower, f.upper = f.upper, 
                tol = tol, maxiter = maxiter)$root
            old.n <- ceiling(old.n)
            new.n.vec <- n.or.n1.min:old.n
            dum.list <- ciBinomHalfWidth(n.or.n1 = new.n.vec, 
                p.hat.or.p1.hat = p.hat.or.p1.hat, conf.level = conf.level, 
                sample.type = "one.sample", ci.method = ci.method, 
                correct = correct, warn = FALSE)
            index <- (dum.list$half.width - half.width) <= tol.half.width & 
                abs(p.hat.or.p1.hat - dum.list$p.hat) <= tol.p.hat
            if (!any(index)) {
                while (!any(index) & max(new.n.vec) < n.or.n1.max) {
                  old.n <- max(new.n.vec) + 1
                  new.n.vec <- old.n:min(old.n + 100, n.or.n1.max)
                  dum.list <- ciBinomHalfWidth(n.or.n1 = new.n.vec, 
                    p.hat.or.p1.hat = p.hat.or.p1.hat, conf.level = conf.level, 
                    sample.type = "one.sample", ci.method = ci.method, 
                    correct = correct, warn = FALSE)
                  index <- (dum.list$half.width - half.width) <= 
                    tol.half.width & abs(p.hat.or.p1.hat - dum.list$p.hat) <= 
                    tol.p.hat
                }
            }
            index <- ((1:length(new.n.vec))[index])[1]
            ret.val <- c(n = dum.list$n[index], p.hat = dum.list$p.hat[index], 
                half.width = dum.list$half.width[index])
        }
    }
    ret.val
}
