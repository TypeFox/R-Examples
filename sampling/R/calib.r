calib<-function (Xs, d, total, q = rep(1, length(d)), method = c("linear", 
    "raking", "truncated", "logit"), bounds = c(low = 0, upp = 10), 
    description = FALSE, max_iter = 500) 
{
    if (any(is.na(Xs)) | any(is.na(d)) | any(is.na(total)) | 
        any(is.na(q))) 
        stop("the input should not contain NAs")
    if (!(is.matrix(Xs) | is.array(Xs))) 
        Xs = as.matrix(Xs)
    if (is.matrix(Xs)) 
        if (length(total) != ncol(Xs)) 
            stop("Xs and total have different dimensions")
    if (is.vector(Xs) & length(total) != 1) 
        stop("Xs and total have different dimensions")
    if (any(is.infinite(q))) 
        stop("there are Inf values in the q vector")
    if (missing(method)) 
        stop("specify a method")
    if (!(method %in% c("linear", "raking", "logit", "truncated"))) 
        stop("the specified method is not in the list")
    if (method %in% c("linear", "raking") & !missing(bounds)) 
        stop("for the linear and raking the bounds are not allowed")
    EPS = .Machine$double.eps
    EPS1 = 1e-06
    n = length(d)
    lambda = as.matrix(rep(0, n))
    lambda1 = ginv(t(Xs * d * q) %*% Xs, tol = EPS) %*% (total - 
        as.vector(t(d) %*% Xs))
    if (method == "linear" | max(abs(lambda1)) < EPS) 
        g = 1 + q * as.vector(Xs %*% lambda1)
    else if (method == "truncated") {
        if (!missing(bounds)) {
            if (bounds[2] <= 1 || bounds[1] >= 1 || bounds[1] > 
                bounds[2]) 
                warning("The conditions low<1<upp are not fulfilled")
        }
        else stop("Specify the bounds")
        g = 1 + q * as.vector(Xs %*% lambda1)
        list = 1:length(g)
        l = 0
        g1 = g
        Xs1 = Xs
        d1 = d
        t2 = total
        list1 = 1:length(g)
        q1 = q
        while (l < max_iter) {
            l = l + 1
            if (any(g < bounds[1]) | any(g > bounds[2])) {
                g[g < bounds[1]] = bounds[1]
                g[g > bounds[2]] = bounds[2]
                list = (1:length(g))[g > bounds[1] & g < bounds[2]]
                if (length(list) != 0) {
                  g1 = g[list]
                  t2 = total - as.vector(t(g[-list] * d[-list]) %*% 
                    Xs[-list, ])
                  Xs1 = Xs[list, ]
                  d1 = d[list]
                  q1 = q[list]
                  list1 = list
                }
            }
            t1 = as.vector(t(d1) %*% Xs1)
            lambda1 = ginv(t(Xs1 * d1 * q1) %*% Xs1, tol = EPS) %*% 
                (t2 - t1)
            if (length(list1) > 1) 
                g1 = 1 + q1 * as.vector(Xs1 %*% lambda1)
            else if (length(list1) == 1) {
                g1 = 1 + q1 * as.vector(as.vector(Xs1) %*% as.vector(lambda1))
            }
            g[list1] = g1
            tr = crossprod(Xs, g * d)
            if (max(abs(tr - total)/total) < EPS1 & all(g >= 
                bounds[1] & g <= bounds[2])) 
                break
        }
        if (l == max_iter) {
            cat("No convergence in", max_iter, "iterations with the given bounds. \n")
            cat("The bounds for the g-weights are:", min(g), 
                " and ", max(g), "\n")
            cat(" and the g-weights are given by g\n")
        }
    }
    else if (method == "raking") {
        lambda = as.matrix(rep(0, ncol(Xs)))
        w1 = as.vector(d * exp(Xs %*% lambda * q))
        for (l in 1:max_iter) {
            phi = t(Xs) %*% w1 - total
            T1 = t(Xs * w1)
            phiprim = T1 %*% Xs
            lambda = lambda - ginv(phiprim, tol = EPS) %*% phi
            w1 = as.vector(d * exp(Xs %*% lambda * q))
            if (any(is.na(w1)) | any(is.infinite(w1))) {
                warning("No convergence")
                g = NULL
                break
            }
            tr = crossprod(Xs, w1)
            if (max(abs(tr - total)/total) < EPS1) 
                break
        }
        if (l == max_iter) {
            warning("No convergence")
            g = NULL
        }
        else g = w1/d
    }
    else if (method == "logit") {
        if (bounds[2] <= 1 || bounds[1] >= 1 || bounds[1] > bounds[2]) 
            stop("The conditions low<1<upp are not fulfilled")
        A = (bounds[2] - bounds[1])/((1 - bounds[1]) * (bounds[2] - 
            1))
        u = rep(1, length(d))
        F = (bounds[1] * (bounds[2] - 1) + bounds[2] * (1 - bounds[1]) * 
            u)/(bounds[2] - 1 + (1 - bounds[1]) * u)
        w1 = as.vector(d * F)
        T = t(Xs * w1)
        phiprim = ginv(T %*% Xs, tol = EPS)
        g = F
        tr = crossprod(Xs, w1)
        if (max(abs(tr - total)/total) > EPS1 | any(g < bounds[1]) | 
            any(g > bounds[2])) {
            lambda1 = rep(0, ncol(Xs))
            list = 1:length(g)
            t2 = total
            Xs1 = Xs
            d1 = d
            g1 = g
            q1 = q
            list1 = 1:length(g)
            for (l in 1:max_iter) {
                if (any(g < bounds[1]) | any(g > bounds[2])) {
                  g[g < bounds[1]] = bounds[1]
                  g[g > bounds[2]] = bounds[2]
                  list = (1:length(g))[g > bounds[1] & g < bounds[2]]
                  if (length(list) != 0) {
                    g1 = g[list]
                    t2 = total - as.vector(t(g[-list] * d[-list]) %*% 
                      Xs[-list, ])
                    Xs1 = Xs[list, ]
                    d1 = d[list]
                    q1 = q[list]
                    list1 = list
                  }
                }
                t1 = as.vector(t(d1) %*% Xs1)
                phi = t(Xs1) %*% as.vector(d1 * g1) - t1
                T = t(Xs1 * as.vector(d1 * g1))
                phiprime = T %*% Xs1
                lambda1 = lambda1 - ginv(phiprime, tol = EPS) %*% 
                  (as.vector(phi) - t2 + t1)
                u = exp(A * (Xs1 %*% lambda1 * q1))
                F = g1 = (bounds[1] * (bounds[2] - 1) + bounds[2] * 
                  (1 - bounds[1]) * u)/(bounds[2] - 1 + (1 - 
                  bounds[1]) * u)
                if (any(is.na(g1))) {
                  warning("no convergence")
                  g1 = g = NULL
                  break
                }
                g[list1] = g1
                tr = crossprod(Xs, g * d)
                if (max(abs(tr - total)/total) < EPS1 & all(g >= 
                  bounds[1] & g <= bounds[2])) 
                  break
            }
            if (l == max_iter) {
                cat("no convergence in", max_iter, "iterations with the given bounds. \n")
                cat("the bounds for the g-weights are:", min(g), 
                  " and ", max(g), "\n")
                cat(" and the g-weights are given by g\n")
                g = NULL
            }
        }
    }
    if (description && !is.null(g)) {
        par(mfrow = c(3, 2), pty = "s")
        hist(g)
        boxplot(g, main = "Boxplot of g")
        hist(d)
        boxplot(d, main = "Boxplot of d")
        hist(g * d)
        boxplot(g * d, main = "Boxplot of w=g*d")
        if (method %in% c("truncated", "raking", "logit")) 
            cat("number of iterations ", l, "\n")
        cat("summary - initial weigths d\n")
        print(summary(d))
        cat("summary - final weigths w=g*d\n")
        print(summary(as.vector(g * d)))
    }
g
}
