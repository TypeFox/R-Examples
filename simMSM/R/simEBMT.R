simEBMT <- function(n){
    bhr <- 0.2
    bhr.11 <- function(t){return(0)}; eta.1.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.12 <- function(t){return(2.0*bhr)}
    bhr.13 <- function(t){return(2.5*bhr)}
    bhr.14 <- function(t){return(0)}; eta.1.4 <- function(x.i, t){return(0*sum(x.i))}
    bhr.15 <- function(t){return(0.5*bhr)}
    bhr.16 <- function(t){return(0.5*bhr)}
    bhr.21 <- function(t){return(0)}; eta.2.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.22 <- function(t){return(0)}; eta.2.2 <- function(x.i, t){return(0*sum(x.i))}
    bhr.23 <- function(t){return(0)}; eta.2.3 <- function(x.i, t){return(0*sum(x.i))}
    bhr.24 <- function(t){return(2.0*bhr)}
    bhr.25 <- function(t){return(1.0*bhr)}
    bhr.26 <- function(t){return(0.5*bhr)}
    bhr.31 <- function(t){return(0)}; eta.3.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.32 <- function(t){return(0)}; eta.3.2 <- function(x.i, t){return(0*sum(x.i))}
    bhr.33 <- function(t){return(0)}; eta.3.3 <- function(x.i, t){return(0*sum(x.i))}
    bhr.34 <- function(t){return(3.0*bhr)}
    bhr.35 <- function(t){return(0.5*bhr)}
    bhr.36 <- function(t){return(2.0*bhr)}
    bhr.41 <- function(t){return(0)}; eta.4.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.42 <- function(t){return(0)}; eta.4.2 <- function(x.i, t){return(0*sum(x.i))}
    bhr.43 <- function(t){return(0)}; eta.4.3 <- function(x.i, t){return(0*sum(x.i))}
    bhr.44 <- function(t){return(0)}; eta.4.4 <- function(x.i, t){return(0*sum(x.i))}
    bhr.45 <- function(t){return(bhr)}
    bhr.46 <- function(t){return(bhr)}
    bhr.51 <- function(t){return(0)}; eta.5.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.52 <- function(t){return(0)}; eta.5.2 <- function(x.i, t){return(0*sum(x.i))}
    bhr.53 <- function(t){return(0)}; eta.5.3 <- function(x.i, t){return(0*sum(x.i))}
    bhr.54 <- function(t){return(0)}; eta.5.4 <- function(x.i, t){return(0*sum(x.i))}
    bhr.55 <- function(t){return(0)}; eta.5.5 <- function(x.i, t){return(0*sum(x.i))}
    bhr.56 <- function(t){return(0)}; eta.5.6 <- function(x.i, t){return(0*sum(x.i))}
    bhr.61 <- function(t){return(0)}; eta.6.1 <- function(x.i, t){return(0*sum(x.i))}
    bhr.62 <- function(t){return(0)}; eta.6.2 <- function(x.i, t){return(0*sum(x.i))}
    bhr.63 <- function(t){return(0)}; eta.6.3 <- function(x.i, t){return(0*sum(x.i))}
    bhr.64 <- function(t){return(0)}; eta.6.4 <- function(x.i, t){return(0*sum(x.i))}
    bhr.65 <- function(t){return(0)}; eta.6.5 <- function(x.i, t){return(0*sum(x.i))}
    bhr.66 <- function(t){return(0)}; eta.6.6 <- function(x.i, t){return(0*sum(x.i))}
    eta.1.2 <- function(x.i, t){
        beta <- c(0, -0.4, 0.4, 0.5, 0, 0)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.1.3 <- function(x.i, t){
        beta <- c(0, -0.3, 0, 0, 0, 0)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.1.5 <- function(x.i, t){
        beta <- c(0, 0.4, 0.4, 0, 0, 0)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.1.6 <- function(x.i, t){
        beta <- c(0, -0.6, -0.4, -0.5, 0.8, 0.9)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.2.4 <- function(x.i, t){
        beta <- c(0, -0.3, 0, 0, 0.3, 0.5)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.2.5 <- function(x.i, t){
        beta <- c(0.4, 0.3, 0, 0, -0.3, 0)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.2.6 <- function(x.i, t){
        beta <- c(0, 0, -0.8, -1, 0, 1.5)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.3.4 <- function(x.i, t){
        beta <- c(0, 0, 0.5, 0.9, -0.4, -0.3)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.3.5 <- function(x.i, t){
        beta <- c(-0.4, 0, -0.3, -0.6, 0, 0.4)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.3.6 <- function(x.i, t){
        beta <- c(0, 0.3, -0.6, 0, 0, 0.5)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.4.5 <- function(x.i, t){
        beta <- c(-0.3, 0, 0, -0.4, 0.4, 0.3)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    eta.4.6 <- function(x.i, t){
        beta <- c(0.6, 0, -0.4, -0.4, 0.8, 1.3)
        eta <- as.numeric(t(beta)%*%x.i)
        return(eta)}
    mpl <- list(from.1 = list(from = 1,
                all.to = c(2, 3, 5, 6),
                all.bhr = list(bhr.11, bhr.12, bhr.13, bhr.14, bhr.15, bhr.16),
                eta = list(to.1 = eta.1.1, to.2 = eta.1.2, to.3 = eta.1.3,
                to.4 = eta.1.4, to.5 = eta.1.5, to.6 = eta.1.6)),
                from.2 = list(from = 2,
                all.to = c(4, 5, 6),
                all.bhr = list(bhr.21, bhr.22, bhr.23, bhr.24, bhr.25, bhr.26),
                eta = list(to.1 = eta.2.1, to.2 = eta.2.2, to.3 = eta.2.3,
                to.4 = eta.2.4, to.5 = eta.2.5, to.6 = eta.2.6)),
                from.3 = list(from = 3,
                all.to = c(4, 5, 6),
                all.bhr = list(bhr.31, bhr.32, bhr.33, bhr.34, bhr.35, bhr.36),
                eta = list(to.1 = eta.3.1, to.2 = eta.3.2, to.3 = eta.3.3,
                to.4 = eta.3.4, to.5 = eta.3.5, to.6 = eta.3.6)),
                from.4 = list(from = 4,
                all.to = c(5, 6),
                all.bhr = list(bhr.41, bhr.42, bhr.43, bhr.44, bhr.45, bhr.46),
                eta = list(to.1 = eta.4.1, to.2 = eta.4.2, to.3 = eta.4.3,
                to.4 = eta.4.4, to.5 = eta.4.5, to.6 = eta.4.6)),
                from.5 = list(from = 5,
                all.to = NULL,
                all.bhr = list(bhr.51, bhr.52, bhr.53, bhr.54, bhr.55, bhr.56),
                eta = list(to.1 = eta.5.1, to.2 = eta.5.2, to.3 = eta.5.3,
                to.4 = eta.5.4, to.5 = eta.5.5, to.6 = eta.5.6)),
                from.6 = list(from = 6,
                all.to = NULL,
                all.bhr = list(bhr.61, bhr.62, bhr.63, bhr.64, bhr.65, bhr.66),
                eta = list(to.1 = eta.6.1, to.2 = eta.6.2, to.3 = eta.6.3,
                to.4 = eta.6.4, to.5 = eta.6.5, to.6 = eta.6.6)))
    x1 <- sample(size=n, c(0, 1), prob=c(0.76, 0.24), replace = T)
    x2 <- sample(size=n, c(0, 1), prob=c(0.76, 0.24), replace = T)
    x3 <- sample(size=n, c(0, 1, 2), prob=c(0.28, 0.39, 0.33), replace = T)
    x4 <- sample(size=n, c(0, 1, 2), prob=c(0.24, 0.53, 0.23), replace = T)
    X <- model.matrix(~ -1 + x1 + x2 + I(x3==1) + I(x3==2) + I(x4==1) + I(x4==2))
    X <- cbind(x1, x2, I(x3==1), I(x3==2), I(x4==1), I(x4==2))
    colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6")
    max.time = 100
    hf <- function(x, k) {
        if (!is.null(x$all.to)) {
            return(x[[k]])
        }
    }
    all.possible.from.states <- as.numeric(do.call(c, lapply(mpl,
                                                             FUN = hf, k = 1)))
    all.first.from <- rep(1, times=n)
    p <- ncol(X)
    histories <- NULL
    for (i in 1:n) {
        history.i <- simsinglehistory(first.entry = 0, first.from = all.first.from[i],
                                      max.time, mpl, x.i = X[i, ])
        if (!is.null(nrow(history.i))) {
            history.i <- cbind(history.i, rep(i, 1))
            for (x.index in 1:p) {
                history.i <- cbind(history.i, rep(X[i, x.index],
                                                  1))
            }
        }
        else {
            history.i <- cbind(history.i, rep(i, nrow(history.i)))
            for (x.index in 1:p) {
                history.i <- cbind(history.i, rep(X[i, x.index],
                                                  nrow(history.i)))
            }
        }
        histories <- rbind(histories, history.i)
        rm(history.i)
    }
    histories.as.list <- list(id = histories[, 5], entry = histories[, 1],
                              exit = histories[, 2], from = histories[, 3],
                              to = histories[, 4])
    for (x.index in 1:p) {
        histories.as.list[[5 + x.index]] <- histories[, 5 + x.index]
        names(histories.as.list)[5 + x.index] <- colnames(X)[x.index]
    }
    histories <- data.frame(histories.as.list)
    rm(histories.as.list)
    ## censoring:
    for(i in 1:nrow(histories)){
        if(histories$from[i] == "1"){
            cens <- sample(c(TRUE, FALSE), size = 1, prob=c(0.15, 0.85))
            if(cens){histories$to[i] <- "cens"}
        }
        if(histories$from[i] == "2"){
            cens <- sample(c(TRUE, FALSE), size = 1, prob=c(0.51, 0.49))
            if(cens){histories$to[i] <- "cens"}
        }
        if(histories$from[i] == "3"){
            cens <- sample(c(TRUE, FALSE), size = 1, prob=c(0.25, 0.75))
            if(cens){histories$to[i] <- "cens"}
        }
        if(histories$from[i] == "4"){
            cens <- sample(c(TRUE, FALSE), size = 1, prob=c(0.63, 0.37))
            if(cens){histories$to[i] <- "cens"}
        }
    }
    ## delete transitions (if there are any) after first censoring (if there is one):
    IDs <- unique(histories$id)
    for(i in IDs){
        hi <- which(histories$id == i)
        ho <- histories[hi, ]$to == "cens"
        if(any(ho)){
            if(min(hi[ho]) < max(hi)){
                histories <- histories[-hi[which(hi > min(hi[ho]))], ]
            }
        }
    }
    histories$event <- as.numeric(histories$to != "cens")
    ## long format:
    all.to.list <- list(from.1 = mpl$from.1$all.to,
                        from.2 = mpl$from.2$all.to,
                        from.3 = mpl$from.3$all.to,
                        from.4 = mpl$from.4$all.to)
    d <- histories[1, ]
    for(i in 1:nrow(histories)){
        all.to <- all.to.list[[as.numeric(histories$from[i])]]
        add <- rbind(histories[i, ], histories[i, ], histories[i, ], histories[i, ])
        add <- add[1:length(all.to), ]
        reminder.to <- histories[i, "to"]
        reminder.event <- histories[i, "event"]
        add[, "to"] <- all.to
        add[, "event"] <- 0
        if(reminder.event > 0.5){
            add[which(all.to == as.numeric(reminder.to)), "event"] <- 1
        }
        d <- rbind(d, add)
    }
    d <- d[-1, ]
    ## transition indicator:
    d$trans <- paste(d$from, d$to, sep="")
    ## reorder columns:
    d <- d[, c("id", "entry", "exit", "from", "to", "trans", "event", "x1", "x2", "x3", "x4", "x5", "x6")]
    if(any(d$exit - d$entry == 0)){
        d$exit[which(d$exit - d$entry == 0)] <- d$exit[which(d$exit - d$entry == 0)] + 0.000001
        print("problem of equal entry and exit times occurred!")
    }
    centers <- rep(0, times = 186); count <- 1
    ## x1:
    x <- "x1"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## x2:
    x <- "x2"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## x3:
    x <- "x3"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## x4:
    x <- "x4"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## x5:
    x <- "x5"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## x6:
    x <- "x6"
    ## from 1:
    ho <- meancentering(d = d, x = x, q = "12"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13"), q.name = "12.13"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15"), q.name = "12.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "16"), q.name = "12.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15"), q.name = "13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "16"), q.name = "13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("15", "16"), q.name = "15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15"), q.name = "12.13.15"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "16"), q.name = "12.13.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "15", "16"), q.name = "12.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("13", "15", "16"), q.name = "13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("12", "13", "15", "16"), q.name = "12.13.15.16"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 2:
    ho <- meancentering(d = d, x = x, q = "24"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25"), q.name = "24.25"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "26"), q.name = "24.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("25", "26"), q.name = "25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("24", "25", "26"), q.name = "24.25.26"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 3:
    ho <- meancentering(d = d, x = x, q = "34"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35"), q.name = "34.35"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "36"), q.name = "34.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("35", "36"), q.name = "35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("34", "35", "36"), q.name = "34.35.36"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ## from 4:
    ho <- meancentering(d = d, x = x, q = "45"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = "46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    ho <- meancentering(d = d, x = x, q = c("45", "46"), q.name = "45.46"); d[, ho$name] <- ho$x.q; centers[count] <- ho$center; count <- count + 1
    rm(count)
    attributes(d)$center <- centers
    return(d)}
