simCR <- function(n, bhr, beta.x.12, beta.x.13, beta.x.14){
  bhr.11 <- function(t){return(0)}; eta.1.1 <- function(x.i, t){return(0*sum(x.i))}
  bhr.12 <- function(t){return(bhr)}
  bhr.13 <- function(t){return(bhr)}
  bhr.21 <- function(t){return(0)}; eta.2.1 <- function(x.i, t){return(0*sum(x.i))}
  bhr.22 <- function(t){return(0)}; eta.2.2 <- function(x.i, t){return(0*sum(x.i))}
  bhr.23 <- function(t){return(0)}; eta.2.3 <- function(x.i, t){return(0*sum(x.i))}
  bhr.31 <- function(t){return(0)}; eta.3.1 <- function(x.i, t){return(0*sum(x.i))}
  bhr.32 <- function(t){return(0)}; eta.3.2 <- function(x.i, t){return(0*sum(x.i))}
  bhr.33 <- function(t){return(0)}; eta.3.3 <- function(x.i, t){return(0*sum(x.i))}
  eta.1.2 <- function(x.i, t){
    eta <- beta.x.12 * x.i
    return(eta)}
  eta.1.3 <- function(x.i, t){
    eta <- beta.x.13 * x.i
    return(eta)}
  bhr.14 <- function(t){return(bhr)}
  bhr.24 <- function(t){return(0)}; eta.2.4 <- function(x.i, t){return(0*sum(x.i))}
  bhr.34 <- function(t){return(0)}; eta.3.4 <- function(x.i, t){return(0*sum(x.i))}
  bhr.41 <- function(t){return(0)}; eta.4.1 <- function(x.i, t){return(0*sum(x.i))}
  bhr.42 <- function(t){return(0)}; eta.4.2 <- function(x.i, t){return(0*sum(x.i))}
  bhr.43 <- function(t){return(0)}; eta.4.3 <- function(x.i, t){return(0*sum(x.i))}
  bhr.44 <- function(t){return(0)}; eta.4.4 <- function(x.i, t){return(0*sum(x.i))}
  eta.1.4 <- function(x.i, t){
    eta <- beta.x.14 * x.i
    return(eta)
  }
  mpl <- list(from.1 = list(from = 1,
                            all.to = c(2, 3, 4),
                            all.bhr = list(bhr.11, bhr.12, bhr.13, bhr.14),
                            eta = list(to.1 = eta.1.1, to.2 = eta.1.2, to.3 = eta.1.3, to.4 = eta.1.4)),
              from.2 = list(from = 2,
                            all.to = NULL,
                            all.bhr = list(bhr.21, bhr.22, bhr.23, bhr.24),
                            eta = list(to.1 = eta.2.1, to.2 = eta.2.2, to.3 = eta.2.3, to.4 = eta.2.4)),
              from.3 = list(from = 3,
                            all.to = NULL,
                            all.bhr = list(bhr.31, bhr.32, bhr.33, bhr.34),
                            eta = list(to.1 = eta.3.1, to.2 = eta.3.2, to.3 = eta.3.3, to.4 = eta.3.4)),
              from.4 = list(from = 4,
                            all.to = NULL,
                            all.bhr = list(bhr.41, bhr.42, bhr.43, bhr.44),
                            eta = list(to.1 = eta.3.1, to.2 = eta.3.2, to.3 = eta.3.3, to.4 = eta.4.4)))
  x <- rnorm(n) + 1
  X <- matrix(ncol = 1, nrow = n, data = x)
  colnames(X) <- c("x")
  max.time = 100
  all.first.from <- rep(1, times=n)
  p <- ncol(X)
  histories <- NULL
  for (i in 1:n) {
    history.i <- simsinglehistory(first.entry = 0, first.from = all.first.from[i],
                                  max.time, mpl, x.i = X[i, ])
    if (!is.null(nrow(history.i))) {
      history.i <- cbind(history.i, rep(i, 1))
      for (x.index in 1:p) {
        history.i <- cbind(history.i, rep(X[i, x.index], 1))
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
    cens <- sample(c(TRUE, FALSE), size = 1, prob = c(0.15, 0.85))
    if(cens){histories$to[i] <- "cens"}
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
  d <- d[, c("id", "entry", "exit", "from", "to", "trans", "event", "x")]
  hi <- which(d$trans == 12)
  x.12.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.12.mean
  ho[which(!(d$trans == 12))] <- 0
  d$x.12 <- ho
  hi <- which(d$trans == 13)
  x.13.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.13.mean
  ho[which(!(d$trans == 13))] <- 0
  d$x.13 <- ho
  hi <- which(d$trans %in% c(12, 13))
  x.12.13.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.12.13.mean
  ho[which(!(d$trans %in% c(12, 13)))] <- 0
  d$x.12.13 <- ho
  hi <- which(d$trans == 14)
  x.14.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.14.mean
  ho[which(!(d$trans == 14))] <- 0
  d$x.14 <- ho
  hi <- which(d$trans %in% c(12, 14))
  x.12.14.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.12.14.mean
  ho[which(!(d$trans %in% c(12, 14)))] <- 0
  d$x.12.14 <- ho
  hi <- which(d$trans %in% c(13, 14))
  x.13.14.mean <- mean(d[hi, "x"])
  ho <- d[, "x"] - x.13.14.mean
  ho[which(!(d$trans %in% c(13, 14)))] <- 0
  d$x.13.14 <- ho
  x.mean <- mean(d$x)
  d$x <- d$x - x.mean 
  attributes(d)$x.means <- c(x.12.mean, x.13.mean, x.14.mean, 
                            x.12.13.mean, x.12.14.mean, x.13.14.mean, x.mean)
  return(d)}