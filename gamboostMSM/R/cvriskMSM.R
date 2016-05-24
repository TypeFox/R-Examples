#formulaMSM <- as.formula("Surv(entry, exit, event) ~ bols(x1.12, intercept=F, df=1) +
#                   bols(x1.13, intercept=F, df=1) + bols(x1.14, intercept=F, df=1) +
#                   bols(x1.12.13, intercept=F, df=1) + bols(x1.12.14, intercept=F, df=1) +
#                   bols(x1.13.14, intercept=F, df=1) + bols(x1, intercept=F, df=1)")
#source("C:/Users/Holger Reulen/Desktop/Dropbox/packages/gamboostMSM/R/plloss.R")
#source("C:/Users/Holger Reulen/Desktop/Dropbox/packages/gamboostMSM/R/helpfunctionmultistate1.R")
#source("C:/Users/Holger Reulen/Desktop/Dropbox/packages/gamboostMSMrevision/R/meancentering.R")
#xlist <- list("x1")
#qlist <- list(12, 13, 14, c(12, 13), c(12, 14), c(13, 14))
#ho <- cvriskMSM(m = bm, d = d, formulaMSM = formulaMSM, xlist = xlist, qlist = qlist, k = rep(1:10, each = 60), riskset = riskset)
cvriskMSM <- function(m, d, id, formulaMSM, xlist, qlist, k, riskset){
  riskset.full.data <- riskset$Ri
  full.data <- d
  full.model <- m
  if(length(k) > 1.5){
    index <- sort(unique(k))
    type <- "cv"
    cvpl.matrix <- matrix(nrow = full.model$control$mstop, ncol = length(index), data = 0)
  }else{
    index <- 1:k
    type <- "subsampling"
    cvpl.matrix <- matrix(nrow = full.model$control$mstop, ncol = k, data = 0)
  }  
  for(hi1 in index){
    cat(paste("subsample ", hi1, "\n", sep = ""))
    if(type == "cv"){
      d <- full.data[k != hi1, ]
    }
    if(type == "subsampling"){
      IDs <- sample(unique(full.data[, id]))
      IDs <- IDs[1:round(length(unique(full.data[, id])))]
      d <- subset(full.data, id %in% IDs)
    }
    for(hi2 in 1:length(xlist)){
      x <- xlist[[hi2]]
      for(hi3 in 1:length(qlist)){
        q <- qlist[[hi3]]
        ho <- meancentering(d = d, x = x, q = q, x.name = NULL, q.name = NULL)
        d[, ho$name] <- ho$x.q
      }
    }
    riskset <- buildrisksets(entry = d$entry, exit = d$exit, trans = d$trans, event = d$event,
                        statusinfo = FALSE)
    m <- gamboost(formulaMSM, data = d, family = multistate(Ri = riskset$Ri, Ci = riskset$Ci),
                  control = boost_control(mstop = full.model$control$mstop, trace = FALSE, 
                                          nu = full.model$control$nu))
    lpl.2 <- -m$risk()
    f.hat <- predict(object = m,  newd = full.data, aggregate = "cumsum")
    lpl.1 <- rep(0, m$control$mstop)
    for(hi2 in 1:m$control$mstop){
      lpl.1[hi2] <- sum(plloss(event = full.data$event, f = f.hat[, hi2], Ri = riskset.full.data))
    }
    cvpl.matrix[, hi1] <- lpl.1 - lpl.2
    rm(d)
  }
  mcvpl <- apply(cvpl.matrix, MARGIN = 1, FUN = mean)
  stop.at <- which.max(mcvpl)
  attributes(stop.at)$cvpl.matrix <- cvpl.matrix
  return(stop.at)
}