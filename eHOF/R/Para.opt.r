Para_opt <- function (
 		resp, 
		model=NULL, 
		punctual = FALSE,
		...) {
  if(is.null(model)) model <- pick.model(resp, gam=FALSE, ...)
  M <- resp$M
  x <- seq(resp$range[1] - diff(resp$range), resp$range[2] + diff(resp$range),length.out=10000)
  pred <- predict.HOF(resp, newdata = x, model = model)
  HOFfun <- function(resp, x, model) predict(resp, newdata = x, M = M, model = model)
  HOFfun3 <- function(x, y, resp) abs(y - predict(resp, newdata = x, M = M, model = model))
  
  if (model == "I") {
      opt <- NA  # alternative: resp$range
      top <- fitted(resp, 'I')[1]
      mini <- top
      pess <- NA
  }

  if (model == "II") {
      tmp <- optimize(HOFfun, resp$range, resp = resp, model = model, maximum = TRUE)
      opt <- as.numeric(tmp$maximum)
      top <- as.numeric(tmp$objective)
      tmp <- optimize(HOFfun, resp$range, resp = resp, model = model, maximum = FALSE)
      pess <- as.numeric(tmp$minimum)
      mini <- as.numeric(tmp$objective)
  }

  if (model == "III") {
    top <- optimize(HOFfun, interval=resp$range, model = model, resp = resp, maximum = TRUE)$objective
    # optimum as range between 9/10 of the highest response plateau and the gradient edge
    opt <- optimize(HOFfun3, resp$range, y = top * 9/10, resp = resp, maximum = FALSE)$minimum
  if(resp$models$III$par['a'] < 0) 
		opt <- if(punctual) opt - (opt - min(resp$range))/2	else c(opt.min = min(resp$range), opt.max = opt) else
	if(resp$models$III$par['a'] >= 0) 
		opt <- if(punctual) opt + (max(resp$range) - opt)/2	else c(opt.min = opt, opt.max = max(resp$range)) 
#  else if(resp$models$III$par['a'] >= 0 & resp$models$III$par['c'] >= 0)  else warning('Check optimum estimation code.')
  pess <- optimize(HOFfun, model=model, resp$range, resp = resp, maximum = FALSE)
  mini <- pess$objective
  if (predict(resp, newdata = min(resp$range), M = M, model = model) > predict(resp, newdata = max(resp$range), M = M, model = model))  pess <- c(pess$minimum, max(resp$range)) else pess <- c(min(resp$range), pess$minimum)
    }

   if (model == "IV") {
      p <- coef(resp, model)
      ranx <- diff(resp$range)
      minx <- min(resp$range)
      opt <- (p[3] - p[1])/p[2]/2
      opt <- as.numeric(ranx * opt + minx)
      if(opt < minx) opt <- minx; if(opt > max(resp$range)) opt <- max(resp$range)
      top <- as.numeric(predict(resp, newdata = opt, M = M, model = model))
      tmp <- optimize(HOFfun, resp$range, resp = resp, model = model, maximum = FALSE)
      mini <- as.numeric(tmp$objective)
      pess <- as.numeric(tmp$minimum)
  }

  if (model == "V") {
      p <- coef(resp, model)
      if (p[2] * p[4] >= 0) {
          tmp <- optimize(HOFfun3, resp$range, y=0, resp = resp, maximum = TRUE)
          opt <- tmp$maximum
          top <- tmp$objective
          pess <- min(predict(resp, x, M = M, model = model))
          if (top < 16 * .Machine$double.eps) {
              tmp <- seq(resp$range[1], resp$range[2], len = 31)
              ytmp <- predict(resp, newdata = tmp, model = "V")
              tmp <- tmp[ytmp > 16 * .Machine$double.eps]
              tmp <- optimize(HOFfun, range(tmp), resp = resp, model=model, maximum = TRUE)
              opt <- as.numeric(tmp$maximum)
              top <- as.numeric(tmp$objective)
              if (tmp$obj <= 0) opt <- top <- NA
          }
      }
      tmp <- optimize(HOFfun, resp$range, resp = resp, model = model, maximum = FALSE)
      mini <- tmp$objective
      pess <- tmp$minimum
  }
  
#   if (model == "VI") {
#       infl <- c(FALSE, diff(diff(pred)>0)!=0)
#       if(sum(infl) == 2) {# one optimum outside range
#         whichopt <- which.max(pred[which(infl)])
#         if(whichopt == 1) {# second optimum on the right side of the gradient
#           opt <- c(opt1 = x[which(infl)[1]], opt2 = resp$range[2]- +diff(resp$range)/1000)
#           top <- c(top1 = pred[which(infl)[1]], top2 = pred[length(pred)])
#           pess <- x[which(infl)[2]]
#           mini <- pred[which(infl)[2]]
#         } else {
#           opt <- c(opt1 = resp$range[1]- -diff(resp$range)/1000, opt2 = x[which(infl)[2]])
#           top <- c(top1 = pred[which(infl)[2]], top2 = pred[1])
#           pess <- x[which(infl)[1]]
#           mini <- pred[which(infl)[1]]
#         }
#       } else {# both optima inside gradient
#       top <- c(top1 = pred[which(infl)[1]], top2 = pred[which(infl)[3]])
#       pess <- x[which(infl)[2]]
#       opt <- c(opt1 = x[which(infl)[1]], opt2 = x[which(infl)[3]])
#       mini <- pred[which(infl)[2]]
#       }
#       new <- seq(resp$range[1], pess, length.out = 5000)
#       pm <- predict.HOF(resp, newdata = new, model = model)
#       expect1 <- sum(pm * new)/sum(pm)
#       new <- seq(pess, resp$range[2], length.out = 5000)
#       pm <- predict.HOF(resp, newdata = new, model = model)
#       expect2 <- sum(pm * new)/sum(pm) 
#       expect <- c(expect1, expect2)
#    }
# 
#   
  if (model %in% c('VI', 'VII')) {
    max1 <- optimize(HOFfun, resp$range, resp = resp, model=model, maximum = TRUE)
    min.left <- optimize(HOFfun, c(resp$range[1], max1$maximum), model=model, resp = resp, maximum = FALSE)
    min.right <- optimize(HOFfun, c(max1$maximum, resp$range[2]), model=model, resp = resp, maximum = FALSE)
    second.opt <- c('left','right')[which.max(c(diff(c(resp$range[1], min.left$minimum)), abs(resp$range[2] - min.right$minimum)))]
    if(second.opt == 'left') {
      max2 <- optimize(HOFfun, lower=resp$range[1], upper=min.left$minimum, resp = resp, model=model, maximum = TRUE)
      opt <- c(opt1 = max2$maximum, opt2 = max1$maximum)
      top <- c(top1 = max2$objective, top2 = max1$objective)
    }
    if(second.opt == 'right') {
      max2 <- optimize(HOFfun, lower=min.right$minimum, upper=resp$range[2], resp = resp, model=model, maximum = TRUE)
      opt <- c(opt1 = max1$maximum, opt2 = max2$maximum)
      top <- c(top1 = max1$objective, top2 = max2$objective)
    }
    min <- optimize(HOFfun, lower=opt['opt1'], upper=opt['opt2'], resp = resp, model=model, maximum = FALSE)
    mini <- min$objective
    pess <- min$minimum
    new <- seq(resp$range[1], pess, length.out = 10000)
    pm <- predict.HOF(resp, newdata = new, model = model)
    expect1 <- sum(pm * new)/sum(pm)
    new <- seq(pess, resp$range[2], length.out = 10000)
    pm <- predict.HOF(resp, newdata = new, model = model)
    expect2 <- sum(pm * new)/sum(pm) 
    expect <- c(expect1, expect2)
  }
  if(model %in% c('I', 'II', 'III', 'IV', 'V')) {      
		expect <- sum(pred * x)/sum(pred)
	}

  list(opt = opt, top = top*M, pess = pess, mini = mini*M, expect = expect)
}
