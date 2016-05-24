
predict.moult <- function(object, newdata = NULL, predict.type = "prob", intervals = 0.1, ...)
{ 		# --- predicted proportions not in moult, in moult, and completed

#   if (!is.element(predict.type, c("prob", "response", "start")))
#     stop("Type of prediction not available") 

   if (missing(newdata)) 
     
    { return(object$fitted.values) }

   switch(predict.type,

          ## predict moult score (fitted values)
     response = {return(object$fitted.values)},

          ## predict probability that in certain moult stage
     prob = { p <- c(object$coefficients, recursive = TRUE)
              p <- as.vector(p)

              data.type <- object$type

              p.dur <- unlist(object$coefficients$duration)
              p.mean <- unlist(object$coefficients$mean)
              p.sd <- unlist(object$coefficients$sd)
      
              day <- newdata[, 1]              
              mm <- model.matrix(object$terms$mean, newdata)
              md <- model.matrix(object$terms$duration, newdata)

              if (!is.null(object$terms$sd))
                  { msd <- model.matrix(object$terms$sd, newdata)}
              else
                  { msd <- model.matrix(~ 1, newdata) }

              no.moult <- 1 - pnorm(day, mean = mm %*% p.mean, sd = msd %*% p.sd)
              
              in.moult <- pnorm(day, mean = mm %*% p.mean, sd = msd %*% p.sd) - 
                  pnorm(day - p[1], mean = mm %*% p.mean, sd = msd %*% p.sd ) 
              compl.moult <- pnorm(day - p[1], mean = mm %*% p.mean, sd = msd %*% p.sd )

              NIC <- cbind(no.moult, in.moult, compl.moult)
              colnames(NIC) <- c("pre-moult", "in moult", "post-moult")

              ints <- seq(0, 1, by = intervals)
           
              IM <- matrix(0, ncol = length(ints) - 1, nrow = length(day))
              rownames(IM) <- day
              colnames(IM) <- paste(ints[1:length(ints) - 1], "-", ints[2:length(ints)])
              for (i in 2:length(ints))
                  { IM[,i-1] <- pnorm(day - md %*% p.dur * ints[i-1], mean = mm %*% p.mean, sd = msd %*% p.sd) - 
                        pnorm(day - md %*% p.dur * ints[i], mean = mm %*% p.mean, sd = msd %*% p.sd) 
                }

              switch(data.type, 
                     { M <- round(cbind(NIC[, 1], IM, NIC[,3]), digits = 3) },  # --- type 1
                     { M <- round(cbind(NIC[, 1], IM, NIC[,3]), digits = 3) }, # --- type 2
                     { M <- round(cbind(NIC[, 1], IM / in.moult , NIC[, 3]), digits = 3) },  # --- type 3
                     { M <- round(cbind(NIC[, 1], IM / ( 1 - no.moult), compl.moult / (1 - no.moult)), digits = 3) },  # --- type 4
                     { M <- round(cbind(no.moult / (1 - compl.moult), IM / ( 1 - compl.moult), compl.moult), digits = 3) },     # --- type 5                       
                     print("predict.moult: not a valid type") )                   # ---- end of switch

              colnames(M) <- c("0", colnames(IM), "1")
              out <- list("M" = M)
              return(out)
          }, 

          ## predict.type = "start"
          start = { VC <- object$vcov

                    p.dur <- unlist(object$coefficients$duration)
                    p.mean <- unlist(object$coefficients$mean)
           
                    l.dur <- length(p.dur)
                    l.mean <- length(p.mean)

                    VC <- VC[-(1:l.dur), -(1:l.dur)]
                    VC <- VC[1:l.mean, 1:l.mean]
           
                    mm <- model.matrix(object$terms$mean, newdata)
           
                    pred.start <- mm %*% p.mean
                    pred.var <- mm %*% VC %*% t(mm) 
                    pred.se <- sqrt(diag(pred.var))

                    out <- data.frame(mean.start = pred.start, SE = pred.se) 
                    return(out)
                },

          ## predict duration
    duration = { VC <- object$vcov

                 p.dur <- unlist(object$coefficients$duration)
                 
                 l.dur <- length(p.dur)
                  
                 VC <- VC[1:l.dur, 1:l.dur]
           
                 mm <- model.matrix(object$terms$duration, newdata)
           
                 pred.duration <- mm %*% p.dur
                    pred.var <- mm %*% VC %*% t(mm) 
                    pred.se <- sqrt(diag(pred.var))

                    out <- data.frame(duration = pred.duration, SE = pred.se) 
                    return(out)
                },

          ## predict both start and duration 
    both = { VC <- object$vcov

             p.dur <- unlist(object$coefficients$duration)
             p.mean <- unlist(object$coefficients$mean)
             p <- c(p.dur, p.mean)
           
             l.dur <- length(p.dur)
             l.mean <- length(p.mean)
             l <- l.dur + l.mean

             VC <- VC[1:l, 1:l]
           
             mm1 <- model.matrix(object$terms$mean, newdata)
             mm2 <- model.matrix(object$terms$duration, newdata)
             mm <- as.matrix(bdiag(mm1, mm2))
             
             pred.est  <- mm %*% p
             pred.var <- mm %*% VC %*% t(mm) 
             pred.se <- sqrt(diag(pred.var))

             nrows <- dim(newdata)[1]

             out <- data.frame(duration.est = pred.est[1:nrows], duration.se = pred.se[1:nrows],
                                     start.est = pred.est[-(1:nrows)], start.se = pred.se[-(1:nrows)],
                                     dur.start.cov = diag(pred.var[(1:nrows), -(1:nrows)]))

             # %out <- list(estimates = data.frame(estimate = pred.est, SE = pred.se), vcov.matrix = pred.var) 
             return(out)
         },

         print("predict.type not valid (predict.moult)"))

}



