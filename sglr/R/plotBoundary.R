if(getRversion() >= "2.15.1")  {
    utils::globalVariables(c("lower", "upper", "label", "y"))
}

plotBoundary <- function(b1, b0, p, glrTables=NULL, tol=1e-7, legend=FALSE, textXOffset=2, textYSkip=2) {
  boundary <- computeBoundary(b1=b1, b0=b0, p=p, glrTables=glrTables)
  estimate <- boundary$estimate
  N <- length(boundary$upper)
  d <- data.frame(x=1:N, lower=boundary$lower, upper=boundary$upper)

  plot <- ggplot()
  plot <- plot + geom_step(aes(x=x, y=lower), data=d, colour="blue") +
      geom_step(aes(x=x, y=upper), data=d, colour = 'red')
  plot <- plot + theme(legend.position = "none")
  plot <- plot + scale_x_continuous('Total No. of AEs') +
      scale_y_continuous('No. of Vaccine AEs')
  x <- floor(sum(range(d$upper, na.rm=TRUE)/2))
  y1 <- d$upper[x] + 2
  plot <- plot + geom_text(aes(x, y, label=label),
                           data=data.frame(x=floor(x/2), y=y1, label="Reject H_0"),
                           color="red", hjust=0.0)
  y2 <- d$lower[x] - 2
  plot <- plot + geom_text(aes(x,y,label=label),
                           data=data.frame(x=ceiling((x + length(d$upper))/2), y=y2, label="Accept H_0"),
                           color="blue", hjust=0.0)

  if (legend) {
      leg1 <- paste("Hypothesis: (p[0] *\",\" * p[1]) * \"=\" * (", p[1], "* \",\" * ", p[2], ")",
                    sep="")

      leg2 <- paste("SGLR *\" \" * Boundaries: (b[0] *\",\" * b[1]) * \"=\" * (", b0, "*\",\" * ", b1, ")",
                    sep="")

      leg3 <- paste("Type *\" \" * I * \", \" * II * \" \" * Errors: (alpha *\", \"* beta) *\"=\" * (",
                    round(estimate[1], 5), "*\", \" * ", round(estimate[2], 5), ")", sep="")

      leg4 <- paste("Max.* \" # of AEs N\"* \"=\" *", N, sep="")

      top <- max(d$upper, na.rm=TRUE)
      plot <- plot + geom_text(aes(x,y,label=label),
                               data=data.frame(x=textXOffset, y=top, label=leg1),
                               parse=TRUE, hjust=0.0)
      plot <- plot + geom_text(aes(x,y,label=label),
                               data=data.frame(x=textXOffset, y=top-textYSkip, label=leg2),
                               parse=TRUE, hjust=0.0)
      plot <- plot + geom_text(aes(x,y,label=label),
                               data=data.frame(x=textXOffset, y=top-2*textYSkip, label=leg3),
                               parse=TRUE, hjust=0.0)
      plot <- plot + geom_text(aes(x,y,label=label),
                               data=data.frame(x=textXOffset, y=top-3*textYSkip, label=leg4),
                               parse=TRUE, hjust=0.0)
  }
  plot
}

