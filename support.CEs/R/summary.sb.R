summary.sb <- function(object, ...)
{
  x <- object

  operation       <- x$sb$input$operation
  effect          <- x$sb$input$effect
  determinant.c   <- x$sb$input$determinant.c
  interaction     <- x$sb$input$interaction
  messages        <- x$sb$messages
  version         <- x$sb$version

  b.chr           <- x$sb$output$b.chr
  l.chr           <- x$sb$output$l.chr
  c.chr           <- x$sb$output$c.chr
  cinv.chr        <- x$sb$output$cinv.chr
  correlation.chr <- x$sb$output$correlation.chr

  dimnames(b.chr)           <- list(c(1:nrow(b.chr)),
                                    c(1:ncol(b.chr)))
  dimnames(l.chr)           <- list(c(1:nrow(l.chr)),
                                    c(1:ncol(l.chr)))
  dimnames(c.chr)           <- list(c(1:nrow(c.chr)),
                                    c(1:ncol(c.chr)))
  dimnames(cinv.chr)        <- list(c(1:nrow(cinv.chr)),
                                    c(1:ncol(cinv.chr)))
  dimnames(correlation.chr) <- list(c(1:nrow(correlation.chr)),
                                    c(1:ncol(correlation.chr)))

  if (x$sb$input$operation == "construct") {
    treatment   <- x$sb$input$treatment
    generators  <- x$sb$input$generators
    choice.sets <- x$sb$output$choice.sets
    dimnames(treatment)   <- list(c(1:nrow(treatment)),
                                  c(1:ncol(treatment)))
    dimnames(generators)  <- list(c(1:nrow(generators)),
                                  c(1:ncol(generators)))
    dimnames(choice.sets) <- list(c(1:nrow(choice.sets)),
                                  c(1:ncol(choice.sets)))
    choice.sets.check <- NULL
  } else {
    choice.sets.check <- x$sb$input$choice.sets.check
    dimnames(choice.sets.check) <- list(c(1:nrow(choice.sets.check)),
                                        c(1:ncol(choice.sets.check)))
    treatment   <- NULL
    generators  <- NULL
    choice.sets <- NULL
  }

  rtn <- list(
    operation         = operation,
    effect            = effect,
    treatment         = treatment, 
    generators        = generators,
    determinant.c     = determinant.c,
    interaction       = interaction,
    choice.sets.check = choice.sets.check,
    messages          = messages,
    version           = version,
    choice.sets       = choice.sets,
    b.chr             = b.chr,
    l.chr             = l.chr,
    c.chr             = c.chr,
    cinv.chr          = cinv.chr,
    correlation.chr   = correlation.chr)

  class(rtn) <- c("summary.sb")
  return(rtn)
}
