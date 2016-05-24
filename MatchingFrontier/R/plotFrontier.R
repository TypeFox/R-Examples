plotFrontier <-
function(frontier.object,
         xlab = 'Number of Observations Pruned',
         ylab = frontier.object$metric,
         main = 'Frontier Plot',
         ...){

    plot(frontier.object$frontier$Xs, frontier.object$frontier$Ys,
         xlab = xlab,
         ylab = ylab,
         main = main,
         ...)
}
