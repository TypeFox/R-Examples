#bstick <- function(n, ...) {
#   UseMethod("bstick")
#}

#bstick.default <- function(n, tot.var=1, ...) rev(cumsum(tot.var/n:1)/n)

bstick.chclust <- function(n, ng=10, plot=TRUE, ...) {
   if (n$method != "coniss")
      stop("bstick cannot display conslink results")
   disp <- rev(n$height)
   tot.disp <- disp[1]
   disp <- abs(diff(disp))
   nobj <- length(n$height)
   bs <- bstick(nobj, tot.disp)
   yR <- range(disp[1:(ng-1)], bs[1:(ng-1)])
   if (plot) {
      plot(2:ng, disp[1:(ng-1)], type="o", ylim=yR, ylab="Sum of Squares", xlab = "Number of groups")
      lines(2:ng, bs[1:(ng-1)], col="red", type="o")
   }
   invisible(data.frame(nGroups = 2:(ng), dispersion=disp[1:(ng-1)], bstick = bs[1:(ng-1)]))
}

