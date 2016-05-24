SOPIE <-
function (data = NULL, h = 1, to = 1, alpha = 0.05, g = 20, r = 10, m = 1, grid = 512)
{
   #Checks validity of parameters entered
   if (is.null(data)) 
       stop("\n", "A data vector is needed")

   if ((trunc(g) != g) || (trunc(r) != r) || (trunc(m) != m) || (trunc(grid) != grid) || is.na(r) || is.na(g))
        stop("\n", "g, r, m and grid must be integers")

   #Checks whether the values of g and r is feasible ito the length of the data
   if (length(data) < (g*r)) 
   {
       g<-1 
       r<-1
   }

   cl<-match.call()

   #Estimate the smoothing parameter and calculate the minimum points to be used in the KDE
   sp <- findh(data, h, to)
   min_points <- circ.kernel(data, sp, to = 1, grid, m)

   #Estimate the left and right end point of the interval of uniformity
   a_hat <- a.estimate(data, to, min_points$minimum, alpha, g, r)
   b_hat <- b.estimate(data, to, min_points$minimum, alpha, g, r)

   #Combine results in different data structures for output
   table <- matrix(c(a_hat$summary, b_hat$summary), ncol = 5, byrow = T, dimnames = list(c("a-hat", "b-hat"), c("Cramer von Mises", "Kolmogorov-Smirnoff",
                   "Anderson-Darling", "Rayleigh","MEDIAN")))

   a_hat$General$call <- cl
   a_hat$General$kernel_bandwidth <- sp
   a_hat$General$kernel_grid <- grid

   result<-list(Summary = table, General = a_hat$General)

   # Construct graph of results
   graphics.off()
   hist(data, nclass = 100, freq = FALSE, col = "grey", main = "Histogram, kernel density estimator and est. off-pulse interval", xlab = "", ylab = "")
   height <- 0.5 * (max(hist(data, nclass = 20, plot = FALSE)$density)  + min(hist(data, nclass = 20, plot = FALSE)$density))
   lines(c((table[1, 5] + 0.05), table[1,5], table[1,5], table[1, 5]  +0.05), c(height, height, -0.015, -0.015), col = "blue", lwd = 2)
   lines(c((table[2, 5] - 0.05), table[2,5], table[2,5], table[2,5] - 0.05), c(height, height, -0.015, -0.015), col = "blue", lwd = 2)
   lines(min_points$x, min_points$y / (sum(min_points$y) * (1 / grid)), col = "red", lwd = 2, lty = 1)

   result
}
