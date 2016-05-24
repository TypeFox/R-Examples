#' @title Empirical Rule

#' @description An app to investigate how the Empirical Rule applies to symmetric data and skewed data.  The user can select
#' is they want to view a histogram of symmetric data or skewed data.  Vertical bars are also plotted to signify 
#' one, two, and three standard deviations from the mean.  Summary data is output to the console giving the proportion of the histogram that falls within one, two, 
#' and three standard deviations of the mean.  
#' 
#' @rdname EmpRule
#' @usage EmpRule()
#' @return Graphical and numerical output
#' @export
#' @author Rebekah Robinson \email{rebekah_robinson@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) EmpRule()
#' }
EmpRule <- function () 
{
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  rpareto <- function(n,alpha,theta) {#random values for Pareto(alpha,theta) distribution
    theta*((1-runif(n))^(-1/alpha)-1)
  }
  
  dpareto <- function(x,alpha,theta) {  #pdf for Pareto(alpha,theta) distribution
    alpha*theta^alpha/(x+theta)^(alpha+1)
  }
  
  
  
  results = matrix(0, 1, 3, dimnames = list("Proportion", 
                                            c("One SD", "Two SD", "Three SD")))
  manipulate(n = slider(5, 1000, initial = 50, label = "Sample Size n"), 
             type = picker("Symmetric", "Skewed", "Super-Skewy", label = "Target Data Shape"),
             showpop=checkbox(FALSE,"Show Population Density Curve"),
{
  if (type == "Symmetric") {
    tally1 = 0
    tally2 = 0
    tally3 = 0
    mu = 20
    stdev = 2
    data = rnorm(n, mean = mu, sd = stdev)
    xbar = mean(data)
    s = sd(data)
    hist(data, freq = FALSE, col = "lightblue", 
         xlim = c(mu - 5 * stdev, mu + 5 * stdev), ylim = c(0, 0.35),
         main = "Empirical Rule with Target Symmetric")
    abline(v = xbar-s, col = "red")
    abline(v = xbar+s, col = "red")
    points(xbar, 0, pch = 20)
    abline(v = xbar - 2 * s, col = "blue")
    abline(v = xbar + 2 * s, col = "blue")
    abline(v = xbar - 3 * s, col = "green")
    abline(v = xbar + 3 * s, col = "green")
    if (showpop) curve(dnorm(x,mean=mu,sd=stdev),col="red",add=TRUE,n=1001)
    for (i in 1:n) {
      if (data[i] < xbar + s & data[i] > xbar - 
            s) {
        tally1 = tally1 + 1
      }
      if (data[i] < xbar + 2 * s & data[i] > xbar - 
            2 * s) {
        tally2 = tally2 + 1
      }
      if (data[i] < xbar + 3 * s & data[i] > xbar - 
            3 * s) {
        tally3 = tally3 + 1
      }
    }
    results[, 1] = round(tally1/n, 4)
    results[, 2] = round(tally2/n, 4)
    results[, 3] = round(tally3/n, 4)
    print(results)
  }
  if (type == "Skewed") {
    tally1 = 0
    tally2 = 0
    tally3 = 0
    alpha = 1.5
    beta = 5
    data = rbeta(n, alpha, beta)
    mu = 1/(1 + (beta/alpha))
    stdev = sqrt((alpha * beta)/((alpha + beta)^2 * 
                                   (alpha + beta + 1)))
    xbar = mean(data)
    s = sd(data)
    breaks.symm = ifelse(n<200,"Sturges","Scott")
    hist(data, freq = FALSE, col = "lightblue", breaks=breaks.symm,
         xlim = c(mu - 3 * stdev, mu + 5 * stdev), ylim = c(0, 5),
         main = "Empirical Rule with Target Skewed")
    if (showpop) curve(dbeta(x,shape1=alpha,shape2=beta),
                       xlim=c(0,1),col="red",add=TRUE,n=1001)
    abline(v = xbar - s, col = "red")
    abline(v = xbar + s, col = "red")
    points(xbar, 0, pch = 20)
    abline(v = xbar - 2 * s, col = "blue")
    abline(v = xbar + 2 * s, col = "blue")
    abline(v = xbar - 3 * s, col = "green")
    abline(v = xbar + 3 * s, col = "green")
    for (i in 1:n) {
      if (data[i] < xbar + s & data[i] > xbar - 
            s) {
        tally1 = tally1 + 1
      }
      if (data[i] < xbar + 2 * s & data[i] > xbar - 
            2 * s) {
        tally2 = tally2 + 1
      }
      if (data[i] < xbar + 3 * s & data[i] > xbar - 
            3 * s) {
        tally3 = tally3 + 1
      }
    }
    results[, 1] = round(tally1/n, 4)
    results[, 2] = round(tally2/n, 4)
    results[, 3] = round(tally3/n, 4)
    print(results)
  }
  if (type == "Super-Skewy") {
    tally1 = 0
    tally2 = 0
    tally3 = 0
    
    p.alpha = 3  #parameter in 2-parameter pareto(alpha,theta) distribution
    p.theta <- 100 #parameter in pareto
    
    rpareto <- function(n,alpha,theta) {
      theta*((1-runif(n))^(-1/alpha)-1)
    }
    
    data = rpareto(n,p.alpha,p.theta)
    mean.par = p.theta/(p.alpha-1)
    sd.par  = mean.par*sqrt(p.alpha/(p.alpha-2))
    xbar=mean(data)
    s=sd(data)
    xmin = mean.par - 3 * sd.par
    xmax = mean.par+10*sd.par
    ymax = 1.3*p.alpha/p.theta
    hist(data, freq = FALSE, col = "lightblue", breaks="FD",
         xlim = c(xmin,xmax), 
         ylim = c(0, ymax),
         main = "Empirical Rule with Target Super-Skewy")
    if (showpop) curve(dpareto(x,alpha=p.alpha,theta=p.theta),
                       xlim=c(0,xmax),col="red",add=TRUE,n=1001)
    abline(v = xbar - s, col = "red")
    abline(v = xbar + s, col = "red")
    points(xbar, 0, pch = 20)
    abline(v = xbar - 2 * s, col = "blue")
    abline(v = xbar + 2 * s, col = "blue")
    abline(v = xbar - 3 * s, col = "green")
    abline(v = xbar + 3 * s, col = "green")
    for (i in 1:n) {
      if (data[i] < xbar + s & data[i] > xbar - 
            s) {
        tally1 = tally1 + 1
      }
      if (data[i] < xbar + 2 * s & data[i] > xbar - 
            2 * s) {
        tally2 = tally2 + 1
      }
      if (data[i] < xbar + 3 * s & data[i] > xbar - 
            3 * s) {
        tally3 = tally3 + 1
      }
    }
    results[, 1] = round(tally1/n, 4)
    results[, 2] = round(tally2/n, 4)
    results[, 3] = round(tally3/n, 4)
    print(results)
  }
})
}

# if(getRversion() >= "2.15.1")  utils::globalVariables(c("type","showpop"))
