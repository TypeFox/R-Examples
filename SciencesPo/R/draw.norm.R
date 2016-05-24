#' @encoding UTF-8
#' @title Shades Normal Distribuion
#'
#' @description Produces a plot of a normal density distribution with shaded areas.
#'
#' @param below sets a lower endpoint.
#' @param above sets an upper endpoint.
#' @param pcts the
#' @param mu the mean.
#' @param sigma standard deviations.
#' @param numpts the number os points/observations to drawn upon.
#' @param color the color of the area.
#' @param dens the density of the color.
#' @param justabove just plots the upper tail.
#' @param justbelow just plots the lower tail.
#' @param lines to draw lines.
#' @param between plots between specified points.
#' @param outside alternative "outside" area.
#'
#' @return A plot with a normal distribution density with shaded areas
#'
#'@examples
#'
#' draw.norm()
#' draw.norm(below=-1.5)
#' draw.norm(below=-1.5,justbelow=TRUE)
#' draw.norm(above=1.5, justabove=TRUE)
#' draw.norm(below=-1.5,above=1.5)
#' draw.norm(between=c(-4,0),color="black")
#' draw.norm(between=c(0,4),color="black")
#' draw.norm(between=c(-1,+1),color="darkgray")
#' title("P[-1 < z < 1] = 68%")
#' draw.norm(between=c(-2,+2),color="darkgray")
#' title("P[-2 < z < 2] = 95%")
#' draw.norm(between=c(-3,+3),color="darkgray")
#' title("P[-3 < z < 3] = 99.7%")
#' draw.norm(between = c(-1.75, 0, 2, 0.5, -1))  ## Plots between specified points
#' draw.norm(below=50,justbelow=TRUE,color="black",mu=47.3,sigma=9.3)
#'
#' ## Can plot one and then another on top of it using lines = TRUE
#' draw.norm(mu=2, sigma=10, outside=c(-3, 12), dens=15)
#' draw.norm(mu=2, sigma=15, between=c(-3, 12),lines=TRUE, col="blue",dens=15)
#' ## Example: Plotting a Hypothesis Test for the mean
#' ## Truth:      mu.true  = 8
#' ## Hypothesis: mu.ho    = 6
#' ## Generate Data Under Truth
#' mu.true = 5 ## Alternative Mean
#' mu.ho   = 6
#' sig     = 8
#' N       = 250 ## Sample Size
#'
#' std.err = sig/sqrt(N)
#' crits = qnorm(c(0.025,0.975),mean=mu.ho, sd = std.err)
#' draw.norm(outside = crits, mu = mu.ho, sigma = std.err,dens=15)
#' draw.norm(between = crits, mu = mu.true, sigma = std.err, lines=TRUE, color="green",dens=15)
#'
#'
#'@export
`draw.norm` <- function(below=NULL, above=NULL, pcts = c(0.025,0.975), mu=0, sigma=1, numpts = 500, color = "gray", dens = 40, justabove= FALSE, justbelow = FALSE, lines=FALSE, between=NULL, outside=NULL) {

  if(is.null(between)){
    below = ifelse(is.null(below), stats::qnorm(pcts[1],mu,sigma), below)
    above = ifelse(is.null(above), stats::qnorm(pcts[2],mu,sigma), above)
  }

  if(is.null(outside)==FALSE){
    below = min(outside)
    above = max(outside)
  }
  lowlim = mu - 4*sigma
  uplim  = mu + 4*sigma

  x.grid = seq(lowlim,uplim, length= numpts)
  dens.all = stats::dnorm(x.grid,mean=mu, sd = sigma)
  if(lines==FALSE){
    graphics::plot(x.grid, dens.all, type="l", xlab="X", ylab="Density")
  }
  if(lines==TRUE){
    graphics::lines(x.grid,dens.all)
  }

  if(justabove==FALSE){
    x.below    = x.grid[x.grid<below]
    dens.below = dens.all[x.grid<below]
    graphics::polygon(c(x.below,rev(x.below)),c(rep(0,length(x.below)),rev(dens.below)),col=color,density=dens)
  }
  if(justbelow==FALSE){
    x.above    = x.grid[x.grid>above]
    dens.above = dens.all[x.grid>above]
    graphics::polygon(c(x.above,rev(x.above)),c(rep(0,length(x.above)),rev(dens.above)),col=color,density=dens)
  }

  if(is.null(between)==FALSE){
    from = min(between)
    to   = max(between)

    x.between    = x.grid[x.grid>from&x.grid<to]
    dens.between = dens.all[x.grid>from&x.grid<to]
    graphics::polygon(c(x.between,rev(x.between)),c(rep(0,length(x.between)),rev(dens.between)),col=color,density=dens)
  }
}
NULL
