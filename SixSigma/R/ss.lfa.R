# Allow compatibility with previous R versions
if(getRversion() >= '2.15.1') utils::globalVariables(c("..density..", "value"))

#' Loss Function Analysis
#' 
#' This function performs a Quality Loss Function Analysis, based in the Taguchi
#' Loss Function for "Nominal-the-Best" characteristics. 
#' 
#' \code{lfa.output} can take the values "text", "plot" or "both".
#' 
#' @param lfa.data Data frame with the sample to get the average loss.
#' @param lfa.ctq Name of the field in the data frame containing the data.
#' @param lfa.Delta Tolerance of the process.
#' @param lfa.Y0 Target of the process (see note).
#' @param lfa.L0 Cost of poor quality at tolerance limit.
#' @param lfa.size Size of the production, batch, etc. to calculate the total loss in a group 
#' (span, batch, period, ...)
#' @param lfa.output Type of output (see details).
#' @param lfa.sub Subtitle for the graphic output.
#' 
#' @return
#'    \item{lfa.k }{Constant k for the loss function}
#'   \item{lfa,lf }{Expression with the loss function}
#'   \item{lfa.MSD}{Mean Squared Differences from the target}
#'   \item{lfa.avLoss}{Average Loss per unit of the process}
#'   \item{lfa.Loss}{Total Loss of the process (if a size is provided)}
#' 
#' @references 
#' Taguchi G, Chowdhury S,Wu Y (2005) \emph{Taguchi's quality engineering handbook}. John
#' Wiley\cr
#' 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.\cr
#' 
#' @note 
#' For smaller-the-better characteristics, the target should be zero (\code{lfa.Y0 = 0}). 
#' For larger-the-better characteristics, the target should be infinity (\code{lfa.Y0 = Inf}).
#' 
#' @seealso 
#' \code{\link{ss.lf}}, \code{\link{ss.data.bolts}}.
#' 
#' @author EL Cano
#' 
#' @examples
#' ss.lfa(ss.data.bolts, "diameter", 0.5, 10, 0.001, 
#' 		lfa.sub = "10 mm. Bolts Project", 
#' 		lfa.size = 100000, lfa.output = "both")
#' @export
#' @keywords loss function Taguchi
ss.lfa <- function(lfa.data, lfa.ctq, lfa.Delta, lfa.Y0, lfa.L0, 
  lfa.size = NA, lfa.output = "both", lfa.sub = "Six Sigma Project"){
  if (missing(lfa.Delta) || !is.numeric(lfa.Delta)){
    stop("Please provide a valid tolerance value (lfa.Delta).")
  }
  if (lfa.Y0 == 0){
    lfa.k <- lfa.L0/(lfa.Delta^2)
    lfa.lf <- bquote(bold(L==.(lfa.k)%.%Y^2))
    lfa.MSD <- with(lfa.data,
      sum((get(lfa.ctq))^2) / length(get(lfa.ctq)))
  } else if (lfa.Y0 == Inf){
    lfa.k <- lfa.L0*(lfa.Delta^2)
    lfa.lf <- bquote(bold(L==.(lfa.k)%.%(1/Y^2)))
    lfa.MSD <- with(lfa.data,
      sum((1/get(lfa.ctq))^2) / length(get(lfa.ctq)))
    
  } else {
    lfa.k <- lfa.L0/(lfa.Delta^2)
    lfa.lf <- bquote(bold(L==.(lfa.k)%.%(Y-.(lfa.Y0))^2))
    lfa.MSD <- with(lfa.data,
      sum((get(lfa.ctq) - lfa.Y0)^2) / length(get(lfa.ctq)))
  } 
  lfa.avLoss <- lfa.k * lfa.MSD
  
  
  if (is.numeric(lfa.size)){
    lfa.Loss <- lfa.size * lfa.avLoss
  }
  else{
    lfa.Loss <- NA
  }
  
  if (lfa.output %in% c("both", "plot")){
    .ss.prepCanvas(main = "Loss Function Analysis", sub = lfa.sub)
    vp.function <- grid::viewport(name = "Taguchi", 
      layout = grid::grid.layout(2, 2,
        heights = c(0.9, 0.1), 
        widths = c(0.8, 0.2)))
    grid::pushViewport(vp.function)
    vp.plot <- grid::viewport(name = "plot", 
      layout.pos.row = 1, 
      layout.pos.col = 1)
    #plot
    grid::pushViewport(vp.plot)
    ggdata <- reshape2::melt(with(lfa.data, get(lfa.ctq)))
    qqp <- ggplot(ggdata, aes(x = value))
    qqp <- qqp + stat_function(fun = function(x) {
          if (lfa.Y0 == 0){
            eval(lfa.k)*(x)^2
          } else if (lfa.Y0 == Inf){
            eval(lfa.k)*(1/x^2)
          } else {
            eval(lfa.k)*(x-eval(lfa.Y0))^2
          }
        },
        size = 1.2) + 
      ylab("Cost of Poor Quality") + 
      xlab("Observed Value")
    if (lfa.Y0 == Inf) qqp <- qqp + scale_y_continuous(limits = c(0, lfa.L0*10))
    qqp <- qqp + geom_vline(xintercept = eval(lfa.Y0), 
      linetype = 3, size = 1.1)
    qqp <- qqp + geom_hline(yintercept = 0, 
      size = 1)
    
    if (lfa.Y0 != 0){
      xpos <- ifelse(lfa.Y0 == Inf, 0 + lfa.Delta, lfa.Y0 - lfa.Delta)
      qqp <- qqp + geom_vline(xintercept = xpos, 
        linetype = 2)
      qqp <- qqp + annotate(geom = "text", 
        x = xpos, 
        y = lfa.avLoss, 
        label = "LSL", 
        hjust = -0.1)
    }
    if (lfa.Y0 != Inf){
      qqp <- qqp + geom_vline(xintercept = eval(lfa.Y0 + lfa.Delta), 
        linetype = 2)
      qqp <- qqp + annotate(geom = "text", 
        x = lfa.Y0 + lfa.Delta, 
        y = lfa.avLoss, 
        label = "USL", 
        hjust = 1.1)
      qqp <- qqp + annotate(geom = "text", 
        x = lfa.Y0, 
        y = ss.lf(lfa.Y0 - lfa.Delta, 
          lfa.Delta, 
          lfa.Y0, 
          lfa.L0), 
        label = "T", 
        hjust = 1.1)
    }
    if (lfa.Y0 == Inf){
      qqp <- qqp + geom_vline(xintercept = 0, 
        linetype = 1, size = 1) 
    }
    print(qqp, newpage = FALSE)
    grid::popViewport()
    
    #function
    vp.fun <- grid::viewport(name = "fun", 
      layout.pos.row = 2, 
      layout.pos.col = 1:2)
    grid::pushViewport(vp.fun)
    grid::grid.rect(width = 0.95, 
      gp = grid::gpar(lty = 0, fill = "#DDDDDD"))
    grid::grid.text(lfa.lf)
    grid::popViewport()
    
    #data
    vp.data <- grid::viewport(name = "data", 
      layout.pos.row = 1:2, 
      layout.pos.col = 2)
    grid::pushViewport(vp.data)
    vp.data.input <- grid::viewport(name="input", 
      layout.pos.row = 2, 
      layout.pos.col = 1)
    grid::pushViewport(vp.data.input)
    grid::grid.rect(y = 0.95, 
      width = 0.99, 
      just = "top", 
      height=0.7)
    my.margin <- 0.9
    grid::grid.text(expression(bold(Data)), 
      y = my.margin, 
      just = "top")
    grid::grid.text(paste("CTQ:", eval(lfa.ctq)),
      y = unit(my.margin, "npc") - unit(2, "lines"),
      just = "top",
      gp = grid::gpar(cex = 0.8))
    grid::grid.text(bquote(Y[0]==.(lfa.Y0)),
      y = unit(my.margin, "npc") - unit(3, "lines"),
      just = "top")
    grid::grid.text(bquote(Delta==.(lfa.Delta)),
      y = unit(my.margin, "npc") - unit(4, "lines"),
      just = "top")
    grid::grid.text(bquote(L[0]==.(lfa.L0)),
      y = unit(my.margin,"npc") - unit(5, "lines"),
      just = "top")
    if (is.numeric(lfa.size)){
      size.exists = 1
      grid::grid.text(bquote(Size==.(lfa.size)),
        y = unit(my.margin, "npc") - unit(6, "lines"),
        just = "top")
    }
    grid::grid.lines(y = unit(my.margin, "npc") - unit(8, "lines"))
    
    grid::grid.text(bquote(Mean==.(round(with(lfa.data, mean(get(lfa.ctq))), 
            digits = 4))),
      y = unit(my.margin,"npc") - unit(10, "lines"),
      just = "top")
    grid::grid.text(bquote(k==.(lfa.k)),
      y = unit(my.margin,"npc") - unit(11, "lines"),
      just = "top")
    grid::grid.text(bquote(MSD==.(round(lfa.MSD, digits = 4))),
      y = unit(my.margin, "npc") - unit(12, "lines"),
      just = "top")
    grid::grid.text(bquote(Av.Loss==.(round(lfa.avLoss, digits = 4))),
      y = unit(my.margin, "npc") - unit(13, "lines"),
      just = "top")
    if (is.numeric(lfa.size)){
      grid::grid.text(bquote(Loss==.(round(lfa.Loss, digits = 4))),
        y = unit(my.margin, "npc") - unit(14, "lines"),
        just = "top")
    }
  }
  if (lfa.output %in% c("both", "text")){
    #pintar valores en output
    return(list(lfa.k = lfa.k, 
        lfa.lf = as.expression(lfa.lf), 
        lfa.MSD = lfa.MSD, 
        lfa.avLoss = lfa.avLoss, 
        lfa.Loss = lfa.Loss))
  }
  else{
    invisible(list(lfa.k = lfa.k, 
        lfa.lf = as.expression(lfa.lf), 
        lfa.MSD = lfa.MSD, 
        lfa.avLoss = lfa.L0, 
        lfa.Loss = lfa.Loss))
  }
}



#' Evaluates the Loss Function for a process.
#' 
#' The quality loss function is one of the tools of the Six Sigma methodology.
#' The function assigns a cost to an observed value, that is larger as far as it
#' is from the target. 
#' 
#' @param lfa.Y1 The observed value of the CTQ (critical to quality) characteristic
#' that will be evaluated.
#' @param lfa.Delta The tolerance for the CTQ.
#' @param lfa.Y0 The target for the CTQ.
#' @param lfa.L0 The cost of poor quality when the characteristic is \eqn{Y_0 + \Delta}.
#' 
#' @return
#' \item{ss.lf}{A number with the evaluated function at \eqn{Y_1}} 
#' 
#' @references 
#' Taguchi G, Chowdhury S,Wu Y (2005) \emph{Taguchi's quality engineering handbook}. John
#' Wiley
#' 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' @seealso \code{\link{ss.lfa}}
#' @author EL Cano
#' 
#' @examples 
#' #Example bolts: evaluate LF at 10.5 if Target=10, Tolerance=0.5, L_0=0.001
#' ss.lf(10.5, 0.5, 10, 0.001)
#' @export
#' @keywords loss function Taguchi
ss.lf <- function(lfa.Y1, lfa.Delta, lfa.Y0, lfa.L0) {
  if (lfa.Delta <= 0){
    stop("The tolerance of the process must be greater than 0")
  }
  if (lfa.L0 <= 0){
    warning("The Cost of poor quality at tolerance limit should be greater than 0")
  }
  lfa.k <- lfa.L0/lfa.Delta
  return(lfa.k*(lfa.Y1-lfa.Y0)^2)
  
}
