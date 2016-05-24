#' Amplitude Plot VIC and FAM Channels of a Droplet Digital PCR Experiment
#' 
#' This function generates an amplitude plot of two fluorescence channels as
#' found in droplet digital PCR.
#' 
#' Droplet digital PCR experiments consist of three steps (droplet generation,
#' clonal amplification, droplet amplitude analysis). Typically 20000
#' nano-sized droplets are analyzed and separated into amplification-positive
#' and amplification-negative droplets. An example of such system is the
#' Bio-Rad QX100 and QX200 (Pinheiro et al. 2012). Such systems have
#' applications in the detection of rare DNA target copies, the determination
#' of copy number variations (CNV), detection of mutation, or expression
#' analysis of genes or miRNA. Each droplet is analyzed individually using a
#' virtual two-color detection system. The channels are treated separately but
#' finally aligned (e.g., FAM and VIC or FAM and HEX).
#' 
#' @param vic Amplitudes of the VIC channel - object of class
#' \code{\linkS4class{ddpcr}}.
#' @param fam Amplitudes of the FAM channel - object of class
#' \code{\linkS4class{ddpcr}}.
#' @param col_vic Color of the VIC channel.
#' @param col_fam Color of the FAM channel.
#' @param circle If TRUE circles are drawn, if FALSE not. If "numeric",
#' specifies the radius of circles.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @references Pinheiro, L.B., Coleman, V.A., Hindson, C.M., Herrmann, J.,
#' Hindson, B.J., Bhat, S., and Emslie, K.R. (2012). \emph{Evaluation of a
#' droplet digital polymerase chain reaction format for DNA copy number
#' quantification}. Anal. Chem. 84, 1003 - 1011.
#' @keywords hplot
#' @examples
#' 
#' # Generate an amplitude plot for the first fluorescence channel (e.g., FAM)
#' fluos1 <- sim_ddpcr(m = 16, n = 30, times = 100, pos_sums = FALSE, n_exp = 1, 
#'   fluo = list(0.1, 0))
#' 
#' # Generate an amplitude plot for the second fluorescence channel (e.g., VIC)
#' fluos2 <- sim_ddpcr(m = 16, n = 30, times = 100, pos_sums = FALSE, n_exp = 1, 
#'   fluo = list(0.1, 0))
#' 
#' # Plot the amplitudes of both fluorescence channel in an aligned fashion
#' plot_vic_fam(fam = fluos1, vic = fluos2)
#' 
#' # Same as above but different colors
#' plot_vic_fam(fam = fluos1, vic = fluos2, col_vic = "red", col_fam = "yellow")
#' 
#' # Same as above without circles
#' plot_vic_fam(fam = fluos1, vic = fluos2, col_vic = "red", col_fam = "yellow", circle = FALSE)
#' 
#' # Generate two channels in one object and plot them
#' fluos_both <- sim_ddpcr(m = 16, n = 30, times = 100, pos_sums = FALSE, n_exp = 2, 
#'   fluo = list(0.1, 0))
#' plot_vic_fam(extract_dpcr(fluos_both, 1), extract_dpcr(fluos_both, 2))
#' 
#' @export plot_vic_fam
plot_vic_fam <- function(vic, fam, col_vic = "green", col_fam = "blue", circle = TRUE) {
  if (class(vic) == "ddpcr" && class(fam) == "ddpcr") { 
    if (ncol(vic) > 1 && ncol(fam) > 1)
      stop("Both 'vic' and 'fam' must contain only one panel.", call. = TRUE, domain = NA)    
    if (nrow(vic) == 1 && nrow(fam) == 1)
      stop("Both 'vic' and 'fam' cannot contain total number of positive chambers.", call. = TRUE, 
           domain = NA)    
  } else {
    stop("Both 'vic' and 'fam' must have the 'ddpcr' class", call. = TRUE, domain = NA)
  }
  vic_thr <- slot(vic, "threshold")
  fam_thr <- slot(fam, "threshold")
  vic <- as.vector(vic)
  fam <- as.vector(fam)
  
  y_max_lim <- ifelse(max(vic) > max(fam), max(vic), max(fam))
  y_min_lim <- ifelse(min(vic) > min(fam), min(vic), min(fam))
  
  #center of the circle
  y_points <- -y_max_lim * 0.06
  
  #end of the 'tick'
  y_points2 <- y_points * 0.3
  plot(vic, col = col_vic, type = "l", lwd = 2, ylim = c(y_points, y_max_lim), 
       axes = FALSE, xlab = "", ylab = "Number of molecules")
  axis(2)
  lines(fam, col = col_fam, lty = "dashed")
  abline(h = 0)
  
  #translate color to rgb
  col_vic <- col2rgb(col_vic)/255
  col_fam <- col2rgb(col_fam)/255
  
  vic_pos <- findpeaks(vic, threshold = vic_thr)[, 1:2]
  fam_pos <- findpeaks(fam, threshold = fam_thr)[, 1:2]
  circ <- unique(c(vic_pos[, 2], fam_pos[, 2]))
  
  
  #circle drawing
  if(circle == TRUE) 
    radius <- (y_max_lim - y_min_lim)*5
  if(is.numeric(circle)) {
    if(length(circle) == 1) {
      radius <- circle
    } else {
      stop("If 'circle' has type 'numeric', it must have length 1.", call. = TRUE, domain = NA)
    }
  }
  
  for (i in circ)
    lines(c(i, i), c(0, y_points2))
  
  if (circle != FALSE) {
    plot_vf_circ(circ, y_points, radius, "white") 
    plot_vf_circ(vic_pos[, 2], y_points, radius, 
                 rgb(col_vic[1, 1],
                     col_vic[2, 1], 
                     col_vic[3, 1], 
                     vic_pos[, 1]/y_max_lim))
    plot_vf_circ(fam_pos[, 2], y_points, radius, 
                 rgb(col_fam[1, 1],
                     col_fam[2, 1], 
                     col_fam[3, 1], 
                     fam_pos[, 1]/y_max_lim))
    
    #peaks common for both channels
    common <- fam_pos[fam_pos[, 2] %in% vic_pos[, 2], 2]
    if (length(common) != 0) {
      fracs <- sapply(common, function(x) 
        c(vic[x]/(vic[x] + fam[x]), (vic[x] + fam[x])/2))
      common_cols <- rgb(0, fracs[1, ], 1 - fracs[1, ], alpha = fracs[2, ]/y_max_lim)
      plot_vf_circ(common, y_points, radius, common_cols)
    }
  }
}
