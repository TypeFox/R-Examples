#' Find optimum reparameterization between two images
#'
#' Finds the optimal warping function between two images using the elastic
#' framework
#'
#' @param It template image matrix
#' @param Im test image matrix
#' @param gam initial warping array
#' @param b basis matrix
#' @param stepsize gradient stepsize (default=1e-5)
#' @param itermax maximum number of iterations (default=1000)
#' @param lmark use landmarks (default=FALSE)
#' @return Returns a list containing \item{gamnew}{final warping}
#' \item{Inew}{aligned image}
#' \item{H}{energy}
#' \item{stepsize}{final stepsize}
#' @keywords image alignment
#' @references Q. Xie, S. Kurtek, E. Klassen, G. E. Christensen and A. Srivastava. Metric-based pairwise and multiple image registration. IEEE European Conference on Computer Vision (ECCV), September, 2014
#' @export
reparam_image <- function(It, Im, gam, b, stepsize=1e-5, itermax=1000, lmark=FALSE){
    m = dim(It)[1]
    n = dim(It)[2]
    gamid = makediffeoid(m,n)

    # Main Loop-----------------------------------------------------------------
    H = rep(0,itermax+1)
    Iold = apply_gam_to_imag(Im, gam)
    Iold = Iold - min(Iold)
    Iold = Iold / max(Iold)
    qt = imag_to_q(It)
    qm = imag_to_q(Iold)

    gamold = gam
    gamnew = gamold
    Inew = Iold
    iter = 1
    H[iter] = comp_energy(qt, qm)
    cat(sprintf("Iteration %d, energy %f\n",iter-1,H[iter]))

    gamupdate = update_gam(qt,qm,b)
    cutoff = 1e-3
    for (iter in 2:(itermax+1)) {
        gaminc = gamid + stepsize*gamupdate
        G = check_crossing(gamnew)
        if (!G){
            cat("Possible Crossing!\n")
            gamnew = gamold
            stepsize = 0.67*stepsize
            H[iter] = H[iter+1]
            next
        } else {
            gamnew = apply_gam_to_gam(gamnew, gaminc)
        }

        Inew = apply_gam_to_imag(Im, gamnew)
        Inew = Inew - min(Inew)
        Inew = Inew / max(Inew)
        qm = imag_to_q(Inew)
        H[iter] = comp_energy(qt,qm)
        cat(sprintf("Iteration %d, energy %f\n",iter-1,H[iter]))

        if (iter>4){
            hstop = 1.
            for (i in 1:4) {
                hstop = hstop*(H[iter]>=H[iter-1])
            }
            if (hstop!=0){
                cat("Warning: energy constantly increasing")
            }
        }

        if ((iter>4) && (H[iter]>=H[iter-1]) && (H[iter-1]>=H[iter-2]) && (H[iter-2] >= H[iter-3])) {
            cat("Warning: energy is not changing")
        }

        if (((iter>1) && (H[iter]>H[iter-1])) || ((iter>3) && ((H[iter-1]<=H[iter-2]) && (H[iter-2]>H[iter-3])))) {
            stepsize = 0.9*stepsize
            gamnew = gamold
            H[iter] = H[iter-1]
            next
        }

        gamold = gamnew
        gamupdate = update_gam(qt,qm,b)
    }

    H = H[1:iter]

    return(list(gamnew=gamnew,Inew=Inew,H=H,stepsize=stepsize))
}
