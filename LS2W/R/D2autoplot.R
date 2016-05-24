`D2autoplot` <-
function (J, filter.number = 1, family = "DaubExPhase", direction = 3, 
    main = "2-D Autocorrelation Wavelet", OPLENGTH = 10000, scaling = "other", box=TRUE) 
{
    if (J >= 0) 
        stop("J must be a negative integer")
    if (J - round(J) != 0) 
        stop("J must be an integer")
    if (direction != 1 && direction != 2 && direction != 3) 
        stop("Direction can only take the arguments 1, 2, and 3!")
    if (direction == 3) {
        tmp <- D2ACW(J = J, filter.number = filter.number, family = family, 
            OPLENGTH = OPLENGTH)
        nr <- nrow(tmp[[-3 * J]])
        nc <- ncol(tmp[[-3 * J]])
        m <- matrix(0, nrow = nr, ncol = nc)
        m <- tmp[[-3 * J]]
        if (scaling == "Haar") {
            x <- c(1:nr)/2^(-J) - 1
            y <- c(1:nc)/2^(-J) - 1
        }
        else if (scaling == "other") {
            x <- c(1:nr)
            y <- c(1:nc)
        }
        persp(x, y, m, theta = 30, phi = 30, expand =0.5, xlab = "tau_1", ylab = "tau_2", zlab = "ACW coeff", box = box)
        title(main)
    }
    if (direction == 2) {
        tmp <- D2ACW(J = J, filter.number = filter.number, family = family, 
            OPLENGTH = OPLENGTH)
        nr <- nrow(tmp[[-2 * J]])
        nc <- ncol(tmp[[-2 * J]])
        m <- matrix(0, nrow = nr, ncol = nc)
        m <- tmp[[-2 * J]]
        if (scaling == "Haar") {
            x <- c(1:nr)/2^(-J) - 1
            y <- c(1:nc)/2^(-J) - 1
        }
        else if (scaling == "other") {
            x <- c(1:nr)
            y <- c(1:nc)
        }
        persp(x, y, m, theta = 30, phi = 30, expand =0.5, xlab = "tau_1", ylab = "tau_2", zlab = "ACW coeff", box=box)
        title(main)
    }
    else if (direction == 1) {
        tmp <- D2ACW(J = J, filter.number = filter.number, family = family, 
            OPLENGTH = OPLENGTH)
        nr <- nrow(tmp[[-1 * J]])
        nc <- ncol(tmp[[-1 * J]])
        m <- matrix(0, nrow = nr, ncol = nc)
        m <- tmp[[-1 * J]]
        if (scaling == "Haar") {
            x <- c(1:nr)/2^(-J) - 1
            y <- c(1:nc)/2^(-J) - 1
        }
        else if (scaling == "other") {
            x <- c(1:nr)
            y <- c(1:nc)
        }
        persp(x, y, m, theta = 30, phi = 30, expand =0.5, xlab = "tau_1", ylab = "tau_2", zlab = "ACW coeff", box=box)
        title(main)
    }
}

