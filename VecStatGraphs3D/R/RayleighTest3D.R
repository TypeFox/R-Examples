RayleighTest3D <- function (coord, Alpha = 0.05) 
{
    TableChicuadrado <- c(6.251, 7.815, 9.348, 11.34, 12.84, 16.27)
    x <- coord[, 1]
    y <- coord[, 2]
    z <- coord[, 3]
    n_elements = length(x)
    m_module = MeanModule3D(coord)
    m_module = 3 * m_module
	
    if (sum(Alpha == c(0.10, 0.05, 0.025, 0.01, 0.005, 0.001)) == 0) {
        print("Alpha value is not available, use 0.10, 0.05, 0.025, 0.01, 0.005, or 0.001")
    }
    else {
        Tcol <- (1:6)[Alpha == c(0.10, 0.05, 0.025, 0.01, 0.005, 0.001)]
		print("Rayleigh test of uniformity")
        if (n_elements >= 10) {
            if (m_module < TableChicuadrado[Tcol]) {
                print(paste("Hypothesis of uniformity accepted, P value =", Alpha))
            }
            else {
                print(paste("Hypothesis of uniformity rejected, P value =", Alpha))
            }
        }
        else {
            print("Rayleigh Test cannot be performed: sample size is too small")
        }
    }
}
