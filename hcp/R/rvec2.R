rvec2 <-
function(x, y, n, j, k, e1, e2)
{
	R <- array(0, 5)
	jp1 <- j + 1
	kp1 <- k + 1
	y1j <- sum(y[1:j])
	yjk <- sum(y[jp1:k])
	ykn <- sum(y[kp1:n])
	xy1j <- sum(x[1:j] * y[1:j])
	xyjk <- sum(x[jp1:k] * y[jp1:k])
        x2yjk <- sum(x[jp1:k]^2 * y[jp1:k])
	xykn <- sum(x[kp1:n] * y[kp1:n])
	R[1] <- e1 * (y1j + ykn) + e2 * yjk
	R[2] <- e1 * (xy1j + x[j] * ykn) + e2 * x[j] * yjk
	R[3] <- e1 * (x[k] - x[j]) * ykn + e2 * (xyjk - x[j] * yjk)
        R[4] <- e1 * (x[k]^2 - x[j]^2) * ykn + e2 * (x2yjk - x[j]^2 * yjk) 
	R[5] <- e1 * (xykn - x[k] * ykn)
	list(r = R)
}
