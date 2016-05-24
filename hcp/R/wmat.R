wmat <-
function(x, y, n, j, k, e1, e2)
{
	W <- matrix(0, 4, 4)
	jp1 <- j + 1
	kp1 <- k + 1
	W[1, 1] <- e1 * (n - k + j) + e2 * (k - j)
	W[1, 2] <- e1 * (sum(x[1:j]) + (n - k) * x[j]) + e2 * (k - j) * x[j]
	W[1, 3] <- e1 * (n - k) * (x[k] - x[j]) + e2 * sum(x[jp1:k] - x[j])
	W[1, 4] <- e1 * sum(x[kp1:n] - x[k])
	W[2, 2] <- e1 * (sum(x[1:j] * x[1:j]) + (n - k) * x[j] * x[j]) + e2 * (k - j) *
 		x[j] * x[j]
	W[2, 3] <- e1 * (n - k) * x[j] * (x[k] - x[j]) + e2 * x[j] * sum(x[jp1:k] - x[
		j])
	W[2, 4] <- e1 * x[j] * sum(x[kp1:n] - x[k])
	W[3, 3] <- e1 * (n - k) * (x[k] - x[j]) * (x[k] - x[j]) + e2 * sum((x[jp1:k] - 
		x[j]) * (x[jp1:k] - x[j]))
	W[3, 4] <- e1 * (x[k] - x[j]) * sum(x[kp1:n] - x[k])
	W[4, 4] <- e1 * sum((x[kp1:n] - x[k]) * (x[kp1:n] - x[k]))
	W[2, 1] <- W[1, 2]
	W[3, 1] <- W[1, 3]
	W[4, 1] <- W[1, 4]
	W[3, 2] <- W[2, 3]
	W[4, 2] <- W[2, 4]
	W[4, 3] <- W[3, 4]
	list(w = W)
}
