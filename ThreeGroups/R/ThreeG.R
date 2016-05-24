ThreeG <- function(data, starting = NULL)	{

# check input arguments for data
if(class(data) == "matrix") stop("data cannot be a matrix")
if(class(data) == "data.frame") stop("data cannot be a data frame")
if(mode(data) != "list" & mode(data) != "numeric") stop("mode of data must be list or numeric")

if(length(data) != 10) stop("data must contain 10 elements")
if(any(is.na(data))) stop("data must not contain missing values")
data <- as.numeric(data)
if(any(is.na(data))) stop("data must be numeric")

if(any(data < 0 | data == Inf)) stop("data must contain non-negative integers")
if(any(abs(data - round(data)) > .Machine$double.eps^0.5)) stop("data must contain non-negative integers")

# check input arguments for starting
if(class(starting) == "matrix") {
warning("starting cannot be a matrix; generating new starting values"); starting <- NULL }
if(class(starting) == "data.frame") {
warning("starting cannot be a data frame; generating new starting values"); starting <- NULL }
if(mode(starting) != "numeric" & mode(starting) != "NULL") {
warning("mode of starting must be numeric; generating new starting values"); starting <- NULL }

if(mode(starting) != "NULL")    { 
    if(length(starting) != 4) stop("starting must contain 4 elements")
    starting <- as.numeric(starting)
    if(any(is.na(starting))) stop("starting must not contain missing values")

    if(any(starting[1:3] < 0 | starting[1:3] > 1)) {
    stop("starting values for probabilities are bounded by 0 and 1") }
    if((starting[4] < -1 | starting[4] > 1)) {
    stop("starting values for complier average treatment effect are bounded by -1 and 1") }

                                }

    starting <- c((data[7] + data[8]) / (data[9] + data[10] + data[7] + data[8]),
                (data[3]) / (data[3] + data[4]),
                (data[5] + data[9]) / (data[6] + data[10] + data[5] + data[9]),
                (data[3] / (data[3] + data[4])) - (data[7] / (data[7] + data[8])))

# optim
optim.out <- suppressWarnings(optim(par = starting, .ThreeGLL, method = "L-BFGS-B",
					lower = c(0,0,0,-1), upper = c(1,1,1,1),
					control = list(fnscale = -1, 
                    reltol = .Machine$double.eps, trace = 0, maxit = 100000), hessian = TRUE,
                    s_b    = data[1],
                    f_b    = data[2],
                    s_t_c  = data[3],
                    f_t_c  = data[4],
                    s_t_nc = data[5],
                    f_t_nc = data[6],
                    s_p_c  = data[7],
                    f_p_c  = data[8],
                    s_p_nc = data[9],
                    f_p_nc = data[10]))

if(optim.out$convergence != 0) warning("optimizer failed to converge. Try different starting values")

est <- cbind(optim.out$par, sqrt(diag(solve(-optim.out$hessian))))
colnames(est) <- c("estimates","standard errors")
res <- list(starting = starting, est = est, optim.out = optim.out)

cat("\n","\n")
cat("Three-group estimator","\n")
cat("---------------------","\n","\n")

cat("Est. proportion of compliers (SE).............................. ", round(res$est[1,1],4),
paste("(",round(res$est[1,2],4),")",sep = ""),"\n")
cat("Est. probability of success given treatment for compliers (SE). ", round(res$est[2,1],4),
paste("(",round(res$est[2,2],4),")",sep = ""),"\n")
cat("Est. probability of success for non-compliers (SE)............. ", round(res$est[3,1],4),
paste("(",round(res$est[3,2],4),")",sep = ""),"\n")
cat("Est. complier average treatment effect (SE).................... ", round(res$est[4,1],4),
paste("(",round(res$est[4,2],4),")",sep = ""),"\n")

return(invisible(res))

    }
