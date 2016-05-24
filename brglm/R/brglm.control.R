`brglm.control` <-
function (br.epsilon = 1e-08, br.maxit = 100, br.trace = FALSE,
          br.consts = NULL, ...) 
{
    if (!is.numeric(br.epsilon) || br.epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(br.maxit) || br.maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(br.epsilon = br.epsilon, br.maxit = br.maxit, br.trace = br.trace,
         br.consts = br.consts)
}
