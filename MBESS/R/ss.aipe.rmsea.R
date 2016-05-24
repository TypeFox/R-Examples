ss.aipe.rmsea <- function (RMSEA, df, width, conf.level = 0.95) 
{
    if (conf.level > 50 & conf.level < 100) 
        conf.level = conf.level/100
    if (conf.level < 0.5 | conf.level > 0.9999) 
        stop("The value of 'conf.level' must be between .5 and .999")
    
    omega <- width
    
    width.n <- function(RMSEA, df, omega, conf.level, N){
        ci.i <- ci.rmsea(rmsea = RMSEA, df = df, N = N, 
            conf.level = conf.level)
        w.i <- ci.i$Upper.Conf.Limit - ci.i$Lower.Conf.Limit
        return(abs(omega-w.i))
        }
    
    ci.2000 <- ci.rmsea(rmsea = RMSEA, df = df, N=2000, conf.level=conf.level)
    w.2000 <- ci.2000$Upper.Conf.Limit - ci.2000$Lower.Conf.Limit
    small.N <-ifelse (omega>w.2000, TRUE, FALSE) 
    med.N <- lg.N <- FALSE
    
    if(small.N){
        ss <-optimize(width.n, c(2, 2000), RMSEA=RMSEA, df=df, omega=omega, conf.level=conf.level)
        N <- ceiling(ss$minimum)
        }
    
    if (!small.N){
        ci.5000 <- ci.rmsea(rmsea = RMSEA, df = df, N=5000, conf.level=conf.level)
        w.5000 <- ci.5000$Upper.Conf.Limit - ci.5000$Lower.Conf.Limit
        med.N <- ifelse(omega>w.5000, TRUE, FALSE)
        lg.N <- ifelse(omega<w.5000, TRUE, FALSE)
        }

    if(med.N){
        ss <-optimize(width.n, c(2000,5000), RMSEA=RMSEA, df=df, omega=omega, conf.level=conf.level)
        N <- ceiling(ss$minimum)
        }
        
    if (lg.N){
        ss <- optim(par=6000, fn=width.n, RMSEA=RMSEA, df=df, omega=omega, conf.level=conf.level)
        N <- ceiling(ss$par)
        }
    
    ci.v <- ci.rmsea(rmsea = RMSEA, df = df, N=N, conf.level=conf.level)
    w.v<-ci.v$Upper.Conf.Limit - ci.v$Lower.Conf.Limit
    warn.msg <- paste("The sample size calculated is correct but is probably at the function break point.",
                      "\n","It is recommended to use function 'ci.rmsea' to verify this sample size."
                      )
    
    if(!isTRUE(all.equal(w.v, omega, tolerance=.001))){
        ci.v1 <- ci.rmsea(rmsea = RMSEA, df = df, N=N-1, conf.level=conf.level)
        w.v1<-ci.v1$Upper.Conf.Limit - ci.v1$Lower.Conf.Limit
        ci.v2 <- ci.rmsea(rmsea = RMSEA, df = df, N=N+1, conf.level=conf.level)
        w.v2<-ci.v2$Upper.Conf.Limit - ci.v2$Lower.Conf.Limit
        if(w.v1<omega | w.v2>omega){
            ss <-optimize(width.n, c(2, 20000), RMSEA=RMSEA, df=df, omega=omega, conf.level=conf.level)
            N <- ceiling(ss$minimum)
            }#else warning(warn.msg,immediate. = TRUE)
        }
    
    cat("Necessary sample size so that the expected width of the ", 
        conf.level * 100, "% confidence interval","\n", "is no greater than ", 
        width, ", given a population RMSEA of ", RMSEA, ", is:", "\n",sep="" 
        )
    return(N)
}
