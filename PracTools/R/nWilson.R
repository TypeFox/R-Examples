
nWilson <- function(moe.sw,alpha=0.05,pU,e){
    if (!(moe.sw==1 | moe.sw==2))
        stop("moe.sw must equal 1 or 2.\n")
    if (e <=0 | e >= 1) stop("e must be in (0,1).\n")
    if (pU <=0 | pU >= 1) stop("pU must be in (0,1).\n")
    
    za <- qnorm(1-alpha/2)
    qU <- 1-pU
    if (moe.sw == 1){
        rad <- e^2 - pU*qU * (4*e^2 - pU*qU)
    }
    if (moe.sw == 2){
        e <- e*pU
        rad <- e^2 - pU*qU * (4*e^2 - pU*qU)
    }

    n.sam <- (pU*qU - 2*e^2 + sqrt(rad) ) * (za/e)^2 / 2
    
    d <- za*sqrt(za^2 + 4*n.sam*pU*qU)
    CI.L <- (2*n.sam*pU + za^2 - d)/2/(n.sam + za^2)
    CI.U <- (2*n.sam*pU + za^2 + d)/2/(n.sam + za^2)
    leng.CI <- CI.U - CI.L

    list(n.sam=n.sam, "CI lower limit"=CI.L, "CI upper limit"=CI.U, "length of CI" = leng.CI)
}
