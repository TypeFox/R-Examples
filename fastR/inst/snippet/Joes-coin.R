wald.ci <- function(x,n,level=0.95) {
    alpha = 1 - level;
    pi.hat <- x/n;
    se <- sqrt(pi.hat * (1-pi.hat) / n);
    z.star <- qnorm(1- alpha/2);
    pi.hat + c(-1,1) * z.star * se;
}
wilson.ci<- function(x,n,level=0.95) {
    x = x +2; n = n+4;
    alpha = 1 - level;
    pi.hat <- x/n;
    se <- sqrt(pi.hat * (1-pi.hat) / n);
    z.star <- qnorm(1- alpha/2);
    pi.hat + c(-1,1) * z.star * se;
}
score.ci <- function(x,n,level=0.95) {
    alpha = 1 - level;
    z.star <- qnorm(1 - alpha/2);
    pi.hat <- x/n;
    A <- pi.hat + z.star^2/(2*n);
    B <- z.star * sqrt( ( pi.hat * (1-pi.hat) / n ) 
	                  + ( z.star^2/(4 * n^2) ) );
    D <- 1 + z.star^2/n;
    # interval is ( A +- B) / D
    ( A + c(-1,1) * B ) / D;
}
prop.test(115,200)$conf.int;
prop.test(115,200,correct=FALSE)$conf.int;
wald.ci(115,200);
wilson.ci(115,200);
wald.ci(117,204);
score.ci(115,200);

# score interval using uniroot:
p.hat <- 115/200; n <- 200;
f <- function(p) {
    abs(p.hat - p) / sqrt(p * (1-p) / n) + qnorm(0.025);
}
uniroot(f,c(0,p.hat))$root;
uniroot(f,c(p.hat,1))$root;
uniroot(f,c(0,p.hat))$estim.prec
