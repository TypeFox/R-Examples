##' Two-stage (non-linear) measurement error
##'
##' Two-stage measurement error
##' @param model1 Stage 1 model
##' @param formula Formula specifying observed covariates in stage 2 model
##' @param data data.frame
##' @param predictfun Predictions to be used in stage 2
##' @param id1 Optional id-vector of stage 1
##' @param id2 Optional id-vector of stage 2
##' @param ... Additional arguments to lower level functions
##' @seealso stack.estimate
##' @export
##' @examples
##' m <- lvm(c(y1,y2,y3)~u,c(y3,y4,y5)~v,u~~v,c(u,v)~x)
##' transform(m,u2~u) <- function(x) x^2
##' transform(m,uv~u+v) <- prod
##' regression(m) <- z~u2+u+v+uv+x
##' set.seed(1)
##' d <- sim(m,1000,p=c("u,u"=1))
##' 
##' ## Stage 1
##' m1 <- lvm(c(y1[0:s],y2[0:s],y3[0:s])~1*u,c(y3[0:s],y4[0:s],y5[0:s])~1*v,u~b*x,u~~v)
##' latent(m1) <- ~u+v
##' e1 <- estimate(m1,d)
##' 
##' pp <- function(mu,var,data,...) {
##'     cbind(u=mu[,"u"],u2=mu[,"u"]^2+var["u","u"],v=mu[,"v"],uv=mu[,"u"]*mu[,"v"]+var["u","v"])
##' }
##' (e <- measurement.error(e1, z~1+x, data=d, predictfun=pp))
##' 
##' ## uu <- seq(-1,1,length.out=100)
##' ## pp <- estimate(e,function(p,...) p["(Intercept)"]+p["u"]*uu+p["u2"]*uu^2)$coefmat
##' if (interactive()) {
##'     plot(e,intercept=TRUE,vline=0)
##' 
##'     f <- function(p) p[1]+p["u"]*u+p["u2"]*u^2
##'     u <- seq(-1,1,length.out=100)
##'     plot(e, f, data=data.frame(u), ylim=c(-.5,2.5))
##' }
measurement.error <- function(model1, formula, data=parent.frame(), predictfun=function(mu,var,data,...) mu[,1]^2+var[1], id1, id2, ...) {
    if (!inherits(model1,c("lvmfit","lvm.mixture"))) stop("Expected lava object ('lvmfit','lvm.mixture',...)")
    if (missing(formula)) stop("formula needed for stage two (right-hand side additional covariates)")
    p1 <- coef(model1,full=TRUE)
    uhat <- function(p=p1) {
        P <- predict(model1,p=p,x=manifest(model1))
        cbind(predictfun(P,attributes(P)$cond.var,data))
    }
    if (missing(id1)) id1 <- seq(nrow(model.frame(model1)))
    if (missing(id2)) id2 <- seq(nrow(model.frame(model1)))
    if (!inherits(model1,"estimate"))
        e1 <- estimate(NULL,coef=p1,id=id1,iid=iid(model1))
    u <- uhat()
    X0 <- model.matrix(formula, data)
    Y <- model.frame(formula,data)[,1]
    X <- cbind(X0,u)
    stage.two <- lm(Y~-1+X)
    names(stage.two$coefficients) <- colnames(X)
    if (!inherits(stage.two,"estimate"))
        e2 <- estimate(stage.two, id=id2)
    U <- function(alpha=p1,beta=coef(stage.two)) {
        X <- cbind(X0,uhat(alpha))
        r <- (Y-X%*%beta)/summary(stage.two)$sigma^2
        apply(X,2,function(x) sum(x*r))
    }
    Ia <- -numDeriv::jacobian(function(p) U(p),p1)
    stacked <- stack(e1,e2,Ia)
    res <- c(stacked,list(naive=e2,lm=coef(summary(stage.two)),fun=predictfun))
    ## res <- list(estimate=stacked, naive=e2, lm=coef(summary(stage.two)),
    ##             fun=predictfun)
    structure(res,class=c("measurement.error","estimate"))
}
