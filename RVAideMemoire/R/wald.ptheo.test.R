# lme4 : fixef, glmer, VarCorr

wald.ptheo.test <-
function(y,blocks=NULL,p=0.5) {
  if (p<=0 | p>=1) {stop("wrong p: 0 < p < 1")}
  namey <- deparse(substitute(y))
  nameb <- deparse(substitute(blocks))
  y <- if (is.vector(y)) {
    as.numeric(factor(y))-1
  } else {
    as.matrix(y)
  }
  overdisp <- FALSE
  param <- NULL
  if (!is.null(blocks)) {
    model1 <- lme4::glmer(y~1+(1|blocks),family="binomial")
    vpars <- function(m) {
	nrow(m)*(nrow(m)+1)/2
    }
    model1.df <- sum(sapply(lme4::VarCorr(model1),vpars))+length(lme4::fixef(model1))
    rdf <- nrow(model.frame(model1))-model1.df
    rp <- residuals(model1)
    dev <- sum(rp^2)
    if (dev>rdf) {
	ind <- if (is.vector(y)) {
	  factor(1:length(y))
	} else {
	  factor(1:nrow(y))
	}
	model <- suppressWarnings(lme4::glmer(y~1+(1|blocks)+(1|ind),family="binomial"))
	overdisp <- TRUE
    } else {
	model <- model1
    }
    p.est <- c("probability of success"=unique(predict(model,type="response",re.form=NA)))
    dname <- paste0(namey," and ",nameb)
  } else {
    model1 <- glm(y~1,family="binomial")
    summ.mod1 <- summary(model1)
    if (summ.mod1$deviance>summ.mod1$df.residual) {
	model <- glm(y~1,family="quasibinomial")
	overdisp <- TRUE
    } else {
	model <- model1
    }
    p.est <- c("probability of success"=unique(predict(model,type="response")))
    dname <- namey
  }
  int.theo <- log(p/(1-p))
  summ <- summary(model)
  int.moy <- summ$coefficients["(Intercept)","Estimate"]
  int.se <- summ$coefficients["(Intercept)","Std. Error"]
  if (!is.null(blocks)) {
    stat <- c(z=(int.moy-int.theo)/int.se)
    pval <- min(pnorm(stat),pnorm(stat,lower.tail=FALSE))*2
  } else {
    if (overdisp) {
	stat <- c(t=(int.moy-int.theo)/int.se)
	param <- c(df=summ$df.residual)
	pval <- 2*min(pt(stat,param),pt(stat,param,lower.tail=FALSE))
    } else {
	stat <- c(z=(int.moy-int.theo)/int.se)
	pval <- 2*min(pnorm(stat),pnorm(stat,lower.tail=FALSE))
    }
  }
  if (overdisp) {warning("overdispersion taken into account")}
  null <- c("probability of success"=p)
  res <- list(statistic=stat,p.value=pval,estimate=p.est,null.value=null,
    alternative="two.sided",method="Wald test",data.name=dname)
  if (!is.null(param)) {res$parameter=param}
  class(res) <- "htest"
  return(res)
}
