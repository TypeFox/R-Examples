### R code from vignette source 'mob.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("partykit")
options(prompt = "R> ", continue = "+  ", digits = 4, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: PimaIndiansDiabetes
###################################################
data("PimaIndiansDiabetes", package = "mlbench")


###################################################
### code chunk number 3: PimIndiansDiabetes-formula
###################################################
pid_formula <- diabetes ~ glucose | pregnant + pressure + triceps +
  insulin + mass + pedigree + age


###################################################
### code chunk number 4: logit
###################################################
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}


###################################################
### code chunk number 5: PimaIndiansDiabetes-mob
###################################################
pid_tree <- mob(pid_formula, data = PimaIndiansDiabetes, fit = logit)


###################################################
### code chunk number 6: PimaIndiansDiabetes-print
###################################################
pid_tree


###################################################
### code chunk number 7: PimaIndiansDiabetes-glmtree
###################################################
pid_tree2 <- glmtree(diabetes ~ glucose | pregnant +
  pressure + triceps + insulin + mass + pedigree + age,
  data = PimaIndiansDiabetes, family = binomial)


###################################################
### code chunk number 8: PimaIndiansDiabetes-plot
###################################################
plot(pid_tree)


###################################################
### code chunk number 9: PimaIndiansDiabetes-plot2
###################################################
plot(pid_tree2, tp_args = list(ylines = 1, margins = c(1.5, 1.5, 1.5, 2.5)))


###################################################
### code chunk number 10: PimaIndiansDiabetes-sctest1
###################################################
library("strucchange")
sctest(pid_tree, node = 1)


###################################################
### code chunk number 11: PimaIndiansDiabetes-sctest2
###################################################
sctest(pid_tree, node = 2)


###################################################
### code chunk number 12: PimaIndiansDiabetes-sctest3
###################################################
sctest(pid_tree, node = 3)


###################################################
### code chunk number 13: PimaIndiansDiabetes-prune (eval = FALSE)
###################################################
## pid_tree3 <- mob(pid_formula, data = PimaIndiansDiabetes,
##   fit = logit, control = mob_control(verbose = TRUE,
##     minsize = 50, maxdepth = 4, alpha = 0.9, prune = "BIC"))


###################################################
### code chunk number 14: PimaIndiansDiabetes-info
###################################################
names(pid_tree$info)


###################################################
### code chunk number 15: PimaIndiansDiabetes-info
###################################################
names(pid_tree$node$info)


###################################################
### code chunk number 16: PimaIndiansDiabetes-print3
###################################################
print(pid_tree, node = 3)


###################################################
### code chunk number 17: PimaIndiansDiabetes-coef
###################################################
coef(pid_tree)
coef(pid_tree, node = 1)
summary(pid_tree, node = 1)


###################################################
### code chunk number 18: mob.Rnw:783-784
###################################################
exp(coef(pid_tree)[,2])


###################################################
### code chunk number 19: mob.Rnw:787-788
###################################################
risk <- round(100 * (exp(coef(pid_tree)[,2])-1), digits = 1)


###################################################
### code chunk number 20: PimaIndiansDiabetes-logLik
###################################################
logLik(pid_tree)
AIC(pid_tree)
BIC(pid_tree)


###################################################
### code chunk number 21: PimaIndiansDiabetes-deviance
###################################################
mean(residuals(pid_tree)^2)
deviance(pid_tree)/sum(weights(pid_tree))
deviance(pid_tree)/nobs(pid_tree)


###################################################
### code chunk number 22: PimaIndiansDiabetes-predict
###################################################
pid <- head(PimaIndiansDiabetes)
predict(pid_tree, newdata = pid, type = "node")


###################################################
### code chunk number 23: PimaIndiansDiabetes-width
###################################################
width(pid_tree)
depth(pid_tree)


###################################################
### code chunk number 24: PimaIndiansDiabetes-subset
###################################################
pid_tree[3]


###################################################
### code chunk number 25: mob.Rnw:875-878
###################################################
predict(pid_tree2, newdata = pid, type = "node")
predict(pid_tree2, newdata = pid, type = "response")
predict(pid_tree2, newdata = pid, type = "link")


###################################################
### code chunk number 26: Journals-data
###################################################
data("Journals", package = "AER")
Journals <- transform(Journals,
  age = 2000 - foundingyear,
  chars = charpp * pages)


###################################################
### code chunk number 27: Journals-tree
###################################################
j_tree <- lmtree(log(subs) ~ log(price/citations) | price + citations +
  age + chars + society, data = Journals, minsize = 10, verbose = TRUE)


###################################################
### code chunk number 28: Journals-plot
###################################################
plot(j_tree)


###################################################
### code chunk number 29: Journals-print
###################################################
j_tree


###################################################
### code chunk number 30: Journals-methods
###################################################
coef(j_tree, node = 1:3)
summary(j_tree, node = 1:3)
sctest(j_tree, node = 1:3)


###################################################
### code chunk number 31: BostonHousing-data
###################################################
data("BostonHousing", package = "mlbench")
BostonHousing <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))


###################################################
### code chunk number 32: BostonHousing-tree
###################################################
bh_tree <- lmtree(medv ~ log(lstat) + I(rm^2) | zn + indus + chas + nox +
  age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
bh_tree


###################################################
### code chunk number 33: BostonHousing-plot
###################################################
plot(bh_tree)


###################################################
### code chunk number 34: BostonHousing-AIC
###################################################
mean(residuals(bh_tree)^2)
logLik(bh_tree)
AIC(bh_tree)


###################################################
### code chunk number 35: TeachingRatings-data
###################################################
data("TeachingRatings", package = "AER")
tr <- subset(TeachingRatings, credits == "more")


###################################################
### code chunk number 36: TeachingRatings-lm
###################################################
tr_null <- lm(eval ~ 1, data = tr, weights = students)
tr_lm <- lm(eval ~ beauty + gender + minority + native + tenure + division,
  data = tr, weights = students)


###################################################
### code chunk number 37: TeachingRatings-tree
###################################################
tr_tree <- lmtree(eval ~ beauty | minority + age + gender + division +
  native + tenure, data = tr, weights = students, caseweights = FALSE)


###################################################
### code chunk number 38: TeachingRatings-plot
###################################################
plot(tr_tree)


###################################################
### code chunk number 39: TeachingRatings-coef
###################################################
coef(tr_lm)[2]
coef(tr_tree)[, 2]


###################################################
### code chunk number 40: TeachingRatings-rsquared
###################################################
1 - c(deviance(tr_lm), deviance(tr_tree))/deviance(tr_null)


###################################################
### code chunk number 41: Titanic-data
###################################################
data("Titanic", package = "datasets")
ttnc <- as.data.frame(Titanic)
ttnc <- ttnc[rep(1:nrow(ttnc), ttnc$Freq), 1:4]
names(ttnc)[2] <- "Gender"
ttnc <- transform(ttnc, Treatment = factor(
  Gender == "Female" | Age == "Child", levels = c(FALSE, TRUE),
  labels = c("Male&Adult", "Female|Child")))


###################################################
### code chunk number 42: Titanic-tree
###################################################
ttnc_tree <- glmtree(Survived ~ Treatment | Class + Gender + Age,
  data = ttnc, family = binomial, alpha = 0.01)
ttnc_tree


###################################################
### code chunk number 43: Titanic-plot
###################################################
plot(ttnc_tree, tp_args = list(ylines = 1, margins = c(1.5, 1.5, 1.5, 2.5)))


###################################################
### code chunk number 44: GBSG2
###################################################
data("GBSG2", package = "TH.data")
GBSG2$time <- GBSG2$time/365


###################################################
### code chunk number 45: wbreg
###################################################
library("survival")
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}


###################################################
### code chunk number 46: logLik.survreg
###################################################
logLik.survreg <- function(object, ...)
  structure(object$loglik[2], df = sum(object$df), class = "logLik")


###################################################
### code chunk number 47: gbsg2_tree
###################################################
gbsg2_tree <- mob(Surv(time, cens) ~ horTh + pnodes | age + tsize +
  tgrade + progrec + estrec + menostat, data = GBSG2,
  fit = wbreg, control = mob_control(minsize = 80))


###################################################
### code chunk number 48: GBSG2-plot
###################################################
plot(gbsg2_tree)


###################################################
### code chunk number 49: GBSG2-scatter
###################################################
gbsg2node <- function(mobobj, 
  col = "black", linecol = "red", cex = 0.5, pch = NULL,
  jitter = FALSE, xscale = NULL, yscale = NULL, ylines = 1.5,
  id = TRUE, xlab = FALSE, ylab = FALSE)
{
  ## obtain dependent variable
  mf <- model.frame(mobobj)
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, rhs = 0L)
  if(isTRUE(ylab)) ylab <- names(y)[1L]
  if(identical(ylab, FALSE)) ylab <- ""
  if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
  y <- y[[1L]]

  ## plotting character and response
  if(is.null(pch)) pch <- y[,2] * 18 + 1
  y <- y[,1]
  y <- as.numeric(y)
  pch <- rep(pch, length.out = length(y))
  if(jitter) y <- jitter(y)

  ## obtain explanatory variables
  x <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, rhs = 1L)
  xnam <- colnames(x)
  z <- seq(from = min(x[,2]), to = max(x[,2]), length = 51)
  z <- data.frame(a = rep(sort(x[,1])[c(1, NROW(x))], c(51, 51)), b = z)
  names(z) <- names(x)
  z$x <- model.matrix(~ ., data = z)
  
  ## fitted node ids
  fitted <- mobobj$fitted[["(fitted)"]]
      
  if(is.null(xscale)) xscale <- range(x[,2]) + c(-0.1, 0.1) * diff(range(x[,2]))
  if(is.null(yscale)) yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
       
  ## panel function for scatter plots in nodes
  rval <- function(node) {
  
    ## node index
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)

    ## dependent variable
    y <- y[ix]

    ## predictions
    yhat <- if(is.null(node$info$object)) {
      refit.modelparty(mobobj, node = nid)
    } else {
      node$info$object
    }
    yhat <- predict(yhat, newdata = z, type = "quantile", p = 0.5)
    pch <- pch[ix]

    ## viewport setup
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
        	       widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
  		       heights = unit(c(1, 1), c("lines", "null"))),
        	       width = unit(1, "npc"), 
        	       height = unit(1, "npc") - unit(2, "lines"),
  		       name = paste("node_scatterplot", nid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = "white", col = 0))

    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    mainlab <- paste(ifelse(id, paste("Node", nid, "(n = "), ""),
  		     info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
    grid.text(mainlab)
    popViewport()

    plot_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, xscale = xscale,
  	yscale = yscale, name = paste("node_scatterplot", nid, "plot", sep = ""))
    pushViewport(plot_vp)

    ## scatterplot
    grid.points(x[ix,2], y, gp = gpar(col = col, cex = cex), pch = pch)
    grid.lines(z[1:51,2], yhat[1:51], default.units = "native", gp = gpar(col = linecol))
    grid.lines(z[52:102,2], yhat[52:102], default.units = "native", gp = gpar(col = linecol, lty = 2))

    grid.xaxis(at = c(ceiling(xscale[1]*10), floor(xscale[2]*10))/10)
    grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))

    if(isTRUE(xlab)) xlab <- xnam[2]
    if(!identical(xlab, FALSE)) grid.text(xlab, x = unit(0.5, "npc"), y = unit(-2, "lines"))
    if(!identical(ylab, FALSE)) grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2, "lines"), rot = 90)

    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()

    upViewport()
  }
          
  return(rval)
}
class(gbsg2node) <- "grapcon_generator"

plot(gbsg2_tree, terminal_panel = gbsg2node, tnex = 2, 
  tp_args = list(xscale = c(0, 52), yscale = c(-0.5, 8.7)))


###################################################
### code chunk number 50: gbsg2_tree-methods
###################################################
gbsg2_tree
coef(gbsg2_tree)
logLik(gbsg2_tree)


###################################################
### code chunk number 51: bt-packages
###################################################
data("Topmodel2007", package = "psychotree")
library("psychotools")
estfun.btReg <- function(x, ...) x$estfun


###################################################
### code chunk number 52: btfit1
###################################################
btfit1 <- function(y, x = NULL, start = NULL, weights = NULL,
  offset = NULL, ...) btReg.fit(y, ...)


###################################################
### code chunk number 53: bt1
###################################################
system.time(bt1 <- mob(preference ~ 1 | gender + age + q1 + q2 + q3,
  data = Topmodel2007, fit = btfit1))


###################################################
### code chunk number 54: btfit2
###################################################
btfit2 <- function(y, x = NULL, start = NULL, weights = NULL,
  offset = NULL, ..., estfun = FALSE, object = FALSE) {
  rval <- btReg.fit(y, ..., estfun = estfun, vcov = object)
  list(
    coefficients = rval$coefficients,
    objfun = -rval$loglik,
    estfun = if(estfun) rval$estfun else NULL,
    object = if(object) rval else NULL
  )
}


###################################################
### code chunk number 55: bt2
###################################################
system.time(bt2 <- mob(preference ~ 1 | gender + age + q1 + q2 + q3,
  data = Topmodel2007, fit = btfit2))


###################################################
### code chunk number 56: bt2-print
###################################################
bt2
coef(bt2)


###################################################
### code chunk number 57: bt2-worthf (eval = FALSE)
###################################################
## worthf <- function(info) paste(info$object$labels,
##   format(round(worth(info$object), digits = 3)), sep = ": ")
## plot(bt2, FUN = worthf)


###################################################
### code chunk number 58: bt2-plot
###################################################
plot(bt2)


###################################################
### code chunk number 59: bt2-plot2
###################################################
worthf <- function(info) paste(info$object$labels,
  format(round(worth(info$object), digits = 3)), sep = ": ")
plot(bt2, FUN = worthf)


###################################################
### code chunk number 60: bt2-nodeapply (eval = FALSE)
###################################################
## par(mfrow = c(2, 2))
## nodeapply(bt2, ids = c(3, 5, 6, 7), FUN = function(n)
##   plot(n$info$object, main = n$id, ylim = c(0, 0.4)))


###################################################
### code chunk number 61: bt2-nodeapply-plot
###################################################
par(mfrow = c(2, 2))
nodeapply(bt2, ids = c(3, 5, 6, 7), FUN = function(n)
  plot(n$info$object, main = n$id, ylim = c(0, 0.4)))


###################################################
### code chunk number 62: bt2-plot3 (eval = FALSE)
###################################################
## plot(bt2, drop = TRUE, tnex = 2,
##   terminal_panel = node_btplot(bt2, abbreviate = 1, yscale = c(0, 0.5)))


###################################################
### code chunk number 63: node_btplot
###################################################
## visualization function
node_btplot <- function(mobobj, id = TRUE,
  worth = TRUE, names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## node ids
    node <- nodeids(mobobj, terminal = FALSE)
    
    ## get all coefficients 
    cf <- partykit:::apply_to_models(mobobj, node, FUN = function(z)        
      if(worth) worth(z) else coef(z, all = FALSE, ref = TRUE))
    cf <- do.call("rbind", cf)
    rownames(cf) <- node

    ## get one full model
    mod <- partykit:::apply_to_models(mobobj, node = 1L, FUN = NULL)

    if(!worth) {
      if(is.character(ref) | is.numeric(ref)) {
        reflab <- ref
        ref <- TRUE
      } else {
        reflab <- mod$ref
      }
      if(is.character(reflab)) reflab <- match(reflab, mod$labels)
      cf <- cf - cf[,reflab]
    }

    ## reference
    if(worth) {
      cf_ref <- 1/ncol(cf)
    } else {
      cf_ref <- 0
    }

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)
    
    if(index) {
      x <- 1:NCOL(cf)
      if(is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    } else {
      x <- rep(0, length(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)      
    }
    if(is.null(yscale)) yscale <- range(cf) + c(-0.1, 0.1) * diff(range(cf))
         
    ## panel function for bt plots in nodes
    rval <- function(node) {

      ## node index
      id <- id_node(node)
    
      ## dependent variable setup
      cfi <- cf[id,]

      ## viewport setup
      top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
    			 widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
        		 heights = unit(c(1, 1), c("lines", "null"))),
    			 width = unit(1, "npc"), 
    			 height = unit(1, "npc") - unit(2, "lines"),
        		 name = paste("node_btplot", id, sep = ""))
      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "white", col = 0))

      ## main title
      top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
      pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", id, "(n = "), ""),
        	       info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
      grid.text(mainlab)
      popViewport()

      ## actual plot  
      plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
        xscale = xscale, yscale = yscale, 
        name = paste("node_btplot", id, "plot", sep = ""))
      pushViewport(plot_vpi)

      grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = linecol), default.units = "native")
      if(index) {
        grid.lines(x, cfi, gp = gpar(col = col, lty = 2), default.units = "native")
        grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
        grid.xaxis(at = x, label = if(names) names(cfi) else x)
      } else {  	
        if(names) grid.text(names(cfi), x = x, y = cfi, default.units = "native")
          else grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
      }
      grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
      grid.rect(gp = gpar(fill = "transparent"))

      upViewport(2)
    }
	    
    return(rval)
}
class(node_btplot) <- "grapcon_generator"

plot(bt2, drop = TRUE, tnex = 2,
  terminal_panel = node_btplot(bt2, abbreviate = 1, yscale = c(0, 0.5)))


###################################################
### code chunk number 64: tm
###################################################
tm <- data.frame(age = c(60, 25, 35), gender = c("male", "female", "female"),
  q1 = "no", q2 = c("no", "no", "yes"), q3 = "no")
tm


###################################################
### code chunk number 65: tm-predict
###################################################
tm
predict(bt2, tm, type = "node")
predict(bt2, tm, type = function(object) t(worth(object)))
predict(bt2, tm, type = function(object) t(rank(-worth(object))))


