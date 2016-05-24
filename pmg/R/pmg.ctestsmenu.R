#### Classical testse (ctests)
#### This file defines the list that get fed to genericWidget for
#### making dialog

### sample file for t.test.
### abstract so that you can produce these dialogs from a list

## prop.test
prop.test.list = list(
  title = "prop.test()",
  help = "prop.test",
  type = "text",                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "prop.test(",
    ending = ")"
    ),
  arguments = list(
    variables = list(
      x=EMPTY.list,
      n=EMPTY.list
      ),
    hypotheses = list(
      alternative= alternative.list,
      p = list(
        type = "gedit",
        text = "0.5"
        )
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )

##binom.test
## prop.test
binom.test.list = list(
  title = "binom.test()",
  help = "binom.test",
  type = "text",                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "binom.test(",
    ending = ")"
    ),
  arguments = list(
    variables = list(
      x=EMPTY.list,
      n=EMPTY.list
      ),
    hypotheses = list(
      alternative= alternative.list,
      p = list(
        type = "gedit",
        text = "0"
        )
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )


## t.test
t.test.list = list(
  title = "t.test()",
  help = "t.test",
  type = "text",                      # either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "t.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      mu = list(
        type = "gedit",
        text = "0"
        ),
      alternative = alternative.list,
      paired = FALSE.list,
      var.equal = FALSE.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )

## ilcox test
wilcox.test.list = list(
  title = "wilcox.test()",
  help = "wilcox.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "wilcox.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      mu = list(
        type = "gedit",
        text = "0"
        ),
      alternative = alternative.list,
      paired = FALSE.list,
      var.equal = FALSE.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )


## var.test
var.test.list = list(
  title = "var.test()",
  help = "var.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "var.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      alternative= alternative.list,
      ratio = list(
        type = "gedit",
        text = 1
        )
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )



## bartlett
bartlett.test.list = list(
  title = "bartlett.test()",
  help = "bartlett.test",
  type = "text",                        #either text or graphic
  variableType = "model",
  assignto = NULL,
  action = list(
    beginning = "bartlett.test(",
    ending = ")"
    ),
  arguments = list(
    )
  )



## fligner
fligner.test.list = list(
  title = "fligner.test()",
  help = "fligner.test",
  type = "text",                        #either text or graphic
  variableType = "model",
  assignto = NULL,
  action = list(
    beginning = "fligner.test(",
    ending = ")"
    ),
  arguments = list(
    )
  )


## ansari.test
ansari.test.list = list(
  title = "ansari.test()",
  help = "ansari.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "ansari.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      alternative= alternative.list,
      exact = list(
        type = "gedit",
        text = ""
        )
      ),
    "CI" = list(
      conf.int = FALSE.list,
      conf.level = conf.level.list
      )
    )
  )

## cor.test
cor.test.list = list(
  title = "cor.test()",
  help = "cor.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "cor.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      alternative= alternative.list,
      method= list(
        type="gdroplist",
        items=c("\"pearson\"","\"kendall\"","\"spearman\"")
        ),
      exact = EMPTY.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )


## chisq.test
chisq.test.list = list(
  title = "chisq.test()",
  help = "chisq.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "chisq.test(",
    ending = ")"
    ),
  arguments = list(
    "hypotheses" = list(
      correct = FALSE.list,
      p = list(
        type = "gedit",
        text = "0"
        )
      ),
    "calculate"=list(
      simulate.p.value = FALSE.list
      )
    )
  )

## mcnemar
mcnemar.test.list = list(
  title = "mcnemar.test()",
  help = "mcnemar.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "mcnemar.test(",
    ending = ")"
    ),
  arguments = list(
    hypotheses = list(
      correct = TRUE.list
      )
    ),
  "CI" = list(
    conf.level = conf.level.list
    )
  )



## mantelhaen.tes
mantelhaen.test.list =  list(
  title = "mantelhaen.test()",
  help = "mantelhaen.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "mantelhaen.test(",
    ending = ")"
    ),
  arguments = list(
    "hypotheses" = list(
      alternative= alternative.list,
      correct = TRUE.list,
      exact = FALSE.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )

## kolmogorov smirnov
ks.test.list = list(
  title = "ks.test()",
  help = "ks.test",
  type = "text",                        #either text or graphic
  variableType = "bivariate",
  assignto = NULL,
  action = list(
    beginning = "ks.test(",
    ending = ")"
    ),
  arguments = list(
    "Note:" = list(
      label = list(
        type="glabel",
        text = "Only does two sample ks.test"
        )
      ),
    "hypotheses" = list(
      alternative= alternative.list,
      exact = EMPTY.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )


  )

## shapiro
shapiro.test.list = list(
  title = "shapiro.test()",
  help = "shapiro.test",
  type = "text",                        #either text or graphic
  variableType = "univariate",
  assignto = NULL,
  action = list(
    beginning = "shapiro.test(",
    ending = ")"
    )
  )

## oneway test
oneway.test.list = list(
  title = "oneway.test()",
  help = "oneway.test",
  type = "text",                        #either text or graphic
  variableType = "model",
  assignto = NULL,
  action = list(
    beginning = "oneway.test(",
    ending = ")"
    ),
  arguments = list(
    )
  )


### kruskal test
kruskal.test.list= list(
  title = "kruskal.test()",
  help = "kruskal.test",
  type = "text",                        #either text or graphic
  variableType = "model",
  assignto = NULL,
  action = list(
    beginning = "kruskal.test(",
    ending = ")"
    ),
  arguments = list(
    )
  )

  
##################################################
## add this for student convenience
summarized.t.test = function(xbar, sx, nx,
  ybar = NULL, sy=NULL, ny=NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0, var.equal = FALSE,
  conf.level = 0.95) {
  
  paired = FALSE
  y = ybar
  
  alternative <- match.arg(alternative)
  
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")

    dname = "Summarized data"
  
  ##  if( !is.null(y) ) {
    ## 	dname <- paste(deparse(substitute(x)),"and",
    ## 		       deparse(substitute(y)))
    ## 	if(paired)
    ## 	    xok <- yok <- complete.cases(x,y)
    ## 	else {
    ## 	    yok <- !is.na(y)
    ## 	    xok <- !is.na(x)
    ## 	}
    ## 	y <- y[yok]
    ##     }
    ##     else {
    ## 	dname <- deparse(substitute(x))
    ## 	if( paired ) stop("'y' is missing for paired test")
    ## 	xok <- !is.na(x)
    ## 	yok <- NULL
    ##     }
    ##     x <- x[xok]
    ##     if( paired ) {
    ## 	x <- x-y
## 	y <- NULL
##     }
###    nx <- length(x)
    if(nx < 2) stop("not enough 'x' observations")
    mx <- xbar ## mean(x)
    vx <- sx^2 ## var(x)
    estimate <- mx
    if(is.null(y)) {
	df <- nx-1
	stderr <- sqrt(vx/nx)
        if(stderr < 10 *.Machine$double.eps * abs(mx))
            stop("data are essentially constant")
	tstat <- (mx-mu)/stderr
	method <- ifelse(paired,"Paired t-test","One Sample t-test")
	names(estimate) <- ifelse(paired,"mean of the differences","mean of x")

      } else {
##	ny <- length(y)
	if(ny < 2) stop("not enough 'y' observations")
	my <- ybar ##mean(y)
	vy <- sy^2 ##var(y)
	method <- paste(if(!var.equal)"Welch", "Two Sample t-test")
	estimate <- c(mx,my)
	names(estimate) <- c("mean of x","mean of y")
	if(var.equal) {
	    df <- nx+ny-2
	    v <- ((nx-1)*vx + (ny-1)*vy)/df
	    stderr <- sqrt(v*(1/nx+1/ny))
	} else {
	    stderrx <- sqrt(vx/nx)
	    stderry <- sqrt(vy/ny)
	    stderr <- sqrt(stderrx^2 + stderry^2)
	    df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
	}
        if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
	pval <- pt(tstat, df)
	cint <- c(-Inf, tstat + qt(conf.level, df) )
    }
    else if (alternative == "greater") {
	pval <- pt(tstat, df, lower.tail = FALSE)
	cint <- c(tstat - qt(conf.level, df), Inf)
    }
    else {
	pval <- 2 * pt(-abs(tstat), df)
	alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
	cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
	       conf.int=cint, estimate=estimate, null.value = mu,
	       alternative=alternative,
	       method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
  }

## make a generic widget for this

## t.test
t.test.summaries.list = list(
  title = "t.test() (summarized)",
  help = "t.test",
  type = "text",                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "summarized.t.test(",
    ending = ")"
    ),
  arguments = list(
    data = list(
      xbar=list(type="gedit",text=""),
      ybar=list(type="gedit",text=""),
      sx=list(type="gedit",text=""),
      sy=list(type="gedit",text=""),
      nx=list(type="gedit",text=""),
      ny=list(type="gedit",text="")
      ),
    hypotheses = list(
      mu = list(
        type = "gedit",
        text = "0"
        ),
      alternative = alternative.list,
#      paired = FALSE.list,
      var.equal = FALSE.list
      ),
    "CI" = list(
      conf.level = conf.level.list
      )
    )
  )

