cummeta.summaries<-function (effects,stderrs, conf.level = 0.95,
                             names =NULL, weights=NULL,
                             data = NULL, subset = NULL,
                             na.action = get(getOption("na.action")),
                             method = c("fixed","random"),logscale=TRUE)
{
    if (conf.level > 1 & conf.level < 100)
        conf.level <- conf.level/100
    if (is.null(data))
        data <- sys.frame(sys.parent())
    mf <- match.call()
    mf$data <- NULL
    mf$na.action <- NULL
    mf$subset <- NULL
    mf$conf.level <- NULL
    mf[[1]] <- as.name("data.frame")
    mf <- eval(mf, data)
    if (!is.null(subset))
        mf <- mf[subset, ]
    m <- NROW(mf)
    method <- match.arg(method)
    fn <- meta.summaries
    if (method=="fixed") {
        vals <- lapply(1:m, function(i) fn(effects,stderrs, method=method,
                                       data = mf, subset = 1:i, conf.level = conf.level,
                                           logscale=logscale))
    } else {
        vals <- c(list(fn(effects,stderrs, method="fixed",
                      data = mf, subset = 1, conf.level = conf.level,
                     logscale=logscale)),
                  lapply(2:m, function(i) fn(effects,stderrs, method=method,
                                             data = mf, subset = 1:i, conf.level = conf.level,
                                             logscale=logscale))
                  )
    }

    results <- sapply(vals, function(v) c(v$summary, v$se.summary))
    
    rval <- list(summary = results[1, ], se = results[2, ], names =
                 as.character(mf$names),
                 method = method, method=method, total = results[[m]],
                 logeffect=logscale,
                 call = match.call(), conf.level = conf.level)
    class(rval) <- "meta.cum"
    rval
}


cummeta<- function ( ntrt, nctrl, ptrt, pctrl, conf.level = .95,
		      names = NULL, data = NULL, 
		      subset = NULL, na.action = na.fail,
			method=c("meta.MH","meta.DSL"),statistic="OR") {
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$na.action<-NULL
    mf$subset <- NULL
    mf$statistic<- NULL
    mf$conf.level<-NULL
    mf$method<-NULL 
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    m<-NROW(mf)

    method<-match.arg(method)     
    fn<-get(method)
    vals<-lapply(1:m,function(i) fn(ntrt,nctrl,ptrt,pctrl,data=mf,subset=1:i,
				conf.level=conf.level,statistic=statistic))
    if (method=="meta.DSL")
        results<-sapply(vals,function(v) c(v$logDSL,v$selogDSL))
    else
	results<-sapply(vals,function(v) c(v$logMH,v$selogMH))
    	rval<-list(summary=results[1,],se=results[2,],names=as.character(mf$names),
		method=method,statistic=statistic,
		total=results[[m]],call=match.call(),conf.level=conf.level)
	class(rval)<-"meta.cum"
	rval
 }

print.meta.cum<-function(x,...){
	cat("Cumulative meta-analysis:\n")
	print(x$call)
}

summary.meta.cum<-function(object ,conf.level=NULL,...) {
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    
    m<-outer(object$se, c(0, -ci.value, ci.value), "*")
    m<-m+cbind(object$summary, object$summary, object$summary)
    
    if (is.null(object$statistic)){
        if (object$logeffect)
            m<-exp(m)
        dimnames( m ) <- list( as.character( object$names ),
                              c( "Effect", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
        
    } else if (object$statistic=="OR"){
        m <- exp(m)
        dimnames( m ) <- list( as.character( object$names ),
                              c( "OR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    } else {
	m <- exp(m) 
        dimnames( m ) <- list( as.character( object$names ),
                              c( "RR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    }
    
    rval <- list( stats = m, call = object$call, 
                 conf.level=conf.level, statistic=object$statistic
                 )
    class( rval ) <- "summary.meta.cum"
    rval 
}

print.summary.meta.cum <- function( x, ... ) {
    cat( "Cumulative meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    cat( "------------------------------------\n" )
    print( round( x$stats,2 ) )
    cat( "------------------------------------\n" )
}

plot.meta.cum <- function( x,  conf.level=NULL,colors=meta.colors(),xlab=NULL,
			summary.line=TRUE,summary.conf=FALSE,main="Cumulative meta-analysis",lwd=1,... ){
    if (is.null(conf.level))
        conf.level <- x$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    stats<-x$summary
    ses<-x$se

    if (is.null(x$statistic)){
        logeffect<-x$logeffect
        if (is.null(xlab))
            xlab<-"Effect"
    } else if (x$statistic=="OR") {
        if (is.null(xlab)) xlab<-"Odds Ratio"
        logeffect<-TRUE
    } else if (x$statistic=="RR"){
	if (is.null(xlab)) xlab<-"Relative Risk"
        logeffect<-TRUE
    }
    
    m<-length(stats)
    metaplot( stats, ses, labels=x$names, logeffect = logeffect,colors=colors,xlab=xlab,lwd=lwd,... )
    title(main=main)
    if (summary.line) abline(v=exp(stats[m]),lty=3,col=colors$summary,lwd=lwd)
    if (summary.conf) {
        z<-qnorm((1-conf.level)/2,s=ses[m])
        abline(v=exp(stats[m]+z),lty=3,col=colors$line,lwd=lwd)
        abline(v=exp(stats[m]-z),lty=3,col=colors$line,lwd=lwd)
    }
    
}
 


meta.MH <- function ( ntrt, nctrl, ptrt, pctrl, conf.level = .95,
		      names = NULL, data = NULL, 
		      subset = NULL, na.action = na.fail ,statistic="OR") {
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$na.action<-NULL
    mf$subset <- NULL
    mf$statistic<- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
    A <- as.numeric(mf$ptrt) ## integers would overflow
    B <- as.numeric(mf$pctrl)
    C <- as.numeric(mf$ntrt - mf$ptrt)
    D <- as.numeric(mf$nctrl - mf$pctrl)
    logORs <- log( ( A * D ) / ( B * C ) )
    logRRs <- log( A/(A+C))-log(B/(B+D))
    varORs <- 1/A + 1/B + 1/C + 1/D
    varRRs<- C/(A*(A+C))+D/(B*(B+D))
    Ti <- A + B + C + D


    P <- ( A + D ) / Ti
    Q <- ( B + C ) / Ti
    R <- A * D / Ti
    S <- B * C / Ti
   
    G<- A*(B+D)/Ti
    H<- B*(A+C)/Ti

    logMH <- log( sum( R ) / sum( S ) )
    varMH <- sum( P * R ) / ( 2 * sum( R )^2 ) + 
    	     sum( P * S + Q * R ) / ( 2 * sum( R ) * sum( S ) ) + 
    	     sum( Q * S ) / ( 2 * sum( S ) ^ 2 )
    MHchisq <- ( sum( A - ( A + B ) * ( A + C ) / Ti )^2 ) /
    		 sum( ( A + B ) * ( C + D ) * ( A + C ) * ( B + D ) /
    		      ( Ti * Ti * ( Ti - 1 ) 
    		    ) 
    	       )

    logMHRR<-log(sum(G)/sum(H))
    varMHRR <- (sum(((A + B) * (A + C) * (B + D) - A * B *Ti)/Ti^2)/(sum(A * (B + D)/Ti) * sum(B * (A + C)/Ti)))
    
    ok <- is.finite( logORs )
    okRR <- is.finite(logRRs) 
    heterog <- sum( ( ( logORs[ok] - logMH ) ^ 2 ) / varORs[ok] )
    heterogRR <- sum( ( ( logRRs[ok] - logMHRR ) ^ 2 ) / varRRs[ok] )
    df <- sum( ok ) - 1
    vars<-switch(statistic, OR=varORs, RR=varRRs)
    if (statistic=="OR") { 
         rval <- list( logOR = logORs, selogOR = sqrt( vars ), 
    		  logMH = logMH, selogMH = sqrt( varMH ), 
    		  MHtest = c( MHchisq, 1 - pchisq( MHchisq, 1 ) ), 
    		  het = c( heterog, df, 1 - pchisq( heterog, df ) ),
    		  call = match.call(), names = as.character(mf$names),
    		  conf.level = conf.level, statistic="OR" )
 	class( rval ) <- c("meta.MH.OR","meta.MH")

   } else {
        rval <- list( logRR = logRRs, selogRR = sqrt( vars ), 
    		  logMH = logMHRR, selogMH = sqrt( varMHRR ), 
    		  MHtest = c( MHchisq, 1 - pchisq( MHchisq, 1 ) ), 
    		  het = c( heterogRR, df, 1 - pchisq( heterogRR, df ) ),
    		  call = match.call(), names = as.character(mf$names),
    		  conf.level = conf.level, statistic="RR" )
 	class( rval ) <- c("meta.MH.RR","meta.MH")

 }
       rval
}

# Print out the output of meta.MH function.
print.meta.MH <- function( x, ... ) {
    cat( "Fixed effects ( Mantel-Haenszel ) Meta-Analysis\n" )
    cat( "Call: " )
    print( x$call )
    conf.level <- x$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    ci <- exp( x$logMH + c( -ci.value, 0, ci.value ) * x$selogMH )
    label<-paste("Mantel-Haenszel",x$statistic,"=")
    cat( paste( label, round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ", round( ci[1],2 ), 
    	        ", ", round( ci[3],2 ), " )\n", sep = "" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
    	        round( x$het[1],2 ), " ( p-value ", 
    	        round( x$het[3],4 ), " )\n", sep="" ) )
} 

# Summarise meta.MH function, same with print.meta.MH, but with odds
# ratios for each variable.

summary.meta.MH <- function( object ,conf.level=NULL,...) {
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
   if (object$statistic=="OR"){
    m <- exp( outer( object$selogOR, c( 0, -ci.value, ci.value ), "*" ) +
              cbind( object$logOR, object$logOR, object$logOR ) 
            )
    dimnames( m ) <- list( as.character( object$names ),
                          c( "OR", "(lower ", paste(100*conf.level,"% upper)",
                                                    sep=""))  )
    } else {
	m <- exp( outer( object$selogRR, c( 0, -ci.value, ci.value ), "*" ) +
           	  cbind( object$logRR, object$logRR, object$logRR ) 
        	    )
   	 dimnames( m ) <- list( as.character( object$names ),
                          c( "RR", "(lower ", paste(100*conf.level,"% upper)",
                                                    sep=""))  )
   }
    rval <- list( stats = m, call = object$call, 
    		  MHci = exp( object$logMH + c( -ci.value, 0, ci.value ) *
    		              object$selogMH 
    			    ),
    		  het = object$het,
                  conf.level=conf.level, statistic=object$statistic
    	        )
    class( rval ) <- "summary.meta.MH"
    rval
}

# Print out the summarised meta.MH function.
print.summary.meta.MH <- function( x, ... ) {
    cat( "Fixed effects ( Mantel-Haenszel ) meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    cat( "------------------------------------\n" )
    print( round( x$stats,2 ) )
    cat( "------------------------------------\n" )
    conf.level <- x$conf.level
    label<-paste("Mantel-Haenszel",x$statistic,"=")
    cat( paste( label, round( x$MHci[2], 2 ), " ",
    	        conf.level*100, "% CI ( ", round( x$MHci[1],2 ), ",", 
    	        round( x$MHci[3],2 ), " )\n", sep="" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
       	        round( x$het[1],2 ), " ( p-value ",
       	        round( x$het[3],4 ), " )\n", sep = "" 
       	      ) 
       )
}


# Plot the Odds Ratio of meta.MH
plot.meta.MH <- function( x, summary = TRUE, summlabel = "Summary", conf.level=NULL,colors=meta.colors(),xlab=NULL,... ){
    if (is.null(conf.level))
        conf.level <- x$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if (x$statistic=="OR") {
	stats<-x$logOR
	ses<-x$selogOR
	if (is.null(xlab)) xlab<-"Odds Ratio"
    } else {
        stats<-x$logRR
	ses<-x$selogRR
	if (is.null(xlab)) xlab<-"Relative Risk"
    }

    if ( summary )
      metaplot( stats,ses, labels=x$names,
               summn = x$logMH, sumse = x$selogMH, 
               conf.level = conf.level, sumnn = x$selogMH ^ ( -2 ),
               summlabel = summlabel, logeffect = TRUE,colors=colors,
               xlab=xlab,... )
    else 
      metaplot(stats,ses, labels=x$names, logeffect = TRUE,
               colors=colors,xlab=xlab,... )
  }
 
# Meta Analysis using DerSimonian-Laird Method.
meta.DSL <- function( ntrt, nctrl, ptrt, pctrl, conf.level = .95,
                     names = NULL, data = NULL,
                     subset = NULL, na.action = na.fail,statistic="OR") {
  if (conf.level>1 & conf.level<100)
    conf.level<-conf.level/100
  if ( is.null( data ) ) 
    data <- sys.frame( sys.parent() )
  mf <- match.call()
  mf$data <- NULL 
  mf$subset <- NULL
  mf$na.action<-NULL
  mf$statistic<-NULL
  mf[[1]] <- as.name( "data.frame" )
  mf <- eval( mf,data )
  if ( !is.null( subset ) ) 
        mf <- mf[subset,]
  mf <- na.action( mf )
  A <- as.numeric(mf$ptrt) ## integers might overflow
  B <- as.numeric(mf$pctrl)
  C <- as.numeric(mf$ntrt - mf$ptrt)
  D <- as.numeric(mf$nctrl - mf$pctrl)
  logORs <- log( ( A * D ) / ( B * C ) )
  logRRs <- log( A/(A+C))-log(B/(B+D))
  varORs <- 1/A + 1/B + 1/C + 1/D
  varRRs<- C/(A*(A+C))+D/(B*(B+D))
  if(statistic=="OR"){
    ok <- is.finite( logORs )
    logs<-logORs
    vars<-varORs
  } else    {
    ok <- is.finite(logRRs) 
    vars<-varRRs
    logs<-logRRs
  }
  df <- sum( ok ) - 1
  if ( any( !ok ) )
    warning( "Studies with 0/Inf statistic omitted" )
  ## Nick Barrowman points out this is wrong (or at least non-standard)
  ## tau2 <- max( 0, var( logs[ok] ) - sum( vars[ok] ) / ( length( ok ) ) )
  ## wts <- 1 / ( vars + tau2 )
  
  vwts<-1/vars
  logpooled<-sum(vwts[ok]*logs[ok])/sum(vwts[ok])
  heterog <- sum( ( ( logs[ok] - logpooled )^2 ) / vars[ok] )    
  tau2<- max(0, (heterog-df)/(sum(vwts[ok])-sum(vwts[ok]^2)/sum(vwts[ok])))
  
  wts<-1/(vars+tau2)
  logDSL <- sum( wts[ok] * logs[ok] ) / sum( wts[ok] )
  varDSL <- 1 / sum( wts[ok] )
  
  summary.test <- logDSL / sqrt( varDSL )
  
  
  rval <- list( logs = logs, selogs = sqrt( vars ), 
               logDSL = logDSL, selogDSL = sqrt( varDSL ), 
               test = c( summary.test,
                 1 - pchisq( summary.test ^ 2, 1 ) ), 
               het = c( heterog, df, 1 - pchisq( heterog, df ) ),
               call = match.call(),
               names = as.character(mf$names), conf.level = conf.level,
               omitted = !ok, tau2=tau2, statistic=statistic)
  class( rval ) <- "meta.DSL"
  rval
}

# Print out the output of meta.DSL function.
print.meta.DSL <- function( x, ... ){
    conf.level <- x$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    ci <- exp( x$logDSL + c( -ci.value, 0, ci.value ) * x$selogDSL )
    cat( paste( "Summary ",x$statistic,"= ", round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ",round( ci[1],2 ), ", ",
    	        round( ci[3],2 )," )\n", sep="" ) )
    cat( paste( "Estimated random effects variance:",round( x$tau2,2 ),"\n" ) )
} 

# Summarise meta.DSL function, same with print.meta.DSL, but with odds
# ratios for each variable.
summary.meta.DSL <- function( object,conf.level=NULL, ... ){
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- exp( outer( object$selogs, c( 0,-ci.value,ci.value ), "*" ) +
    		     cbind( object$logs, object$logs, object$logs ) )
    dimnames( m ) <- list( as.character( object$names ),
                          c( object$statistic, "(lower ",
                            paste(100*conf.level,"% upper)",sep=""))  )
    rval <- list( ors = m, call = object$call, 
                 ci = exp( object$logDSL + c( -ci.value, 0, ci.value ) *
                   object$selogDSL ),
                 het = object$het, 
                 omitted = as.character( object$names )[object$omitted],
                 conf.level=conf.level,
                 tau2 = object$tau2,statistic=object$statistic )
    class( rval ) <- "summary.meta.DSL"
    rval
  }

## Print out the summarised meta.MH function.
print.summary.meta.DSL <- function( x, ... ) {
  conf.level <- x$conf.level
  cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
  cat( "Call: " )
  print( x$call )
  cat( "------------------------------------\n" )
  print( round( x$ors,2 ) )
  cat( "------------------------------------\n" )
  cat( paste( "Summary",x$statistic,"= ",round( x$ci[2],2 ), "  ",
             conf.level * 100, "% CI ( ",round( x$ci[1],2 ), ",",
             round( x$ci[3],2 )," )\n",sep="" ) )
  cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
             round( x$het[1],2 ),
             " ( p-value ", round( x$het[3],4 ), " )\n", sep="" ) )
  cat( paste( "Estimated random effects variance:",
             round( x$tau2,2 ), "\n" ) )
  if ( length( x$omitted )>0 ){
    cat( paste( "( ", length( x$omitted ),
               "studies with zero or infinite odds ratio omitted )\n" ) )
  }
}

# Plot the Odds Ratio of meta.DSL
plot.meta.DSL <- function( x,summary=TRUE,summlabel="Summary",
                          conf.level=NULL,colors=meta.colors(),xlab=NULL,... ){
  if (is.null(conf.level))
    conf.level <- x$conf.level
  if (conf.level>1 & conf.level<100)
    conf.level<-conf.level/100
  if (is.null(xlab)){
    if (x$statistic=="OR") 
      xlab<-"Odds Ratio"
    else
      xlab<-"Relative Risk"
  }
  if ( summary )
    metaplot( x$logs, x$selogs, conf.level = conf.level,
             labels = x$names, summn = x$logDSL,
             sumse = x$selogDSL, sumnn = x$selogDSL^( -2 ),
             summlabel = summlabel, logeffect = TRUE,colors=colors,
             xlab=xlab,... )
  else 
    metaplot( x$logs, x$selogs, labels = x$names,
             logeffect = TRUE,colors=colors,conf.level=conf.level,
             xlab=xlab,... )
}


meta.summaries<-function(d,se,method=c("fixed","random"),
                         weights=NULL,logscale=FALSE,names=NULL,data=NULL,
                         conf.level=.95,subset=NULL,na.action=na.fail)
{
  if (conf.level>1 & conf.level<100)
    conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$method<-NULL  
    mf$logscale<-NULL
    mf$subset <- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
  
    if (is.null(mf$names)){
	if (is.null(mf$d) || is.null(names(mf$d)))
	   mf$names<-seq(along=mf$d)
	else
	  mf$names<-names(mf$d)
    }
  mf$names<-as.character(mf$names)
  method<-match.arg(method)
  vars<-mf$se^2
  
  vwts<-1/vars
  fixedsumm<-sum(vwts*mf$d)/sum(vwts)
  Q<-sum( ( ( mf$d - fixedsumm )^2 ) / vars ) 
  df<-NROW(mf)-1
  
  tau2<-max(0, (Q-df)/(sum(vwts)-sum(vwts^2)/sum(vwts)))
  
  if(is.null(mf$weights)){
    if (method=="fixed"){
      wt<-1/vars
    } else {
      wt<-1/(vars+tau2)
    }
  } else
  wt<-mf$weights
  
  summ<-sum(wt*mf$d)/sum(wt)
  if (method=="fixed")
	varsum<-sum(wt*wt*vars)/(sum(wt)^2)
  else
	varsum<-sum(wt*wt*(vars+tau2))/(sum(wt)^2)
  
  summtest<-summ/sqrt(varsum)

  df<-length(vars)-1
  rval<-list(effects=mf$d, stderrs=mf$se, summary=summ,se.summary=sqrt(varsum),
	     test=c(summtest,1-pchisq(summtest^2,1)),
	     het=c(Q,df,1-pchisq(Q,df)),
	     call=match.call(), names=mf$names,tau2=tau2,
	     variance.method=method, weights=wt, 
	     weight.method=if(is.null(mf$weights)) method else "user",
	     conf.level=conf.level,logscale=logscale)
  class(rval)<-"meta.summaries"
  rval
}

summary.meta.summaries <- function( object ,conf.level=NULL, ...) {
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- outer( object$stderrs, c( 0, -ci.value, ci.value ), "*" ) +
              cbind( object$effects, object$effects, object$effects )
    label<-"Effect"
    if (object$logscale){
        m<-exp(m)
        label<-"exp(Effect)"
    }
    wt<-object$weights
    m<-cbind(m, round(wt*length(wt)/sum(wt),1))
    dimnames( m ) <- list( as.character( object$names ),
                          c( label, "(lower ", paste(100*conf.level,"% upper)",
                                                     sep=""),"weights")  )
    summci<-object$summary + c( -ci.value, 0, ci.value ) *object$se.summary
    if (object$logscale) summci<-exp(summci)
    rval <- list( studies = m, call = object$call, 
                 summci = summci,
                 het = object$het,tau2=object$tau2,
                 conf.level=conf.level,logscale=object$logscale,
                 weight.method=object$weight.method,
                 variance.method=object$variance.method
                 )
    class( rval ) <- "summary.meta.summaries"
    rval
  }

## Print out the summarised meta.summaries function.
print.summary.meta.summaries <- function( x, ... ) {
  if (x$weight.method=="user")
    cat(paste("Weighted meta-analysis with ",x$variance.method,
              "-effects standard errors\n",sep=""))
  else
    cat(paste(switch(x$variance.method,fixed="Fixed",random="Random"),
              "-effects meta-analysis\n",sep=""))
  
  cat( "Call: " )
  print( x$call )
  cat( "----------------------------------------------------\n" )
  print( round( x$studies,2 ) )
  cat( "----------------------------------------------------\n" )
  conf.level <- x$conf.level
  if (x$logscale)
    cat("Summary exp(effect): ")
  else
    cat("Summary effect: ")
  cat( paste(  round( x$summci[2], 2 ), " ",
             conf.level*100, "% CI ( ", round( x$summci[1],2 ), ",", 
             round( x$summci[3],2 ), " )\n", sep="" 
             ) 
      )
  if (x$variance.method=="fixed")
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
               round( x$het[1],2 ), " ( p-value ",
               round( x$het[3],4 ), " )\n", sep = "" 
               ) 
        )
  else
    cat("Estimated heterogeneity variance:",signif(x$tau2,2)," p=",round(x$het[3],3),"\n")
}

print.meta.summaries<-function(x, ...){
  conf.level<-x$conf.level
  ci.value<- -qnorm((1-conf.level)/2)
  if (x$variance.method=="fixed")
    cat("Fixed-effects meta-analysis")
  else	
    cat("Random-effects meta-analysis")
  if(x$weight.method!="user")
    cat("\n")
  else
    cat(" with user-supplied weights\n")
  cat("Call: ")
  print(x$call)
  ci<-x$summary+c(-ci.value,0,ci.value)*x$se.summary
  if (x$logscale) {
    ci<-exp(ci)
    label<-"exp(summary effect)="
  }
  else
    label<-"Summary effect="
  cat(paste(label,signif(ci[2],3), "   ",
            conf.level*100,"% CI (",signif(ci[1],3),", ",signif(ci[3],3),
	")\n",sep=""))
  cat("Estimated heterogeneity variance:",signif(x$tau2,2)," p=",round(x$het[3],3),"\n")
}

plot.meta.summaries<-function(x,summary=TRUE,summlabel="Summary",
                              conf.level=NULL,colors=meta.colors(),
                              xlab=NULL,logscale=NULL,... )
{
  if (is.null(conf.level))
    conf.level <- x$conf.level
  if (conf.level>1 & conf.level<100)
    conf.level<-conf.level/100
  if(is.null(logscale)) logscale<-x$logscale
  if (is.null(xlab)){
    if (logscale) xlab<-"exp(Effect)" else xlab<-"Effect"
  }
  if ( summary )
    metaplot( x$effects, x$stderrs, conf.level = conf.level,
             labels = x$names, summn = x$summary,
             sumse = x$se.summary, sumnn = sum(x$weights),
             summlabel = summlabel, logeffect = logscale,
             colors=colors,nn=x$weights,xlab,... )
  else 
    metaplot( x$effects, x$stderrs, labels = x$names,
             logeffect = logscale,colors=colors,
             nn=sum(x$weights),conf.level=conf.level,xlab,... )
  
}



funnelplot<-function(x,...)
    UseMethod("funnelplot")

funnelplot.meta.MH<-function(x,...){
    funnelplot.default(x$logOR,x$selogOR,summ=x$logMH,...)
}
funnelplot.meta.DSL<-function(x,...){
    funnelplot.default(x$logOR,x$selogOR,summ=x$logDSL,...)
}

funnelplot.meta.summaries<-function(x,...){
    funnelplot.default(x$effects,x$stderrs,summ=x$summary,...)
}

