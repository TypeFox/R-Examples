poisson.gui <- function(base.txt) {
 
    # fits a poisson distribution with
    # a trend if desired

    diags.value <- tclVar(0)
    
    ########################################
    #  Internal functions
    ########################################
    
    submit <- function() {
    
    	# Obtain selected data.
    	data.select <- as.numeric( tkcurselection( data.listbox))+1
            if( is.nothing) dd.cmd <- "dd <- in2extRemesData"
            else dd.cmd <- paste( "dd <- get( \"", full.list[ data.select], "\")", sep="")
    	eval( parse( text=dd.cmd))
    	write( dd.cmd, file="in2extRemes.log", append=TRUE)
    
        # fit the poisson distribution
    
        resp.select<-as.numeric(tkcurselection(resp.listbox))+1
    
        if (is.na(resp.select))
          return()
    
    	cat( "\n", "\n", "Poisson Fit\n")
    	cat( "-------------------------\n", "\n")
    
    
        if (tclvalue( tkcurselection(trend.list)) =="") {

	    cmd <- paste("print(fpois(dd[[\"data\"]], which.col = ", resp.select, "))", sep = "")
	    eval(parse(text = cmd))
	    write(cmd, file = "in2extRemes.log", append = TRUE)

            # no trend component indicated
        
    	    # gsum.cmd <- paste( "gsum <- sum( dd[[\"data\"]][, ", resp.select, "])", sep="")
    	    # eval( parse( text=gsum.cmd))
    	    # write( gsum.cmd, file="in2extRemes.log", append=TRUE)
    
    	    # m1.cmd <- paste( "m1 <- length( dd[[\"data\"]][, ", resp.select, "])", sep="")
    	    # eval( parse( text=m1.cmd))
    	    # write( m1.cmd, file="in2extRemes.log", append=TRUE)

	    # lambda.cmd <- "lambda <- gsum/m1"
	    # eval( parse( text=lambda.cmd))
	    # write( lambda.cmd, file="in2extRemes.log", append=TRUE)
    
    	    # sigma2.cmd <- paste( "sigma2 <- var( dd[[\"data\"]][, ", resp.select, "], na.rm=TRUE)", sep="")
    	    # eval( parse( text=sigma2.cmd))
    	    # write( sigma2.cmd, file="in2extRemes.log", append=TRUE)

    	    # chisquaretest.cmd <- "chisqtest <- (m1-1)*sigma2/lambda"
    	    # eval( parse( text=chisquaretest.cmd))
    	    # write( chisquaretest.cmd, file="in2extRemes.log", append=TRUE)
    
    	    # pvalCMD <- "pval <- pchisq( chisqtest, df=m1-1, lower.tail=FALSE)"
    	    # eval( parse( text=pvalCMD))
    	    # write( pvalCMD, file="in2extRemes.log", append=TRUE)

    	    # nl1 <- paste("**********", " ", sep="\n")
    	    # cat( nl1)
    	    # nl2 <- paste(" ", " ", sep="\n")	

            # mess <- paste("  Lambda: ",round(lambda,3), sep="")
    	    # cat( nl2)
    	    # cat( mess)

    	    # mess2 <- paste(" Variance: ", round( sigma2,3), sep="")
    	    # cat( mess2)

    	    # cat( nl2)

    	    # mess3 <- paste(" Chi-square statistic: ", round( chisqtest, 3), " with ", m1-1, " degrees of freedom",
	# 		sep="")
    	 #    cat(mess3)
    	  #   cat(nl2)
    	   #  mess4 <- paste(" p-value for Chi-square test of equality of mean and variance: ", round( pval, 3),
	# 		sep="")
    	 #    cat( mess4)
    	  #   cat( nl2)
    	   #  cat( nl1)

        } else {

	    # trend component indicated
    
	    trend.select <- as.numeric( tkcurselection( trend.list))+1 
	    mess <- paste("Trend variable: ", colnames( dd$data)[trend.select],"\n")
    	    cat( mess)

            mess <- paste( "Response variable: ", colnames(dd$data)[resp.select],"\n\n")
            cat( mess)

    	    response.cmd <- paste( "response <- dd[[\"data\"]][, ", resp.select, "]", sep="")
    	    eval( parse( text=response.cmd))
    	    write( response.cmd, file="in2extRemes.log", append=TRUE)

    	    cnamesCMD <- paste( "cnames <- colnames( dd[[\"data\"]])[-", resp.select, "]", sep="")
    	    eval( parse( text=cnamesCMD))
    	    write( cnamesCMD, file="in2extRemes.log", append=TRUE)

    	    trendNameCMD <- "trendName <- character(0)"
    	    eval( parse( text=trendNameCMD))
    	    write( trendNameCMD, file="in2extRemes.log", append=TRUE)

    	    for( i in 1:length( trend.select)) {

    	 	trendNameCMD <- paste( "trendName <- c( trendName, \"", cnames[trend.select[i]], "\")", sep="")
    		eval( parse( text=trendNameCMD))
    	 	write( trendNameCMD, file="in2extRemes.log", append=TRUE)

    	    } # end of for 'i' loop.

    	    trendCMD <- "trend <- dd[[\"data\"]][, trendName]"
    	    eval( parse( text=trendCMD))
    	    write( trendCMD, file="in2extRemes.log", append=TRUE)
    
    	    # Put fit in 'models' component of list.
            number.of.models <- length( dd$models)
            names.of.models <- names( dd$models)
            if( is.null( names.of.models)) names.of.models <- character(0)
            jj <- 0
            if(number.of.models > 0) for(i in 1:number.of.models) if(class(dd$models[[i]])[1] == "glm") jj <- jj+1
            names.of.models <- c( names.of.models, paste( "poisson.fit", jj+1, sep=""))
    
            # fit the Poisson distribution.
    	    mod.fit.cmd <- "mod.fit <- glm( response~trend, family=poisson())"
    	    eval( parse( text=mod.fit.cmd))
    	    write( mod.fit.cmd, file="in2extRemes.log", append=TRUE)
    
    	    # Put the fit in the 'models' component of 'in2extRemesDataObject' object.
    	    putfitCMD <- paste( "dd[[\"models\"]][[\"poisson.fit", jj+1, "\"]] <- mod.fit", sep="")
    	    eval( parse( text=putfitCMD))
    	    write( putfitCMD, file="in2extRemes.log", append=TRUE)
    
    	    if( tclvalue( diags.value) == 1) {

    	        plotCMD <- "plot( mod.fit[[\"residuals\"]], type=\"l\")"
    	        eval( parse( text=plotCMD))
    	        write( plotCMD, file="in2extRemes.log", append=TRUE)
    	        ablineCMD <- "abline(h=0, lty=2)"
    	        eval( parse( text=ablineCMD))
    	        write( ablineCMD, file="in2extRemes.log", append=TRUE)

	    }

    	    nl1 <- paste(" ", "**********", " ", sep="\n")
    	    nl2 <- paste( " ", " ", sep="\n")
            no.of.pars <- length( mod.fit$coef)
            summ.fit <- summary( mod.fit)$coef
            coef.names <- c("Intercept", cnames[trend.select])
            nchar1 <- nchar( coef.names[2])
            if( nchar1 < 9) for( j in 1:(9-nchar1)) coef.names[2] <- paste( coef.names[2], ".", sep="")
            else if( nchar1 > 9) for( j in 1:(nchar1-9)) coef.names[1] <- paste( coef.names[1], ".", sep="")
    
            summaryCMD <- "print( summary( mod.fit))"
            eval( parse( text=summaryCMD))
            write( summaryCMD, file="in2extRemes.log", append=TRUE)
    
            msg3 <- paste( "Residual degrees of freedom: ", mod.fit$df.residual, sep="")
            lr <- mod.fit$null.deviance - mod.fit$deviance
            pval <- pchisq( lr, 1, lower.tail=FALSE)
            msg4 <- paste( "negative log-likelihood (proportional): ",
    		round( mod.fit$deviance, digits=4), sep="")
            msg5 <- paste( "likelihood-ratio (against null model): ", round( lr, digits=4), 
    		" (p-value = ", round( pval, digits=4), ")", sep="")
    
            cat( nl2)
            cat( nl2)
            cat( msg5)
            cat( nl2)
            cat( nl1)
    	
        } # end of if else trend stmts.

    	tkyview.moveto(base.txt,1.0) 
        tkdestroy(base)
    	if( is.nothing) assignCMD <- "assign( \"in2extRemesData\", dd, pos=\".GlobalEnv\")"
        else assignCMD <- paste( "assign( \"", full.list[ data.select], "\", dd, pos=\".GlobalEnv\")", sep="") 
    	eval( parse( text=assignCMD))
    	write( assignCMD, file="in2extRemes.log", append=TRUE)

        invisible()

    } # end of internal 'submit' fcn
    
    refresh <- function() {

    	tkdelete( resp.listbox, 0.0, "end")
    	tkdelete( trend.list, 0.0, "end")
    	data.select <- as.numeric( tkcurselection( data.listbox))+1
    	if( is.nothing) dd <- in2extRemesData
    	else dd <- get( full.list[ data.select])
    	
    	for( i in colnames( dd$data)) {
    		tkinsert( resp.listbox, "end", i)
    		tkinsert( trend.list, "end", i)
    		} # end of for i loop

    	invisible()

    	} # end of refresh fcn
    
    redolists <- function() {
    
    data.select <- as.numeric( tkcurselection( data.listbox))+1
            if( is.nothing) dd <- in2extRemesData
            else dd <- get( full.list[ data.select])
    	
        resp.name<-colnames(dd$data)[as.numeric(tkcurselection(resp.listbox))+1]
     
        # put the correct eligible covariates in the other list boxes
        tkdelete(trend.list, 0.0,"end")
     
        for (i in colnames(dd$data)) {
          if (i != resp.name) {
            tkinsert(trend.list,"end",i)
          }
        }
    } # end of redolists fcn
    
    poissonhelp <- function() {
    	# tkconfigure( base.txt, state="normal")
    	cat("\n", "For fits without a trend, this function simply returns a mean, and\n")
    	cat( "no new objects (or components of objects) are returned.\n")
    	cat( "\n", "For fits with a trend, the R function \'glm\' is invoked.\n")
    	help.msg <- paste( "Use \'help( glm)\' for more help.", " ", sep="\n")
    	# tkinsert( base.txt, "end", help.msg)
    	cat( help.msg)
    	# tkconfigure( base.txt, state="disabled")
    
    	help( glm)
    	invisible()
    	} # end of poissonhelp fcn
    
    endprog <- function() {

        tkdestroy(base)
    
      }
    
    ##################################
    # Frame/ button setup
    ##################################
    
    
    base<-tktoplevel()
    tkwm.title(base,"Fit Poisson Distribution")
    
    data.frm <- tkframe( base, borderwidth=2, relief="groove") 
    top.frm <- tkframe( base, borderwidth=2, relief="groove")
    mid.frm <- tkframe( base, borderwidth=2, relief="groove")
    trend.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    diag.frm <- tkframe( mid.frm, borderwidth=2, relief="flat")
    bot.frm <- tkframe( base, borderwidth=2, relief="groove")
    
    # data frame to select a data object.
    
    data.listbox <- tklistbox( data.frm,
    			yscrollcommand=function(...) tkset( data.scroll, ...),
    			selectmode="single",
    			width=20,
    			height=5,
    			exportselection=0)
    
    data.scroll <- tkscrollbar( data.frm, orient="vert",
    			command=function(...) tkyview( data.listbox, ...))
    
    temp <- ls( all.names=TRUE, name=".GlobalEnv")
    full.list <- character(0)
    is.nothing <- TRUE
    for( i in 1:length( temp)) {
    	if( is.null( class( get( temp[i])))) next
    	if( (class( get( temp[i]))[1] == "in2extRemesDataObject")) {
    		tkinsert( data.listbox, "end", paste( temp[i]))
    		full.list <- c( full.list, temp[i])
    		is.nothing <- FALSE
    		} # end of if class stmt
    	} # end of for i loop
    
    tkpack( tklabel( data.frm, text="Data Object:  ", padx=4), side="left")
    tkpack( data.listbox, data.scroll, side="left", fill="y")
    tkpack( data.frm, fill="x")
    
    # place bindings on data.listbox in order to update response info.
    tkbind( data.listbox, "<Button-1>", "")
    tkbind( data.listbox, "<ButtonRelease-1>", refresh)
    
    # top frame for response variable
     
    top.l <- tkframe( top.frm,borderwidth=2)
    resp.listbox <- tklistbox( top.l,
    			yscrollcommand=function(...)tkset(resp.scroll,...),
    			selectmode="single",
    			width=35,
    			height=6,
    			exportselection=0)
    resp.scroll <- tkscrollbar( top.l,orient="vert",
    			command=function(...)tkyview(resp.listbox,...))
    
    if( is.nothing) { 
      for (i in 1:ncol(in2extRemesData$data)) {
        tkinsert(resp.listbox,"end",paste(colnames(in2extRemesData$data)[i]))
      }
    	} else tkinsert( resp.listbox, "end", "")
    
     
      tkpack(tklabel(top.l,text="Response:      ",padx=4), side="left")
      tkpack(resp.listbox, resp.scroll, side="left", fill="y")
    
      tkpack(top.l,side="left",fill="x") 
    
      # place binding on resp.listbox to eliminate the
      # response from the lists of covs
    tkbind( resp.listbox, "<Button-1>", "")
      tkbind(resp.listbox,"<ButtonRelease-1>",redolists)
    
      # choose the trend variable
      trend.l<-tkframe( trend.frm, borderwidth=2)
      trend.list<-tklistbox(trend.l,yscrollcommand=function(...)tkset(trend.covscr,...),selectmode="multiple",width=15,height=6,exportselection=0)
      trend.covscr <- tkscrollbar(trend.l,orient="vert",command=function(...)tkyview(trend.list,...))
    
    if( is.nothing) {
      for (i in 1:ncol(in2extRemesData$data)) {
        tkinsert(trend.list,"end",paste(colnames(in2extRemesData$data)[i]))
      }
    	} else tkinsert( trend.list, "end", "")
    
      tkpack(tklabel(trend.l,text="Covariate (log link):",padx=4), side="left")
      tkpack(trend.list,side="left")
      tkpack(trend.covscr,side="right",fill="y")
     
      tkpack(trend.l,side="right")
    
    # plot diagnostics frame
    diags.button <- tkcheckbutton( diag.frm, text="plot diagnostics",
    				variable=diags.value)
    tkpack( diags.button)
    tkpack( trend.frm, diag.frm, side="left")
    
      # bottom frame
    	ok.but <- tkbutton( bot.frm, text="OK", command=submit)
    	cancel.but <- tkbutton( bot.frm, text="Cancel", command=endprog)
    	help.but <- tkbutton( bot.frm, text="Help", command=poissonhelp)
     
    	tkpack( ok.but, cancel.but, side="left")
    	tkpack( help.but, side="right")
    
    # place bindings on buttons.
    	tkbind( ok.but, "<Return>", submit)
    	tkbind( cancel.but, "<Return>", endprog)
    	tkbind( help.but, "<Return>", poissonhelp)
    
      tkpack(top.frm, fill="x")
      tkpack(mid.frm)
      tkpack(bot.frm)

} # end of 'poisson.gui' function.
