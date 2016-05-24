# maRketSim - Object functions

# The heirarchy:
# core functions calculate values from values
# market objects calculates from core functions
# bond object and functions calculate from market object
# portfolio calculates from bond object and market object
# account calculates from bond object and market objects held within a portfolio


# --- Object fundamental functions (creation and display of objects) --- #

# - Bond object functions - #

print.bond <- function(x,...) {
	cat("Bond issued at time ",x[["t.issue"]]," with interest of ",x[["i"]]*100,"%, maturity of ",x[["mat"]]," years, price of $",x[["par"]],
		", coupon of $",round(x[["cpn"]],2)," paid every ",x[["f"]]," year.\n",sep="")
	invisible(x) # Return object invisibly
}

# Create a summary bond object which gets passed to print.sum.bond
	# Pass a market object to the function and take the time from there
summary.bond <- function(object,mkt,...) { # object is the bond object, market.rate = current market yield, time.elapsed = years since bond was issued
	# - Check parameters - #
	if(class(mkt)!="market") { stop("market must be a market object") }
	# - Get t and market.rate from the market object - #
	t = mkt$t
	market.rate <- market.rate(mkt,mat=object$mat,f=object$f)
	# - Update, extract, and calculate intermediate values - #
	time.elapsed = t - object$t.issue
	effective.mat <- object[["mat"]] - time.elapsed
	cpn <- object[["cpn"]]
	i <- market.rate*object[["f"]] # Periodic interest rate (market yield divided by number of payments per year)
	y <- market.rate*object[["f"]]
	n <- effective.mat / object[["f"]] # Number of periods
	# - Calculate final numbers - #
	pv <- cpn*(1-(1/(1+i)^n))/i + object[["par"]]/(1+i)^n
	# Duration calculation  #### change this over to use duration() - or should we, since this does fractional periods? or does it?
	periods <- seq(1,effective.mat/object[["f"]])
	pvcf <- cpn*1/(1+market.rate*object[["f"]])^periods #coupon present values of cash flow
	pvcf[length(pvcf)] <- pvcf[length(pvcf)]+object[["par"]]*1/(1+market.rate*object[["f"]])^periods[length(pvcf)] #now add in the PV of the maturity value to the final period
	dur.Macaulay <- object[["f"]] * (sum(periods * pvcf) / sum(pvcf) ) # Macaulay duration, in years
	dur <- dur.Macaulay / (1+y)
	cnvx <- object$f^2 * sum(periods*(periods+1) * pvcf) / ( (1+y)^2 * sum(pvcf) )
	# - Return results - #
	res <- list(mat=object[["mat"]],i=object[["i"]],market.rate=market.rate,effective.mat=effective.mat,par=object[["par"]],
		cpn=object[["cpn"]],f=object[["f"]],dur=dur,dur.Macaulay=dur.Macaulay,pv=pv,cnvx=cnvx)
	
	structure(res,class=c("sum.bond",class(object))) # Return the summary calculated values for the print.sum.bond function to display
}
       
# Print a sum.bond object
print.sum.bond <- function(x,...) {
	cat(x[["mat"]],"-year bond with interest of ",x[["i"]]*100,"% under a current market yield of ",x[["market.rate"]]*100,"%.\n",sep="")
	cat("Current maturity of ",x[["effective.mat"]]," years.\n",sep="")
	cat("Par value is $",x[["par"]],", paying $",round(x[["cpn"]],2)," every ",x[["f"]]," year.\n",sep="")
	cat("Present value (price) of the bond under current market conditions is $",round(x[["pv"]],2),"\n",sep="")
	cat("Macaulay duration is ",round(x$dur.Macaulay,2)," years, modified duration is ",round(x[["dur"]],2)," years, with a convexity of ",round(x$cnvx,2),".\n",sep="")
	invisible(x)
}

# Create an individual bond object which gets its issue date and interest rate from the market object
bond <- function(mkt,mat=NA,par=1000,f=.5,dur=NA) { # i = interest rate; mat = maturity, in years; par = par; f = frequency of payouts, in fractions of a year (.5 = semiannual)
	# Get market rate for a bond of this duration
	if(!is.na(dur)) i <- market.rate(mkt,dur=dur,f=f)
	if(!is.na(mat))	i <- market.rate(mkt,mat=mat,f=f)
	# Error check inputs
	if(is.na(mat)&is.na(dur)) stop("Must specify either maturity or duration.\n")
	if(class(mkt)!="market") stop("Market must be a market object.\n")
	if(f!=.5) stop("Only semi-annual (f=.5) bonds are supported for now.\n")
	if(i>=1) stop("Your interest rate should be expressed as a proportion (e.g. for a 5% yield, i=.05)\n")
	# Calculate coupon
	cpn <- i*par*f
	# Calculate maturity from duration rather than the (default) other way around, if specified
	if(!is.na(dur) & is.na(mat)) {
		mat <- findMat(dur=dur,i=i,f=f)
	}
	# Return results
	res <- list(i=i,mat=mat,par=par,f=f,cpn=cpn,t.issue=mkt$t)
	class(res) <- "bond"
	res
}

# - Cash object function - #
# compound should either be "continuous" or a fractional frequency of compounding (e.g. 0.5 for twice-yearly compounding)
cash <- function(value=0,mkt,compound="continuous",name=NA) {
	MMrate <- market.rate(mkt=mkt,mat=10^-5)
	if(is.na(name)) { name <- as.character(runif(1)) }
	res <- list(name=name,value=value,MMrate=MMrate,t.issue=mkt$t)
	class(res) <- "cash"
	res
}

print.cash <- function(x,...) {
	cat("Cash portfolio \'",x$name,"\' containing $",round(x$value,2)," created at time ",x$t.issue,".\n",sep="")
	invisible(x)
}


# - Portfolio object functions - #

# Print basic bond portfolio characteristics
print.portfolio.bond <- function(x,...) {
	num.bonds <- length(x$bonds)
	# Get the max number of digits for outputting par and coupon values (since it varies over such a wide range)
	digits.par <- nchar(format(max(sapply(x$bonds,function(x) { invisible(x$par) })),scientific=FALSE))
	digits.cpn <- nchar(format(max(sapply(x$bonds,function(x) { invisible(x$cpn) })),scientific=FALSE))
	cat("Portfolio \'",x$name,"\' created at time ",x$orig.mkt$t," containing ",num.bonds," bonds.\n",sep="")
	# Print table of bonds
	padding.par <- rep(" ",digits.par-2)
	padding.cpn <- rep(" ",ifelse(digits.cpn+1>6,digits.cpn-5,0))
	cat("i (%)    Maturity  ",padding.par,"Par   ",padding.cpn,"Coupon        t\n",sep="")
	for(n in seq(num.bonds)) {
		b <- x$bonds[[n]]
		cat(sprintf("%5.2f      %03.1f yr  $%*.0d   $%*.0d       %03.1f ",b$i*100,b$mat,digits.par,b$par,digits.cpn,as.integer(b$cpn),b$t),"\n")
	}
	invisible(x)
}

# Print bond portfolio details, potentially under new market conditions
print.sum.portfolio.bond <- function(x,...) {
	mkt <- x$portfolio.sum$mkt
	cat("Portfolio '",x$portfolio.bonds$name,"' of ",length(x$portfolio.bonds$bonds)," bonds created at time ",x$portfolio.bonds$orig.mkt$t,". Current time is ",mkt$t,".\n",sep="")
	cat("PV of portfolio is $",x$portfolio.sum$pv,"\n",sep="")
	cat("Duration of portfolio is ",x$portfolio.sum$dur,"\n",sep="")
	cat("Coupon yield is ",round(x$portfolio.sum$yield.current*100,2),"%\n",sep="")
	cat("Coupon total is $",x$portfolio.sum$annual.cpn," per year.\n",sep="")
	cat("N.B. Portfolio summaries do not include reinvested coupons or maturing securities.  For that, place the portfolio in an account.\n")
	invisible(x)
}

# Summarize portfolio.bond under future market conditions
summary.portfolio.bond <- function(object,mkt=object$orig.mkt, ...) {
	# Exclude bonds which have start t after current market time (e.g. which have not yet been purchased), and which have already matured
	object <- current.portfolio.bond(object,mkt)
	# Summarize remaining (purchased) bonds
	s <- lapply(object$bonds,summary,mkt=mkt)
	bonds.pv <- unlist(s)[names(unlist(s))=="pv"]
	pv <- sum(bonds.pv)
	bonds.dur <- unlist(s)[names(unlist(s))=="dur"]
	dur <- sum((bonds.pv/pv) * bonds.dur)
	bonds.cpn <- unlist(s)[names(unlist(s))=="cpn"]
	bonds.f <- unlist(s)[names(unlist(s))=="f"]
	annual.cpn <- sum(bonds.cpn/bonds.f)
	yield.current <- annual.cpn / pv
	res <- list(
		portfolio.sum=list(pv=pv,dur=dur,yield.current=yield.current,annual.cpn=annual.cpn,mkt=mkt),
		portfolio.bonds=object
		)
	structure(res,class=c("sum.portfolio.bond",class(object)))
}

# Create a bond portfolio object from a list of bonds
portfolio <- function(bonds,mkt,name=NA) {
	if(class(bonds)!="list" && class(bonds[[1]]!="bond")) { stop("Bonds must be a list of bond objects.\n") }
	if(is.na(name)) { name <- as.character(runif(1)) }
	port <- list(name=name,bonds=bonds,orig.mkt=mkt,t.issue=mkt$t)
	structure(port,class=c("portfolio.bond","bond"))
}


# - Market object functions - #

# -- Market container -- #
market <- function(mkts,t=0) {
	# Validation
	if(substr(class(mkts),1,6) == "market") { # If we've just been fed a raw market object, make it a list
		mkts <- list(mkts)
	}
	if(class(mkts) != "list" ) { stop("Markets list must be a list of market objects") }
	classes1 <- as.data.frame(sapply(mkts,class))[1,]
	if(any(table(classes1)>1)) { stop("Only one market per type!") }
	# Return object
	structure(list(mkts=mkts,t=t),class="market")
}

# Print market container
print.market <- function(x,...) {
	cat("Market time is:",x$t,"\n")
	sapply(x$mkts,"print")
}

# Summary market container
summary.market <- function(object,...) {
	sapply(object$mkts,"summary")
}

# -- Bond market -- #
# Create a market object which gives market (interest rate and yield curve characteristics, as well as the money-market rate) characteristics to bond objects
# if i is specified, it's a flat yield curve with interest rate i (e.g. i=0.05)
# if yield.curve is specified, it should be a quoted formula relating duration to i, e.g.: quote(mat/100+.01)
# MMrate = money market rate; MM.frequency = compounding frequency on cash
market.bond <- function(i=NA,yield.curve=NA,MMrate=NA,MM.frequency=.5) {
	# Error check inputs
	mat <- 5 # specify a maturity to test the yield curve equation
	if(is.na(i)&is.na(eval(yield.curve))) stop("Specify either a flat interest rate i or a yield.curve, not both")	
	if(!is.na(eval(yield.curve)) && (eval(yield.curve)>1 | eval(yield.curve) <0)) stop("Something's funny with your yield curve.")
	if(!is.na(i)&!is.na(eval(yield.curve))) stop("You've specified both i and yield.curve; one or the other please.\n")
	# Specify money-market rate if none is specified
	if(is.na(MMrate)) {
		if(!is.na(i)) {
			MMrate <- i
		} else {
			mat <- 1/12 # One-month maturity
			MMrate <- eval(yield.curve)
			cat("No money market rate specified; defaulting to the ",round(MMrate,3)*100,"% yield of a 1-month maturity issue.\n",sep="")
		}
	}
	# Convert i to a yield.curve if only i is specified (this causes problems!)
	#if(!is.na(i)&is.na(yield.curve)) {
	#	yield.curve=i
	#}
	# Build market object and return it                
	res <- list(i=i,yield.curve=yield.curve,MMrate=MMrate)
	class(res) <- "market.bond"
	invisible(res)
}

# Print market characteristics
print.market.bond <- function(x,...) {
	if(is.na(x$i)) {
		cat("Yield curve relating interest rate to maturity is: ")
		print(x$yield.curve)
		cat("Money market rate is:",x$MMrate,"\n")
	} else {
		cat("Yield curve is flat, with interest rate",x$i,"\nMoney market rate is:",x$MMrate,"\n")
	}
}

# Market summary
# i is the coupon rate of a non-mentioned bond you wish to analyze max duration under
summary.market.bond <- function( object, i=NA, ... ) {
  if(is.na(r)) { r <- object$i } #Assume at-par unless otherwise specified
	max.dur <- findMaxDur( market.rate=market.rate( mat= ), i=object$i )
  max.mat <- findMaxMat( r=r, i=object$i )
	# Output results
	res <- list(i = object$i, yield.curve = object$yield.curve, MMrate = object$MMrate, max.dur = max.dur, max.mat = max.mat, r=r )
	structure(res,class=c("sum.market",class(object)))
}
print.sum.market.bond <- function(x,...) {
	#cat("Maximum duration of a bond with coupon rate ",round(x$r,3)*100,"% is ",round(x$max.dur,2),"\n",sep="")
	invisible(x)
}

# - history.market functions - #
# Define a history.market option
history.market <- function(mkts,...) {
	# Validate input
	if(class(mkts)!="list" | any(sapply(mkts,"class")!="market")) { stop("mkts must be a list of market objects.\n") }
	# Sort
	times <- sapply(mkts,function(x) { invisible(x$t) } )
	mkts <- mkts[order(times)]
	# More input validation
	if( any(as.numeric(table(times))>1) ) { stop("You have more than one market for a given point in time.\n") }
	# Return
	structure(mkts,class=c("history.market","market"))
}

# Print history.market
print.history.market <- function(x,...) {
	cat("Market history\n--------------------------------------------\n")
	sapply(x,"print")
	invisible(x)
}

# Plot history.market
plot.history.market <- function(x,plot.MMrate=TRUE,plot.mats=c(1,2,5,10),cols=rainbow(length(plot.mats)),start.t=x[[1]]$t,end.t=x[[length(x)]]$t+0.5,xlab="Year",ylab="Interest Rate (%)",...) {
  require(taRifx)
	# Validate inputs
	if(length(x)<2) {
		cat("Cannot plot a market history of a single period.\n")
		return(x)
	}
	if(plot.MMrate) {
		plot.mats <- c(0.0001,plot.mats)
	}
	# Extract bond market history
	df <- data.frame(t=NA,i=NA,mat=NA)
	df <- df[0,]
	for(n in seq(length(x))) {
		for(mat in plot.mats) {
			df.row <- data.frame(
				t=x[[n]]$t,
				i=market.rate(mkt=x[[n]],mat=mat),
				mat=mat
				)
			df <- rbind(df,df.row)
		}
	}
	# Preprocess data for plotting (interest rates as %, stagger x coordinates)
	df$i <- df$i*100
	stagger.scale <- 35 # scaling constant
	stagger.vec <- (seq(length(plot.mats))-median(seq(length(plot.mats))))/stagger.scale
	df <- merge(df,data.frame(mat=plot.mats,stagger=stagger.vec))
	df <- sort.data.frame(df,formula=~t+mat)
	df$t <- df$t + df$stagger
	# Plot
	cols <- cols
	plot(i~t,data=df,col=cols,xlim=c(start.t,end.t),xlab=xlab,ylab=ylab,...)
	for(n in seq(length(plot.mats))) {
		df.line <- subset(df,mat==plot.mats[n])
		# Create a dataframe to plot a constant rate until the rate changes
		df.line <- df.line[rep(seq(nrow(df.line)),each=2),] # replicate each line
		new.ts <- df.line[seq(1,nrow(df.line),2),"t"] # Grab the set of times
		new.ts <- c(new.ts[seq(2,length(new.ts))],end.t) # Shift them by one
		df.line[seq(2,nrow(df.line),2),"t"] <- new.ts # Replace the duplicated rows with our shifted times
		# Now plot
		lines(i~t,data=df.line,col=cols[n],...)
	}
	invisible(df)
}

# - Account functions (Accounts hold portfolios and reinvest the proceeds into cash portfolios) - #
# Make account object from its components
account <- function(prts,hist.mkt,t.issue=hist.mkt[[1]]$t) {
	# Validate inputs
	prt.classes <- as.data.frame(sapply(prts,class))[1,]
	if(any(prt.classes!="portfolio.bond" & prt.classes!="cash")) { stop("Currently supported objects: portfolio.bond, cash.\n")}
	if(any(class(hist.mkt)!=c("history.market","market"))) { stop("Market must be a history.market object. \n Stock objects not yet supported.\n") }
	# Calculate original/issue market and time
	mkt.orig <- current.market(hist.mkt,t.issue)
	# Add in a cash account if one isn't specified yet
	if(all(prt.classes!="cash")) {
		prts[[length(prts)+1]] <- cash(value=0,mkt=mkt.orig,name="cash")
	}
	# Name elements of portfolio list
	names(prts) <- unlist(lapply(prts,function(x) x$name))
	# Return
	res <- list(prts=prts,hist.mkt=hist.mkt,mkt.orig=mkt.orig,t.issue=t.issue,pv.issue=sum(sapply(prts,pv)))
	structure(res,class=c("account"))
}

# Print basic account characteristics and portfolio components
print.account <- function(x,...) {
	cat("$",x$pv.issue," invested into account in year ",x$t.issue," as follows:\n",sep="")
	portfolio.types <- table(unlist(lapply(sapply(x$prts,class),last)))
	for(p in seq(length(portfolio.types))) {
		if(p!=1) { cat(" and ",sep="") }
		cat(portfolio.types[p]," ",names(portfolio.types)[p],sep="")
	}
	cat(" portfolio(s).\n",sep="")
	for(n in seq(length(x$prts))) {
		cat(names(x$prts)[n],": ",sep="")
		print(x$prts[[n]])
	}
	invisible(x)
}

# - History.account functions - #
# Create history.account object from accounts
history.account <- function(accts,...) {
	# Validate input
	if(class(accts)!="list" | any(sapply(accts,"class")!="account")) { stop("accts must be a list of account objects.\n") }
	# Return object
	structure(accts,class=c("history.account","account"))
}
# Print history.account
print.history.account <- function(x,...) {
	cat("Account history\n----------------------------------------------------------------------\n")
	print(as.data.frame(x))
	invisible(x)
}
# Plot history.account
plot.history.account <- function(x,ylab="PV ($)",xlab="Time",...) {
	ts.issue <- sapply(x,function(x) x$t.issue)
	pvs.issue <-  sapply(x,function(x) x$pv.issue)
	
	plot(pvs.issue~ts.issue,xlab=xlab,ylab=ylab,...)
}

# Extract/reinterpret other objects as history.accounts
as.history.account <- function(x,...) {
	UseMethod("as.history.account")
}
as.history.account.sum.account <- function(x,...) {
	return(x$history.account)
}

# - Funds - #
# Create fund object
 # Holdings must be a list of other objects, or NA
 # Can store price.history either as a straight history or as a deviation from an index+ER
fund <- function(holdings,price.history,div.history,expense.ratio=0,index="",name=NA,...) {
	# Validate
	if(!is.na(holdings) | class(holdings)!="list") { stop("Holdings must be a list of other objects (bonds, stocks, etc.).\n") }
	if(is.na(name)) { name <- as.character(runif(1)) }
	# Return our object
	structure(holdings=holdings,price.history=price.history,div.history=div.history,name=name)
}
