# maRketSim - Object manipulation functions

# The heirarchy:
# core functions calculate values from values
# market objects calculates from core functions
# bond object and functions calculate from market object
# portfolio calculates from bond object and market object
# account calculates from bond object and market objects held within a portfolio



# --- Object functions --- #

# Extract market interest rate from a market object
market.rate <- function(mkt,mat=NA,dur=NA,f=.5) {
  require(taRifx)
	# Error check inputs
	if(is.na(mat)&is.na(dur)) stop("Must specify either mat or dur.\n")
	if(class(mkt)!="market.bond"&class(mkt)!="market")  stop("Market must be a market object, created using market().\n")
	# Extract our bond market
	if(class(mkt)=="market.bond") {
		mkt.b <- mkt
	} else {
		mkt.b <- mkt$mkts[[sapply(mkt$mkts,class)=="market.bond"]]
	}
	# Return appropriate rates if they're easy to calculate
	if(!is.na(mat) && mat <= 1/12) {
		market.rate <- mkt.b$MMrate
		if(market.rate>=1) stop("Your market interest rate should be expressed as a proportion (e.g. for a 5% yield, i=.05)\n")
		return(market.rate)

	}
	if(!is.na(mkt.b$i)) { 
		if(!is.na(mat) && mat <=1/12) {
			market.rate <- mkt.b$MMrate
		} else {  ### This still doesn't catch maturities of less than 1/12yr that are specified by e.g. mat=NA,dur=.05, but it's not a big problem because that situation should never arise elsewhere in code
			market.rate <- mkt.b$i
		}
		if(market.rate>=1) stop("Your market interest rate should be expressed as a proportion (e.g. for a 5% yield, i=.05)\n")
		return(market.rate)
	}
	# Calculate intermediate values to use later
	tol <- f*1.01 # tolerance to narrow in on our answer
	# Calculate maturity from duration if no maturity is specified
	if(is.na(mat)&&is.na(mkt.b$i)) {
		### This will be difficult, since the market rate depends on the maturity, and the duration is calculated from the market rate and the maturity.
		### We will need to interate towards the answer.
		# example: dur=6.231105, mat=10,i=.1, mkt=market(yield.curve=quote(mat/100),MMrate=.02)
		# example: dur=4.832336,mat=5,i=.0125,mkt=market(yield.curve=quote(mat^2/10000+.01),MMrate=.02)
		# example for iyc: mat=12,i=.18, dur=4.853306, mkt=market(yield.curve=quote(.3-mat/100),MMrate=.02)
		
		# Determine if the yield curve is inverted
		mat <- 2 # Short-term
		st <- eval(mkt.b$yield.curve)
		mat <- 20 # Long-term
		lt <- eval(mkt.b$yield.curve)
		if(st>lt)  iyc <- -1 ###Do we not need a IYC converter at all?
		if(lt>st)  iyc <- 1
		if(st==lt) stop("If you want a flat yield curve, specify i instead of yield.curve when you make your market.\n")
		# Initialize
		is <- c(mkt.b$MMrate*4)
		err <- c(f)
		mats <- c(findMat(dur=dur,i=last(is),f=f))
		n <- length(mats)
		# Loop
		while (TRUE) {
			mat <- last(mats)
			is <- c(is,eval(mkt.b$yield.curve))
			err <- c(err,(findDur(mat=mat,i=last(is)) - dur) * iyc) # the error of our guess; multiply by the inverted yield curve 1/-1 variable to invert if inverted    # *iyc
			if(sign(last(err))==sign(err[length(err)-1])) { # keep heading in the same direction
				mat <- mat - f*sign(last(err))
			}
			mats <- c(mats,mat)
			N <- length(mats)
			if(N>6) {
				if(diff(range(mats[(length(mats)-6):length(mats)]))<tol)   break
				if(N>150) {
					cat("Loop appears to have run away on us.\n")
          cat("Likely cause is that you have selected a duration above the maximum duration.\n")
          cat("Either select a lower duration or a higher interest rate.\n")
					return(list(is=is,mats=mats,err=err))
				}
			}
		}
		mat <- mats[which.min(abs(err))-1] # the -1 is to shift it, since the structure is a bit messed up ;-)
	}
	# Return the yield curve rate for the given maturity
	if(!is.na(eval(mkt.b$yield.curve)) & is.na(mkt.b$i)) {
		market.rate <- eval(mkt.b$yield.curve)
		if(market.rate>=1) stop("Your market interest rate should be expressed as a proportion (e.g. for a 5% yield, i=.05)")
		return(market.rate)
	} else {
		stop("Something went wrong.\n")
	}
}

# Functions to return cash-flows from bonds and portfolio objects
cashflow <- function(x,mkt,...) {
	UseMethod("cashflow",x)
}
cashflow.default <- function(x,mkt,...) {
	stop("There is no default cash flow function.\n")
}
cashflow.bond <- function(x,mkt,future=TRUE,f=x$f,...) { ### Add in future=FALSE (past cash flows) for use with fv()
	if(future) {
		t.elapsed <- mkt$t - x$t.issue
		if(t.elapsed<0) { stop("Bond was issued after current market time!\n") }
		effective.mat <- x$mat - t.elapsed
		if(effective.mat<0) { stop("Bond expired.  Specify future=FALSE if you want to examine past cashflows.\n") }
		start.t <- ifelse(t.elapsed==0, x$t.issue+.5, ceiling(mkt$t/f)*f) # return either the first coupon if we are at issue or the next coupon date
		end.t <- ifelse(t.elapsed==0, x$mat, start.t+effective.mat)
	}
	else { # past
		start.t <- x$t.issue
		end.t <- ifelse(mkt$t < (x$mat+x$t.issue), mkt$t, (x$mat+x$t.issue))
	}
	ts <- seq(start.t,end.t,x$f)
	cashflows <- rep(x$cpn,length(ts))
	cashflows[length(cashflows)] <- x$cpn + x$par
	
	res <- data.frame(t=ts,cashflow=cashflows)
	return(res)
}
cashflow.portfolio.bond <- function(x,mkt=x$orig.mkt,sort=FALSE,condense=TRUE,future=TRUE,...) {
	if(future) {
		x <- current.portfolio.bond(object=x,mkt=mkt,...)
	}
	else { # past
		x <- current.portfolio.bond(object=x,mkt=mkt,drop.future=TRUE,drop.expired=FALSE,...)
		
	}

		l.cf <- lapply(x$bonds,cashflow,mkt=mkt,future=future,...)
		u.cf <- unlist(l.cf)
		cf.df <- data.frame(
			t=u.cf[substr(names(u.cf),1,1)=="t"],
			cashflow=u.cf[substr(names(u.cf),1,8)=="cashflow"]
			)
		if(sort) {
			cf.df <- sort.data.frame(~t+cashflow,cf.df)
		}
		res <- cf.df
		if(condense) {
			condensed <- tapply(cf.df$cashflow,cf.df$t,sum)
			res <- data.frame(t=as.numeric(names(condensed)), cashflow=condensed)
			rownames(res) <- seq(nrow(res))
			
		}
		return(res)
}
cashflow.cash <- function(x,mkt,future=TRUE,...) {
	res <- data.frame(t=0,cashflow=0)
	res <- res[0,]
	
	stop("make cashflow.cash return actual values, using history.market as well as mkt$MM.frequency)")
	return(res)
}
cashflow.account <- function(x,mkt,...) {
	cf.df <- data.frame(t=NA,cashflow=NA)
	cf.df <- cf.df[0,]
	for(n in seq(length(x$prts))) {
		cf.df <- rbind(cf.df,cashflow(x$prts[[n]],mkt=mkt,future=FALSE))
	}
	cf.df <- subset(cf.df, t >= x$t.issue & t <= mkt$t )
	return(cf.df)
}

# Exclude bonds which have start t after current market time (e.g. which have not yet been purchased), and which have already matured
current.portfolio.bond <- function(object,mkt,drop.expired=TRUE,drop.future=TRUE,...) {
	t.issues <- sapply(object$bonds,function(x) {invisible(x$t.issue)} )
	mats <- sapply(object$bonds,function(x) {invisible(x$mat)} )
	expirations <- t.issues + mats
	drop.bnds <- rep(FALSE,length(t.issues))
	if(drop.expired) { drop.bnds <- (drop.bnds | (expirations<mkt$t)) }
	if(drop.future) { drop.bnds <- (drop.bnds | (t.issues>mkt$t)) }
	select.bnds <- !drop.bnds
	if(sum(select.bnds)!=length(select.bnds)) {
	cat("Keeping ",sum(select.bnds)," of ",length(select.bnds)," bonds. Excluded bonds have",
		ifelse(drop.expired,"expired ","")," ",ifelse(drop.expired&drop.future,"or ",""),ifelse(drop.future,"not yet been purchased",""),".\n",sep="")
	}
	object$bonds <- object$bonds[select.bnds]
	return(object)
}

# Functions to return the PV of different kinds of objects
pv <- function(x,...) {
	UseMethod("pv",x)
}
pv.default <- function(x,...) {
	stop("Only bond, portfolio, cash, and account objects currently supported by pv().\n")
}
pv.bond <- function(x,mkt,...) {
	mat <- x$mat - (mkt$t - x$t.issue)
	findPV(i=x$i,mat=mat,market.rate=market.rate(mkt=mkt,mat=mat,f=x$f),par=x$par,f=x$f)$pv
}
pv.portfolio.bond <- function(x,mkt=x$orig.mkt,...) {
	pv <- sum(sapply(x$bonds,pv,mkt=mkt))
	pv
}
pv.cash <- function(x,...) {
	x$value
}
pv.account <- function(x,mkt=x$orig.mkt,...) {
	pv <- sum(sapply(x$prts,pv))
	pv
}
pv.history.account <- function(x,type="history.account",...) { # returns a vector unlike most pv() functions
	if(type=="history.account") { 
		f <- function(y) pv(y,mkt=y$mkt.orig)
	} else if(type=="bond"|type=="portfolio.bond") {
		f <- function(y) {
			prts.bond <- y$prts[sapply(y$prts,function(z) class(z)[1])=="portfolio.bond"]
			sum(sapply(prts.bond,pv))
		}
	} else if(type=="cash") {
		f <- function(y) {
			prts.cash <- y$prts[sapply(y$prts,function(z) class(z)[1])=="cash"]
			sum(sapply(prts.cash,pv))
		}
	}
	pvs <- sapply(x,f)
	return(pvs)
}

# Functions to return the future value of various objects
fv <- function(x, mkt, ...) {
	UseMethod("fv",x)
}
fv.default <- function(x, mkt, ...) {
	stop("Object type not yet supported.")
}
fv.cash <- function(x,mkt,...) {   ### update to use a market.history object instead
	if(class(mkt)=="market") {
		t.elapsed <- mkt$t - x$t.issue
		fv <- findFV(P=x$value,i=mkt$MMrate,t.elapsed=t.elapsed,compound=x$compound)
	}
	fv
}
fv.bond <- function(x,mkt,compound="continuous",...) {
	pv <- pv(x,mkt)
	cflow <- cashflow(x,mkt)
	t.elapseds <- mkt$t - cflow$t
	cpn <- cflow$cashflow 
	fv.cpn <- mapply(findFV,P=cflow$cashflow,t.elapsed=t.elapseds,MoreArgs=list(i=mkt$MMrate,compound=compound))
	res <- list(price=pv,cash=fv.cpn)
	res
}

# Functions to return the duration of a bond or portfolio object
duration <- function(x,type="modified",...) {
	UseMethod("duration",x)
}
duration.default <- function(x,type="modified",mkt,...) {
	stop("There is no default duration function.\n")
}
duration.bond <- function(x,type="modified",mkt,...) {  ### should this call summary.bond() instead so that it properly handles dates in between coupon payments?
	t.elapsed <- mkt$t - x$t.issue
	effective.mat <- x$mat - t.elapsed
	market.rate <- market.rate(mkt=mkt,mat=effective.mat,f=x$f)
	return(findDur(mat=effective.mat,i=x$i,market.rate=market.rate,f=x$f,type=type))
}
duration.cash <- function(x,type="modified",mkt,...) {
	return(0) # This is a strong assumption (that cash is treated as 0 duration
}
duration.portfolio.bond <- function(x,type="modified",mkt,...) {
	# Error check inputs
	if(type!="modified"&type!="Macaulay") stop("Must specify modified or Macaulay for record type.\n")
	# Duration as weighted average of bond durations
	pvs <- sapply(x$bonds,pv,mkt=mkt)
	pv.total <- pv(x,mkt=mkt)
	durs <- sapply(x$bonds,duration,mkt=mkt,type=type)
	return(sum((durs*pvs)/pv.total))
}
duration.account <- function(x,type="modified",mkt,...) { # This will be mainly called by duration.summary.account, since usually we don't think of accounts as needing market arguments
	# Error check inputs
	if(type!="modified"&type!="Macaulay") stop("Must specify modified or Macaulay for record type.\n")
	# Duration as weighted average of bond durations
	pvs <- sapply(x$prts,pv,mkt=mkt)
	pv.total <- pv(x,mkt=mkt)
	durs <- sapply(x$prts,duration,mkt=mkt,type=type)
	return(sum((durs*pvs)/pv.total))
}
duration.history.account <- function(x,type="modified",...) { # returns a vector, unlike most duration() functions
	durs <- sapply(x,function(y) duration(y,mkt=y$mkt.orig,type=type))
	return(durs)
}
duration.sum.account <- function(x,type="modified",...) {
	if(!("history.account" %in% names(x))) {
		stop("Must set return.history.account=TRUE when running summary.account if you want a duration history\n")
	}
	return(duration(x$history.account))
}

# Functions to restructure objects into data.frames
as.data.frame.bond <- function(x,...) {
	res <- data.frame(i=x$i,mat=x$mat,par=x$par,f=x$f,cpn=x$cpn,t.issue=x$t.issue)
	res
}
as.data.frame.portfolio.bond <- function(x,...) {
	res <- data.frame(i=NA,mat=NA,par=NA,f=NA,cpn=NA,t.issue=NA)
	ll <- lapply(x$bonds,as.data.frame)
	for(z in seq(length(ll))) {
		res <- rbind(res,ll[[z]])
	}
	res <- res[-1,] # drop off our NAs
	res
}
as.data.frame.history.account <- function(x,...) {
	# t, pv, n.bonds, pv.bonds, pv.cash, dur
	t <- mktTime(x)
	pv <- pv(x)
	pv.bonds <- pv(x,type="bond")
	pv.cash <- pv(x,type="cash")
	dur <- duration(x)
	n.prt.bond <- count(x,type="portfolio.bond")
	n.bonds <- count(x,type="bond")
	return(data.frame(t=t,pv=pv,pv.bonds=pv.bonds,pv.cash=pv.cash,dur=dur,n.bonds=n.bonds,n.prt.bond=n.prt.bond,...))
}

# Function to return the proportion of the portfolio in each asset class
aa <- function(x,...) {
	UseMethod("aa")
}
aa.default <- function(x,...) {
	stop("No default method for asset allocation function yet.\n")
}
aa.account <- function(x,sort=TRUE,force.cash=TRUE,...) {
	aa.vec <- unlist(lapply(x$prts,pv)) / pv(x)
	names(aa.vec) <- unlist(lapply(lapply(x$prts,class),last))
	# If force.cash, then we must have a cash category, even if it's 0%
	if(force.cash) {
		if(all(names(aa.vec)!="cash")) {
			aa.vec <- c(aa.vec,0)
			names(rt)[length(rt)] <- "cash"
		}
	}
	# Sort if specified
	if(sort) { aa.vec <- aa.vec[order(names(aa.vec))] }
	
	return(aa.vec)
}

# Function to return the issue time of bonds
issueTime <- function(x,...) {
	UseMethod("issueTime")
}
issueTime.default <- function(x,...) {
	stop("No default method for issueTime function.\n")
}
issueTime.bond <- function(x,...) {
	return(x$t.issue)
}
issueTime.portfolio.bond <- function(x,...) {
	return(unlist(lapply(x$bonds,issueTime)))
}

# Function to return the maturity time of bonds
matTime <- function(x,...) {
	UseMethod("matTime")
}
matTime.default <- function(x,...) {
	stop("No default method for matTime function yet.\n")
}
matTime.bond <- function(x,...) {
	return( x$mat + issueTime(x) )
}
matTime.portfolio.bond <- function(x,...) {
	return(unlist(lapply(x$bonds,matTime)))
}
# Function to return the current time of objects.  Use carefully!
mktTime <- function(x,...) {
	UseMethod("mktTime")
}
mktTime.default <- function(x,...) {
	stop("No default method for market time function\n")
}
mktTime.market <- function(x,...) {
	return(x$t)
}
mktTime.account <- function(x,...) {
	return(x$t.issue)
}
mktTime.history.account <- function(x,...) {
	return(sapply(x,mktTime))
}

# Functions to return the number of objects of different types in a history.account object
count <- function(x,...) { # this is equivalent to length() but length is a primitive and cannot be extended
	UseMethod("count")
}
count.history.account <- function(x,type="history.account",...) {
	if(type=="history.account" | type=="account") { # default behavior for length() 
		res <- length(structure(x,class="list"))
	} else { # Return vector counts (each element is a particular time period) of the requested object type
		if(type=="portfolio") {
			f <- function(y) {length(y$prts)}
		} else if(type=="portfolio.bond") {
			f <- function(y) {
				sum(sapply(y$prts,function(z) class(z)[1])=="portfolio.bond")
			}
		} else if(type=="cash") {
			f <- function(y) {
				sum(sapply(y$prts,function(z) class(z)[1])=="cash")
			}
		} else if(type=="bond") {
			f <- function(y) {
				bond.prts <- sapply(y$prts,function(z) class(z)[1])=="portfolio.bond"
				sum(sapply(y$prts[bond.prts],function(z) length(z$bonds)))
			}
		}
		res <- sapply(x,f)
	}
	return(res)
}


# - History.market handy functions - #

# Function to return the market object at a given time from a history.market object
	# If no time point exists, interpolate by grabbing the previous market and returning a market object with those characteristics but the current time
current.market <- function(hist.mkt,t,...) {
	ts <- unlist(lapply(hist.mkt,function(x) x$t))
	sel.mkt <- seq(length(hist.mkt))[t==ts]
	if(length(sel.mkt)!=1) { # If there's no exact match, interpolate
		sel.mkt <- last(seq(length(hist.mkt))[ts<t])
		mkt <- hist.mkt[[sel.mkt]]
		mkt$t <- t
	} else {
		mkt <- hist.mkt[[sel.mkt]]
	}
	invisible(mkt)
}

# Interpolate between points on a yield curve and return a maturity.  Used for the quoted yield.curve statements
	# constant.max.mat=TRUE if you want it to return interest rates equal to the maximum maturity available if a mat is requested which exceeds the max maturity
connectTheDots <- function(mat,df,constant.max.mat=TRUE,...) {
	if(class(df$mat)!="numeric" | class(df$i)!="numeric") { stop("Data.frame not formatted properly in call to connectTheDots\n") }
	df <- subset(df,!is.na(df$i))
	df <- sort.data.frame(~mat,df)
	if(mat %in% df$mat) {
		i.ret <- df$i[mat==df$mat]
	} else {
		l.mat <- last(df$mat[df$mat < mat])
		if(is.na(l.mat)) { stop("Requested maturity ",mat," is lower than lowest value in data.frame ",first(df$mat),"\n") }
		l.i <- df$i[df$mat==l.mat]
		u.mat <- first(df$mat[df$mat > mat])
		if(is.na(u.mat) & constant.max.mat==FALSE ) { 
			stop("Requested maturity ",mat," is higher than highest value in data.frame ",last(df$mat),"\n") 
		} else if(is.na(u.mat) & constant.max.mat==TRUE) {
			i.ret <- last(df$i)+.000000001*(mat-last(df$mat)) # Add in a tiny tiny fraction to help optimization functions find the right direction
		} else {
			u.i <- df$i[df$mat==u.mat]
			# Interpolate
			m <- (u.i-l.i)/(u.mat-l.mat)
			i.ret <- m*(mat - u.mat) + u.i # point-slope form: y – y1 = m(x – x1)
		}
	}
	if(is.na(i.ret)) { stop("NA's produced by connectTheDots\n") }
	invisible(i.ret)
}

# Function to take a data.frame of columns (i,mat,t) and turn it into a market.history with yield curves
read.yield.curve <- function(df,drop.if.no.MMrate=FALSE,MM.mat=1/12,...) {
	mkts <- list()
	ts <- as.numeric(names(table(df$t)))
	for(tm in ts) {
		cat(tm,"\n")
		df.t <- subset(df,t==tm)
		if( any(df.t$mat<=MM.mat) | drop.if.no.MMrate==FALSE ) { # exclude those without a money market rate (1 month maturity)
			mat.vars <- c('m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12','m13','m14','m15')
			i.vars <- c('i1','i2','i3','i4','i5','i6','i7','i8','i9','i10','i11','i12','i13','i14','i15')
			# Construct the list of interest rates and maturities to substitute in
			yc.list <- list()
			for(n in seq(length(i.vars))) {
				yc.list[[length(yc.list)+1]] <- df.t$i[n]
				names(yc.list)[length(yc.list)] <- i.vars[n]
			}
			for(n in seq(length(mat.vars))) {
				yc.list[[length(yc.list)+1]] <- df.t$mat[n]
				names(yc.list)[length(yc.list)] <- mat.vars[n]
			}
			# Construct our market
			mkts[[length(mkts)+1]] <- market(market.bond(
				yield.curve=substitute(
					connectTheDots(mat,data.frame(mat=c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15),i=c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15)))
				,yc.list)
				),t=tm
			)
			
		}
	}
	hm <- history.market(mkts)
	return(hm)
}
