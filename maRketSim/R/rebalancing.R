# maRketSim - Functions to rebalance accounts

# The heirarchy:
# core functions calculate values from values
# market objects calculates from core functions
# bond object and functions calculate from market object
# portfolio calculates from bond object and market object
# account calculates from bond object and market objects held within a portfolio


# Rebalancing schema:
# summary.account loops through each time period, from t.issue to current t.  Each loop, it calls
#  rebalance.account, which takes the account from the last period, updates using the rebal.function passed to summary.account, and
#  returns an account object to be used in the next loop.


# - Rebalancing functions: Intended to be user-replaceable - #
# Returns a list of portfolios
# Since this is intended to be user-replaceable, we make strong assumptions about how things should be done, like only one cash account
# rebal.function.args set either mat or dur, just like in bond()
# sell.mat sells bonds of <= X maturity. Leave as NA to allow bonds to mature.
rebal.function.default <- function(prts,start.mkt,end.mkt,rebal.function.args=list(min.bond.size=1000,new.bond.dur=numeric(),new.bond.mat=5,sell.mat=numeric(),sell.dur=numeric()),...) {
	# TODO: moving from NA to numeric() in arguments made specifying new.bond.dur not work

	min.bond.size <- rebal.function.args$min.bond.size
	new.bond.dur <- rebal.function.args$new.bond.dur
	new.bond.mat <- rebal.function.args$new.bond.mat
	sell.mat <- rebal.function.args$sell.mat
	sell.dur <- rebal.function.args$sell.dur
	if(length(sell.mat)==0) { sell.mat <- 0 } # If no sell.mat is specified, "sell" bonds at 0 when they mature
	if(length(sell.dur)!=0) { sell.mat <- findMat(dur=sell.dur,i=market.rate(mkt=end.mkt,dur=sell.dur,...),...) }
	if(length(new.bond.dur==0)) { new.bond.dur <- NA} # so that it gets passed along properly
	if(is.na(new.bond.dur) & is.na(new.bond.mat)) { stop("Must specify new.bond.dur or new.bond.mat\n") }
	
	end.t <- end.mkt$t
	prt.types <- sapply(lapply(prts,class),function(x) x[length(x)])
	# Condense cash accounts by getting their value at the start of the period
	prts.cash.idx <- seq(length(prt.types))[prt.types=="cash"]
	prts.cash <- prts[prts.cash.idx]
	cash.value <- sum(unlist(sapply(prts.cash,pv)))

	# - Check each bond portfolio to see if any bonds have come due or need to be sold, adding any coupon payments to cash.value, and buying a new bond whenever min.bond.size has been accumulated
	prts.bond <- prts[seq(length(prt.types))[prt.types=="bond"] ]
	for(p in seq(length(prts.bond))) {
		bnds.del <- c()
		bnds <- prts.bond[[p]]$bonds
		for(b in seq(length(bnds))) {
			t.remaining <- (bnds[[b]]$t.issue + bnds[[b]]$mat) - end.t
			# Expired bonds and bonds to be sold
			if(t.remaining>sell.mat) {
				cash.value <- cash.value + bnds[[b]]$cpn
			} else if(t.remaining<=sell.mat) {
				cash.value <- cash.value + bnds[[b]]$cpn + bnds[[b]]$par
				bnds.del <- c(bnds.del,b)
			}
		}
		if(length(bnds.del)>0) {
			#cat("Dropping",length(bnds.del),"expired or expiring bonds from portfolio",prts.bond[[p]]$name,"\n")
			prts.bond[[p]]$bonds <- prts.bond[[p]]$bonds[-bnds.del]
		}
	}
	
	# Purchase new bonds
	if(cash.value>=min.bond.size) {
		# Deduct the new bond from cash
		bnd.par <- floor(cash.value/100)*100
		cash.value <- cash.value - bnd.par
		bnd <- bond(mkt=end.mkt,dur=new.bond.dur,mat=new.bond.mat,par=bnd.par)
		# Add the bond to a bond.rebal portfolio
		if("bond.rebal" %in% names(prts.bond)) {
			prts.bond$bond.rebal$bonds[[length(prts.bond$bond.rebal$bonds)+1]] <- bnd 
		} else {
			prts.bond[[length(prts.bond)+1]] <- portfolio(bonds=list(bnd),mkt=end.mkt,name="bond.rebal")
		}
	}
	
	# Pass other portfolio types through for now
	prts.other <- prts[seq(length(prt.types))[prt.types!="bond" & prt.types!="cash"] ]
	
	# -- Build the portfolio to return -- #
	
	# Delete any empty portfolios if we've created any with our rampant bond deletion
	numbonds <- sapply(prts.bond,function(x) length(x$bonds))
	prts.bond <- prts.bond[numbonds>0]
	# Combine and return
	prts.u <- prts.bond
	prts.u[[length(prts.u)+1]] <- cash(value=cash.value,mkt=end.mkt,name="cash")
	if(length(prts.other)>0) {
		for(o in seq(length(prts.other))) {
			prts.u[[length(prts.u)+1]] <- prts.other[[o]]
		}
	}
	names(prts.u) <- unlist(lapply(prts.u,function(x) x$name))
	
	return(prts.u)
}


# - rebalance() - #
rebalance <- function(obj,start.t,...) {
	UseMethod("rebalance",obj)
}
rebalance.default <- function(obj,start.t,...) {
	stop("Sorry, rebalancing only supported for accounts at the moment.\n")
}
rebalance.account <- function(obj,start.t,rebal.function=rebal.function.default,f,rebal.function.args=list(),...) {
	acct <- obj
	end.t <- start.t+f # Since we're only analyzing one period per iteration of rebalance()
	prts <- acct$prts
	names(prts) <- unlist(lapply(prts,function(x) x$name))
	start.mkt <- current.market(acct$hist.mkt,t=start.t)
	end.mkt <- current.market(acct$hist.mkt,t=end.t)
	
	# Call the rebal.function to update the portfolios
	prts.u <- rebal.function(prts=prts,start.mkt,end.mkt,rebal.function.args,...) 

	acct.u <- account(prts=prts.u,hist.mkt=acct$hist.mkt,t.issue=end.t) #Re-issue the account at the new period 
	invisible(acct.u)
}


# - Account summary - #
# Account summary - Intended to be called by the user
	# advance account into the future, rolling over according to rebalance rules
summary.account <- function(object,t,rebal.function=rebal.function.default,return.history.account=TRUE,f=.5,rebal.function.args=list(),...) {
	acct <- object
	# Store values for output into summary object
	account.start <- acct
	pv.start <- pv(acct)
	t.start <- acct$t.issue
	t.end <- t
	# Loop through periods, rebalancing
	acct.updated <- acct
	acct.hist <- list(acct)
	ts <- seq(acct$t.issue,t-f,f)
	for(tm in ts) {
		acct.updated <- rebalance(acct.updated,rebal.function=rebal.function,start.t=tm,f=f,rebal.function.args=rebal.function.args,...)
		acct.hist[[length(acct.hist)+1]] <- acct.updated 
	}
	# Store final values for output into summary object
	pv.end <- pv(acct.updated,mkt=current.market(acct$hist.mkt,t))
	
	# If specified, return a history.account object containing all the iterated account objects.
	if(return.history.account==TRUE) {
		history.account <- history.account(accts=acct.hist)
	} else {
		history.account <- NA
	}
	
	# Return our summary account object
	res <- list(account.end=acct.updated,pv.end=pv.end,t.end=t.end,account.start=account.start,pv.start=pv.start,t.start=t.start,history.account=history.account)
	structure(res,class=c("sum.account",class(acct)))
}

# Print account summary characteristics, updated to a new time and market
print.sum.account <- function(x,...) {
	cat("Account created with base year ",x$t.start," with an initial value of $",round(x$pv.start,2),"\n",sep="")
	cat("Analysis ended in year ",x$t.end," with a final value of $",round(x$pv.end,2),"\n",sep="")
	cat("Percentage gained:",round(100*(x$pv.end-x$pv.start)/x$pv.start,1),"%\n",sep="")
	cat("Final account contents\n--------------------------------\n")
	print(x$account.end)
	invisible(x)
}


# Plot account summary characteristics over time
plot.sum.account <- function(x,which=c("pv","duration"),n.ticks=7,...) {
	ha <- x$history.account
	ts <- mktTime(ha)
	if("pv" %in% which) {
		pvs <- pv(ha)
	}
	if("duration" %in% which) {
		durs <- duration(ha)
	}
	if(all(c("pv","duration") %in% which)) {
		plot(pvs~ts,...)
		durs.scaled <- (durs-min(durs)) * (max(pvs)-min(pvs))/(max(durs)-min(durs)) + min(pvs)
		dur.scaled.seq <- seq(min(durs.scaled),max(durs.scaled),diff(range(durs.scaled))/n.ticks )
		dur.seq <- seq(min(durs),max(durs),diff(range(durs))/n.ticks)
		axis(4,at=dur.scaled.seq,labels=round(dur.seq,2))
		mtext("duration",side=4)
		points(durs.scaled~ts,col="red",pch=2)
	}
	if("pv" %in% which & !("duration" %in% which) ) {
		plot(pvs~ts,...)
	}
	if("duration" %in% which & !("pv" %in% which) ) {
		plot(durs~ts,...)
	}
}




















# --- Update portfolio and account objects from underlying parameters --- #
	# min.cash - Reinvest all but $min.cash cash into bonds.  If specified, ignores reinvest.bonds
	# Reinvest.bonds - What to do with maturing bonds and coupon payments.  Add to "cash" portfolio, "roll" to new bonds are the options; ignored if min.cash option is specified. 
	# Reinvest.dividends - What to do with stock dividends (future option)
	#,min.cash=NA,reinvest.bonds="cash",reinvest.dividends=TRUE
	#if(is.na(min.cash)&reinvest.bonds!="cash"&reinvest.bonds!="roll") { stop("Bad option:  Reinvest.bonds - What to do with maturing bonds and coupon payments.  Add to 'cash' portfolio, 'roll' to new bonds are the options; ignored if min.cash option is specified. \n") }
	#if(!is.na(min.cash)&!is.numeric(min.cash)) { stop("min.cash must be numeric or NA\n") } 

	#cashflows
	###call cashflow.account
	
#rebal.rules=list(max.cash=1000,reinvest.cash="target",bond.rolling.ladder=TRUE,
#	set.aa.func=function(aa.old,rebal.tgt,max.jump=0.025,...) { # user-supplied AA rebalancing function
#		# First redistribute those that cap (bigjump.select==TRUE), 
#			# looping to the next one that caps (but not redistributing any of its points to ones that already capped)
#		
#		#
#		
#		return(aa.new)
#	}
#	)
#update.account <- function(object,rebal.tgt,rebal.rules,...) {
#	# Config and confirm
#	rt <- rebal.tgt[rebal.tgt>0&rebal.tgt<=1] # Exclude all crazy values, including 0% categories
#	if(all(names(rt)!="cash")) { # Must have a cash category, even if it's 0%
#		rt <- c(rt,0)
#		names(rt)[length(rt)] <- "cash"
#	}
#	rt <- rt[order(names(rt))] # alphabetical order, so it matches up with AA
#	if(sum(rt)!=1) { stop("Rebalancing targets must sum to 100%.\n") }
#	# Initial AA
#	aa.old <- aa(object,force.cash=TRUE)
#	aa.tgt <- rt
#	if(any(names(rt)!=names(aa.old))) { stop("rebalancing target and acct's AA must be the same") }
#	#bigjump.select <- abs(aa.old-aa.tgt)>rebal.rules$max.jump
#	#plusminus <- ifelse(aa.tgt[bigjump.select]>aa.old[bigjump.select],1,-1)
#	#aa.tgt[bigjump.select] <- aa.old[bigjump.select]*(1+plusminus*rebal.rules$max.jump)
#	
#	
#	if(sum(aa.tgt)!=1) { stop("BUG! aa.tgt doesn't sum to 100%") }
#	
#	##### proceed to rebalance according to params
#	
#	return(acct.new)
#}
