# maRketSim - Load various publically available bond data sets and compare a bond fund to a portfolio of bonds

# - Federal Reserve Yield Curve for Treasuries on the secondary market - #
URLdir <- "http://federalreserve.gov/releases/h15/data/Monthly/"
datadir <- "/tmp/"
files <- c('M1','M3','M6','Y1','Y2','Y3','Y5','Y7','Y10','Y20','Y30')
for(f in files) {
	download.file(paste(URLdir,"H15_TCMNOM_",f,".txt",sep=""),paste(datadir,"H15_TCMNOM_",f,".txt",sep=""))
	df <- read.csv(paste(datadir,"H15_TCMNOM_",f,".txt",sep=""),skip=7)
	colnames(df) <- c("date","i")
	df$mat=switch(substr(f,1,1),M=as.integer(substr(f,2,3))/12,Y=as.integer(substr(f,2,3)))
	# Subset out only the half-year periods
	mth <- substr(df$date,1,2)
	df <- df[mth=="12"|mth=="06",]
	
	if(f==files[1]) {
		TreasuryCurve.df <- df
	} else {
		TreasuryCurve.df <- rbind(TreasuryCurve.df,df)
	}
}
splt <- strsplit(as.character(TreasuryCurve.df$date),"/")
TreasuryCurve.df$t <- sapply(splt,function(x) as.integer(x[2])+(as.integer(x[1]))/12)
TreasuryCurve.df <- subset(TreasuryCurve.df,select=c(-date))
TreasuryCurve.df$i <- as.numeric(TreasuryCurve.df$i)/100 # "ND" for no data will be coerced to NAs


# Convert to history.market object
TreasuryCurve.hm <- read.yield.curve(TreasuryCurve.df,drop.if.no.MMrate=TRUE)
plot(TreasuryCurve.hm,main="Treasury Yield Data")

# Compare performance of a rolling portfolio to a bond fund of similar duration
fund.dur <- 5.2 # VFITX ave duration in 2010
fund.weights <- c(1003791750,633466400,476050460,416342700,390484680,355996830,313689000,293081600,264752800,250087200,215642640,210906020,203584120,147038700,136450020,134760780,131523750,122942160,105412320,100536600,91462500,74025360,62334600,54714020,53347530,48630240,42979800,36484500,21205400,14360594,10504410,9177500,7037170,1142886,361589) # market values of the bonds in the fund
fund.weights <- fund.weights/min(fund.weights) # standardaize
fund.mats <- c(65,104,41,66,73,101,45,76,43,51,56,53,75,42,74,89,37,98,131,107,40,67,77,49,7,30,125,128,15,38,119,27,16,66,55)/12 # monthly maturities of bonds in fund, converted to annual
fund.dur.sd <- sd( rep(fund.mats,fund.weights) ) #weighted SD of the maturities of the holdings of VFITX (source https://personal.vanguard.com/us/FundsAllHoldings?FundId=0035&FundIntExt=INT&tableName=Bond&tableIndex=0&sort=marketValue&sortOrder=desc)
prt <- genPortfolio.bond(
	n=50,
	mkt=TreasuryCurve.hm[[1]],
	dur=fund.dur,
	dur.sd=fund.dur.sd,
	name="Rolling Ladder",
	type="random.constrained"
	)
cat("Our generated portfolio has a duration of",round(summary(prt)$portfolio.sum$dur,3),"compared to the fund's duration of",fund.dur,"\n")
acct <- account(prts=list(prt),hist.mkt=TreasuryCurve.hm)
sum.acct <- summary(acct,t=2010.5,rebal.function=rebal.function.default,rebal.function.args=list(min.bond.size=1000,new.bond.dur=fund.dur) )
plot(sum.acct$history.account)

# - Vanguard Intermediate Treasury Fund - #
# -- Price data -- #
pr.df <- VFITX_prices
pr.df$Date <- as.Date(pr.df$Date)
pr.df$month <- format(pr.df$Date,format="%m")
pr.df <- subset(pr.df,month=="06"|month=="12",select=c("Date","Close","month"))
pr.df$month <- as.numeric(pr.df$month)
pr.df$year <- as.numeric(format(pr.df$Date,format="%Y"))
pr.df$day <- as.numeric(format(pr.df$Date,format="%d"))
# Select last available day of each month
by.res <- by(pr.df,list(pr.df$month,pr.df$year),function(x) x[x$day==max(x$day),] )
pr.df <- by.res[[1]]
for(i in seq(2,length(by.res))) {
	if(!is.null(by.res[[i]])) {
		pr.df <- rbind(pr.df,by.res[[i]])
	}
}
pr.df <- subset(pr.df,select=c("Close","month","year"))
names(pr.df)[names(pr.df)=="Close"] <- "p"
# -- Dividend data -- # 
div.df <- VFITX_div
div.df$Date <- as.Date(div.df$Date)
div.df$month <- as.numeric(format(div.df$Date,format="%m"))
div.df$year <- as.numeric(format(div.df$Date,format="%Y"))
# Aggregate into 6-month dividends
div.df$month[div.df$month>=1&div.df$month<=6] <- 6
div.df$month[div.df$month>=7&div.df$month<=12] <- 12
by.res <- by(div.df,list(div.df$month,div.df$year),
	function(x) {
		return(data.frame(div=sum(x$Dividends),N=length(x$Dividends)))
	})
div.df <- by.res[[1]]
for(i in seq(2,length(by.res))) {
	if(!is.null(by.res[[i]])) {
		div.df <- rbind(div.df,by.res[[i]])
	}
}
div.df$month <- as.numeric(rep(dimnames(by.res)[[1]],length(dimnames(by.res)[[2]])))
div.df$year <- as.numeric(rep(dimnames(by.res)[[2]],each=length(dimnames(by.res)[[1]])))
div.df <- subset(div.df,N>=6,select=c(-N)) # Exclude partial half-years
# -- Combine dividend and price data -- #
vfitx.df <- merge(div.df,pr.df)
# -- Calculate comparable performance -- #
vfitx.df <- subset(vfitx.df,year>=2002)
vfitx.df <- sort.data.frame(vfitx.df,~year+month)
N <- nrow(vfitx.df)
shares <- rep(NA,N)
shares[1] <- 50000/vfitx.df$p[1]
for(n in seq(2,N)) {
	# Add dividend payment reinvested into fund
	shares[n] <- (vfitx.df$div[n]*shares[n-1] / vfitx.df$p[n]) + shares[n-1]
}
cat("Total invesment value in 2010 is ",shares[length(shares)]*vfitx.df$p[nrow(vfitx.df)],"\n")
