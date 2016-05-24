# Demonstration of maRketSim capabilities on the bond side


# - Simple bonds - #
mkt1 <- market(market.bond(i=.05),t=0) # All that is required to specify a bond market is an interest rate
mkt1C <- market(market.bond(i=.1),t=0)

bnd.A <- bond(mkt=mkt1,mat=5) # Bonds can be specified by maturity
bnd.B <- bond(mkt=mkt1,dur=2.5) # or duration, in which case the maturity under prevailing interest rates is calculated
bnd.C <- bond(mkt=mkt1C,mat=15)

bnd.A # You can display the basic characteristics of a bond
summary(bnd.A,mkt1) # Or more sophisticated information like duration and convexity

# - Bonds in time and yield curves - #
mkt2 <- market(market.bond(yield.curve=quote(0.01 + log10( mat + 1 )/ 20),MMrate=.01),t=2) #yield curve must be in format specified here.  t=2 implies this is a rate change in the future
mkt1B <- market(market.bond(yield.curve=quote(0.01 + log10( mat + 1 )/ 20),MMrate=.01),t=0) # we'll need this guy later to demonstrate automatic portfolio generation

sum.bnd.A <- summary(bnd.A,mkt2) # Now we're evaluating the same bond two years later, with the intervening coupon payments disappearing into the ether (accounts will address that)
str(sum.bnd.A) # The summary object has structure with useful quantities to be extracted

# Example of extracting duration
durs <- c()
ts <- seq(0,15,.5)
for(t in ts) {
	d <- summary(bnd.C,market(market.bond(i=.1),t=t))$dur
	durs <- c(durs,d)
}
plot(durs~ts,main="Duration vs. the passage of time",xlab="Time",ylab="Duration")

# - Portfolios of bonds - #
prt1 <- portfolio(name="bond1",bonds=list(bnd.A,bnd.B,bnd.C),mkt=mkt1) 

prt1 # Display the bonds in the portfolio
summary(prt1) # Display the portfolio's characteristics under its original market conditions
summary(prt1,mkt=mkt2) # Display the portfolio's characteristics under new market conditions
as.data.frame(prt1) #Another way of looking at the portfolio.  Useful for exporting to user-written functions or spreadsheets.

# Create random portfolios of bonds with certain portfolio characteristics
cat("Portfolio generation to a specified duration is currently very fragile due to a lack of a closed-form solution.\n You may have to try it a few times to make it work.\n")
prt2 <- genPortfolio.bond(10,mkt=mkt1B,dur=5,dur.sd=2,name="bond2")
prt2
summary(prt2)
cat("Duration of our generated portfolio is",round(abs(5-summary(prt2)$portfolio.sum$dur),2),"away from 5.\n") #Actually, the algorithm for generating random portfolios is pretty crude.  Please e-mail the maintainer with suggestions or code for ways you want it done.

# - Market histories - #
mkt3 <- market(market.bond(yield.curve=quote(mat/75+0.02),MMrate=.02),t=3)
h.mkt.simple <- history.market(list(mkt1,mkt2,mkt3))
plot(h.mkt.simple) # Shows how yield curves are handled
plot(h.mkt.simple,plot.mats=c(1,3,5),start.t=0,end.t=5,plot.MMrate=FALSE) # Shows how to change time period plotted,  how to change the maturities plotted, and how to turn off plotting the money market rate

h.mkt.updown <- genHistory.market(
	i.fxn=quote(1/(10*exp(t))*t^2+.02),
	start.t=0,end.t=5,f=.5
)
plot(h.mkt.updown)

# - Accounts - #
# Creating accounts
prts <- list(prt1,prt2,cash(name="cash1",value=300,mkt=mkt1))
acct1 <- account(prts=prts,hist.mkt=h.mkt.updown)
acct2 <- account(prts=prts,hist.mkt=h.mkt.simple)

# Looking at account
acct1
cat("Observe that the value invested doesn't equal the sum of the par values of the bonds plus the cash holdings!\n")
cat("What happened?  Recall that bnd.C we created somewhat nonsensibly with a different prevailing interest rate.\n")
cat("Therefore, although its par is $1000, its coupon is higher, and it immediately became worth more.\n")
cat("pv() will calculate the present value of a particular object, such as our bond: $",pv(bnd.C,mkt=mkt1),"\n")
cat("It works for other objects too.  Here's the present value of prt1: $",pv(prt1,mkt=mkt1),"\n")

# Seeing what happens in the future
summary(acct1,t=5,rebal.function.args=list(min.bond.size=1000,new.bond.dur=5,new.bond.mat=NA)) # this is using the default rebalance function



### More examples ###
# Flat yield curve
mkt1 <- market(market.bond(i=.05),t=0)
mkt2 <- market(market.bond(i=.07),t=0)
prt.ladder <- portfolio(name="Ladder",bonds=list(
		bond(mkt=mkt1,mat=1),
		bond(mkt=mkt1,mat=2),
		bond(mkt=mkt1,mat=3),
		bond(mkt=mkt1,mat=4),
		bond(mkt=mkt1,mat=5),
		bond(mkt=mkt1,mat=6),
		bond(mkt=mkt1,mat=7),
		bond(mkt=mkt1,mat=8),
		bond(mkt=mkt1,mat=9),
		bond(mkt=mkt1,mat=10)
	),mkt=mkt1)
prt.bul <- portfolio(name="Bullet",bonds=list(
		bond(mkt=mkt1,mat=4),
		bond(mkt=mkt1,mat=4),
		bond(mkt=mkt1,mat=5),
		bond(mkt=mkt1,mat=5),
		bond(mkt=mkt1,mat=6),
		bond(mkt=mkt1,mat=6),
		bond(mkt=mkt1,mat=7),
		bond(mkt=mkt1,mat=7),
		bond(mkt=mkt1,mat=4),
		bond(mkt=mkt1,mat=5)
	),mkt=mkt1)
prt.ladder
prt.bul
summary(prt.ladder)
summary(prt.bul)
cat("After changing interest rates from 5% to 7%:\n")
summary(prt.ladder,mkt=mkt2)
summary(prt.bul,mkt=mkt2)

# Sharply upward-sloping yield curve (MMrate=.026, 1-year 5%, 5-year 9.5%, 10-year 12%, 20-year 15%)
mkt.bond.up <- market.bond(yield.curve=quote(log(mat+1.7)/20))
##mkt.bond.up <- market.bond(i=.05)
mkt.up0 <- market(mkt.bond.up,t=0)
prt.ladder <- portfolio(name="Ladder",bonds=list(
		bond(mkt=mkt.up0,mat=1),
		bond(mkt=mkt.up0,mat=2),
		bond(mkt=mkt.up0,mat=3),
		bond(mkt=mkt.up0,mat=4),
		bond(mkt=mkt.up0,mat=5),
		bond(mkt=mkt.up0,mat=6),
		bond(mkt=mkt.up0,mat=7),
		bond(mkt=mkt.up0,mat=8),
		bond(mkt=mkt.up0,mat=9),
		bond(mkt=mkt.up0,mat=10)
	),mkt=mkt.up0)
prt.bul <- portfolio(name="Bullet",bonds=list(
		bond(mkt=mkt.up0,mat=4),
		bond(mkt=mkt.up0,mat=4),
		bond(mkt=mkt.up0,mat=5),
		bond(mkt=mkt.up0,mat=5),
		bond(mkt=mkt.up0,mat=6),
		bond(mkt=mkt.up0,mat=6),
		bond(mkt=mkt.up0,mat=5),
		bond(mkt=mkt.up0,mat=5),
		bond(mkt=mkt.up0,mat=4),
		bond(mkt=mkt.up0,mat=5)
	),mkt=mkt.up0)
mkt.up.hist <- history.market(list(
		market(mkt.bond.up,t=0),
		market(mkt.bond.up,t=40)
	))
#plot(mkt.up.hist)
acct.ladder <- account(prts=list(prt.ladder),hist.mkt=mkt.up.hist)
acct.bul <- account(prts=list(prt.bul),hist.mkt=mkt.up.hist)

sum.ladder <- summary(acct.ladder,t=20,rebal.function.args=list(min.bond.size=1000,new.bond.mat=10,new.bond.dur=NA,sell.mat=NA))
sum.bul <- summary(acct.bul,t=20,rebal.function.args=list(min.bond.size=1000,new.bond.mat=NA,new.bond.dur=3.8,sell.mat=0))
plot(sum.ladder,main="Ladder")
plot(sum.bul,main="Bullet")
