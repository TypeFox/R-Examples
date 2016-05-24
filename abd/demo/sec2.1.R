data(TeenDeaths)
# almost Figure 2.1-1, but the order of the causes is alphabetical
barchart(deaths ~ cause, TeenDeaths, scales=list(x=list(rot=45)))
# 
# Figure 2.1-1
barchart(deaths ~ cause, TeenDeaths, scales=list(x=list(rot=45)))
# Table 2.1-3
table(cut(DesertBirds$count, seq(0,650,by=50)))
# Figure 2.1-2
histogram(~count, DesertBirds,n=12)
densityplot(~count, DesertBirds)
# Figure 2.1-4

histogram(~mass, SockeyeFemales, breaks=seq(1,4,by=0.1))
histogram(~mass, SockeyeFemales, breaks=seq(1,4,by=0.3))
histogram(~mass, SockeyeFemales, breaks=seq(1,4,by=0.5))
plots <- list()
for (b in c(0.1, 0.3, 0.5)) {
  p <- histogram(~mass, data=SockeyeFemales, 
  		breaks = seq(1,4,by=b),
		type='count',
   	 	xlab = "Body mass (kg)"
	)
	plots <- c(plots,list(p))
}
for (i in 1:3)  {
	print(plots[[i]], split=c(i,1,3,1), more=(i<3))
}
