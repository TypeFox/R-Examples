require(PairViz)
source(system.file("demo","demo_fns.R", package = "PairViz"))

require(irr)
data(diagnoses)

# Note re data. When I checked the original paper by Fleiss if seems that 
#different raters diagnosed each patient. 
#The paper says
#"A total of 43 psychiatrists provided diagnoses. In the actual study, between 6 and 1p psychiatrists per patient from the pool of 42 were unsystematically selected to diagnose a subject. For the illistrative purpose of this paper, randomly selected diagnoses were dropped to bring the number of assignments per patient down to a constant of 6."

#"The data as analysed here is artificial as we assume there are 6 psychiatrists each rating all 30 patients. In the data listed in the {\tt irr} package, the six ratings on a patient are sorted by their nominal values. This explains the upward trend of ratings and the apparent agreement between adjacent raters.



lev <- levels(diagnoses[,1])
d <-as.data.frame(lapply(diagnoses,match,lev))


lev <- as.numeric(levels(factor(d[,1])))
names(lev) <- c("D","PD","S","N","O")
nlev <- length(lev)
nraters <- ncol(d)

colnames(d) <- as.character(1:nraters)

col <- rainbow(nlev,alpha=0.6)
cols <- col[as.numeric(d[,1])]


cols <- colour_var(d[,1])

bar.col <- "#FDCDAC"
pcp.main <- "Rating psychiatric patients"


#cols <- "grey50"
#bar.col <-"grey80"
#-------rater agreement-----


a <-outer(d, d, Vectorize(function(x,y)sum(x==y)))
a <- a/nrow(d) # the proportion of patients for which two raters agree

#a <-outer(d, d, Vectorize(function(x,y) kappa2(cbind(x,y))$value))


ad <- dist2edge(a)

#-------
#standard PCP, with bars showing rater agreement
o <- eulerian(1-a)

dev.new(width=3.5,height=2)
par(tcl = -.2, cex.axis=.5,mgp=c(3,.2,0),cex=.8,cex.main=.8)


guided_pcp(d,edgew=ad,path=o,pcp.scale=FALSE,pcp.col=cols,lwd=1,bar.axes=TRUE,bar.col=bar.col,
pcp.mar=c(1,1,1,1),main= "pcp.main")



#---spread out the data---------


ds <- rater_spreadout(d,lev)


o <- 1:nraters

#---guided pcp with standard path---------


dev.new(width=3.5,height=2)
par(tcl = -.2, cex.axis=.5,mgp=c(3,.2,0),cex=.8,cex.main=.8)

guided_pcp(ds$data ,edgew=ad,path=o,pcpbars=ds$bars,pcp.col=cols,lwd=2,pcp.scale=FALSE,bar.col="grey80",pcpfn=catpcp, pcpbars.axis.at=1:nlev, pcpbars.axis.labels=names(lev),main= pcp.main,pcp.mar=c(1,1,1,1),bar.axes=TRUE)




#---guided pcp with eulerian path---------

o <- eulerian(1-a)


dev.new(width=6.5,height=2)
par(tcl = -.2, cex.axis=.5,mgp=c(3,.2,0),cex=.8,cex.main=.8)


guided_pcp(ds$data ,edgew=ad,path=o,pcpbars=ds$bars,pcp.col=cols,lwd=2,pcp.scale=FALSE,bar.col="grey80",pcpfn=catpcp, pcpbars.axis.at=1:nlev, pcpbars.axis.labels=names(lev),main= pcp.main,pcp.mar=c(1,1,1,1),bar.axes=TRUE)




#-----------------


#Find the best ordering of diagnoses


diag_confusion <- function(d,lev){
	nlev <- length(lev)	
	nraters <- ncol(d)
	rpairs <- combinations(nraters,2)
    sim <- matrix(0,nlev,nlev)
    colnames(sim) <- rownames(sim)<- names(lev)
    for (i in 1:nrow(rpairs)){
	  r <- d[,rpairs[i,]]
	  new <- table(factor(r[,1],levels=lev),factor(r[,2],levels=lev))
      sim <- sim+ new
	 }
   sim1 <- t(sim)
   diag(sim1)<- 0
   sim <- sim+sim1
   return(sim)
	}

sim <- diag_confusion(d,lev)
# sim[i,j]   gives the number of patients who were given diagnosis i by one rater and j by another


best <- order_best(max(sim) - sim) # the best permutation

sim[best,best] # with this permutation, frequently confused diagnoses are adjacent
  
 ds <- rater_spreadout(d,best)

o <- eulerian(1-a)


dev.new(width=6.5,height=2)
par(tcl = -.2, cex.axis=.5,mgp=c(3,.2,0),cex=.8,cex.main=.8)

guided_pcp(ds$data ,edgew=ad,path=o,pcpbars=ds$bars,pcp.col=cols,lwd=1,pcp.scale=FALSE,bar.col="grey80",pcpfn=catpcp, pcpbars.axis.at=1:nlev, pcpbars.axis.labels=names(lev)[best],main= pcp.main,pcp.mar=c(1,1,1,1),bar.axes=TRUE)



#----------draw the similarity matrices- for the diagnoses
dcols <- matrix("wheat",nlev,nlev)
diag(dcols) <- "paleturquoise"
dcols[abs(row(dcols) - col(dcols)) == 1] <- "springgreen1"    

dcols <- dcols[,nlev:1]

dev.new(width=4,height=4)
par(mar=c(2.5,2.5,1,1))
s <- sim
s <- s[,nlev:1]  
  
table_plot(sqrt(s),sqrt(s),col=dcols,spacex=.1,spacey=.1)
# Area of off diagnonal cells is proportional to the total number of times that the row and column diagnoses are given to a patient.Area of diagonal cells is proportional to the total number of times that two differnt raters give the row (column) diagnosis to a patient. 



dev.new(width=4,height=4)
par(mar=c(2.5,2.5,1,1))
s <- sim[best,best] # with this permutation, frequently confused diagnoses are adjacent
s <- s[,nlev:1]  
  
table_plot(sqrt(s),sqrt(s),col=dcols,spacex=.1,spacey=.1)
# In this version, frequently confused diagnoses are adjacent so bigger cells are placed on the off diagonal.


#----------draw the similarity matrices- for the raters

rcols <- matrix("wheat",nraters,nraters)
diag(rcols) <- "paleturquoise"

rcols <- rcols[,nraters:1]
dev.new(width=4,height=4)
par(mar=c(2.5,2.5,1,1))
s <- a[nraters:1,]  
  
table_plot(sqrt(s),sqrt(s),col=rcols,spacex=.1,spacey=.1)
# Area of diagonal cells is proportional to the total number of patients rated by each rater. As each rater rated all 30 patients, diagonal cell areas are identical.  Area of off diagnonal cells is proportional to the total number of times that the row and column raters gave the a patient the same rating.




bestr <- order_best(1 - a) # the best permutation
s <- a[bestr,bestr] # with this permutation,  raters with high agreement are adjacent

dev.new(width=4,height=4)
par(mar=c(2.5,2.5,1,1))
s <- s[nraters:1,]  
  
table_plot(sqrt(s),sqrt(s),col=rcols,spacex=.1,spacey=.1)


#-------displays of patient profiles

#  patient plot


	
	
dheights <- apply(as.matrix(d),1, function(di) table(factor(di,levels=lev)))
colnames(dheights) <- 1:ncol(dheights)
rownames(dheights) <- 1:nrow(dheights)


library(RColorBrewer)
  
pcol <- brewer.pal(nlev, "Set3")


dev.new(width=6,height=3)
par(mar=c(2,2,2,2))
par(tcl = -.2, cex.axis=.45,mgp=c(3,.3,0))

barplot(dheights,col=pcol)  


#- next reorder patients so most similar are adjacent
pdist <- dist(t(dheights),method="manhattan")


require("seriation")
po	<- get_order(seriate(pdist,method="OLO",control=(list(method="average"))),1)



dev.new(width=6,height=3)
par(mar=c(2,2,2,2))
par(tcl = -.2, cex.axis=.45,mgp=c(3,.3,0))

barplot(dheights[,po],col=pcol)  

#- and the legend..

dev.new(width=1.2,height=3)
par(mar=c(0,3,0,.5))

barplot(rep(1,nlev),col=pcol,horiz=TRUE,space=0, axes=FALSE,names.arg=names(lev)[best],las=2,cex.names=.8)

