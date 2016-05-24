###################


library(PairViz)
source(system.file("demo","demo_fns.R", package = "PairViz"))


data <- mtcars[,c(1,3:7)]





cols <- colour_var(mtcars[,9])
# transmission type, green=automatic


colnames(data) <- c("Mpg" , "Disp", "Hp"  , "Drat" ,"Wt" ,  "Qsec")

o <- hpaths(1:6,matrix=FALSE)


dev.new(width=7.5,height=2.5)


pc <- rep("grey93",17) 
pc[c(6,12)] <- 0


par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0),cex.main=.8)

pcp(data,order=o,horizontal=TRUE,col=cols,lwd=2,
main = "Hamiltonian decomposition",panel.colors=pc)


# add a correlation guide and find "better" hamiltonians...


corw <- as.dist(cor(data))

o <- weighted_hpaths(-corw, matrix=FALSE)

corw <- dist2edge(corw)

edgew <- cbind(corw*(corw>0), corw*(corw<0))
dev.new(width=7.5,height=3)
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))

guided_pcp(data,edgew, path=o,pcp.col=cols,panel.colors=pc,lwd=2,
         main="Correlation guided Hamiltonian decomposition",bar.col = desaturate_color(c("blue","purple"),.3))


#================================================================


# Sleep data setup

require(alr3)
data(sleep1)
data <- na.omit(sleep1)

# logging the brain and body weights
data[,4:5] <- log(data[,4:5])

colnames(data) <- c("SW","PS" ,"TS" ,"Bd", "Br","L","GP","P" ,"SE" , "D"  )



cols <- desaturate_color(c("red","navy","lightblue3"   ),.6)[cut(rank(data[,6]),3,labels=FALSE)] # colours for pc 


#cols <- colour_var(cut(rank(data[,6]),3))


corw <- as.dist(cor(data))

o <- eulerian(-corw)

#o <- eulerian(-abs(corw))


core <- dist2edge(corw)
edgew <- cbind(core*(core>0), core*(core<0))

dev.new(width=9.5,height=3)
par(tcl = -.2, cex.axis=.45,mgp=c(3,.3,0),cex.main=.8)
guided_pcp(data,edgew, path=o,pcp.col=cols,lwd=1.4,
         bar.col = desaturate_color(c("blue","orange"),.4),main="Eulerian on correlation",bar.ylim=c(-1,1),bar.axes=TRUE)



dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.45,mgp=c(3,.3,0),cex.main=.8)
guided_pcp(data,edgew, path=o,pcp.col=cols,lwd=1.4,
         bar.col = desaturate_color(c("blue","purple"),.5),main="Eulerian on correlation",bar.ylim=c(-1,1),bar.axes=TRUE,zoom=1:26)



#o <- find_path(-corw,   order_best,maxexact=10)
# answer is ...
o <-  c(8, 10,  9,  4 , 5,  7,  6,  2,  3,  1) 
dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.45,mgp=c(3,.3,0))
guided_pcp(data,edgew, path=o,pcp.col=cols,lwd=1.4,
         bar.col = desaturate_color(c("blue","purple"),.5),main="Best Hamiltonian on correlation",bar.ylim=c(-1,1),bar.axes=TRUE
)



#_________utility functions__________

pc_scag_title <- function(title="Parallel Coordinates",sel_scag)
 {if (length(sel_scag)>1) 
   {if (length(sel_scag)==length(scags))
      {title <- paste(title, 
      	               "on all scagnostics.", 
      	               sep=" ")
   	} else {title <- paste(title, 
      	               "on scagnostics:", 
      	               sep=" ")
   		     for (scag in sel_scag) {
	         title <- paste(title,scag, sep=" ")}}
    } else {title <-  	paste(title, 
      	               "on scagnostics:", 
      	               sep=" ")
      	      title <- paste(title,sel_scag, sep=" ")}

	 title}



select_scagnostics <- function(sc,names){
	sc1 <- sc[,names]
	class(sc1) <- class(sc)
	return(sc1)
	}

#____________________________________________
	 
	
library(scagnostics)
library(RColorBrewer)

		 


# The base data for finding weights.
sc <- t(scagnostics(data))

# Scag names and colours
# These are preserved throughout in this order.
scags <- colnames(sc)


scag_cols <- rev(brewer.pal(9, "Pastel1"))
#scag_cols[5] <- scag_cols[6]
names(scag_cols) <- scags


#____________________________________________
sel_scag <- c("Outlying")

sc1 <- select_scagnostics(sc,sel_scag)
# o <- find_path(-sc1,   order_best,maxexact=10)
o <-c( 2 , 4 , 1,  5 , 6,  7 , 3 , 8,  9, 10) #most outlying

dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=o,pcp.col=cols,lwd=1.4,
         main=pc_scag_title("Best Hamiltonian",sel_scag),bar.col = scag_cols[sel_scag],legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,.6))

#____________________________________________

sel_scag <- c("Clumpy")

sc1 <- select_scagnostics(sc,sel_scag)
# o <- find_path(-sc1,   order_best,maxexact=10)
o <- c(5 , 9,  2,  7 , 3 , 1 , 8,  6, 10,  4) # "clumpy"

dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=o,pcp.col=cols,lwd=1.4,
         main=pc_scag_title("Best Hamiltonian",sel_scag),bar.col = scag_cols[sel_scag],legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,.6))
#____________________________________________


sel_scag <- c("Striated","Sparse")

sc1 <- select_scagnostics(sc,sel_scag)
# o <- find_path(-sc1,   order_best,maxexact=10)
o <-  c(4, 10 , 2 , 9 , 1,  7 , 8,  6 , 5 ,3) #most"Sparse"+ Striated"


dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=o,pcp.col=cols,lwd=1.4,
         main=pc_scag_title("Best Hamiltonian",sel_scag),bar.col = scag_cols[sel_scag],legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,.6))
#____________________________________________
sel_scag <- c("Outlying")
sc1 <- select_scagnostics(sc,sel_scag)


# find first 10 axis of eulerian

dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=eulerian, zoom=1:10,pcp.col=cols, lwd=2,
         main=pc_scag_title("Eulerian",sel_scag),bar.col = scag_cols[sel_scag],bar.axes=TRUE,bar.ylim=c(0,.6))
      
#____________________________________________



# find best hamiltonian on selected scagnostics

# o <- find_path(-sc1,   order_best,maxexact=10)
# answer is

o <-  c(3 , 7 , 6,  5 , 9  ,2 , 8,  1,  4 ,10)  #most outlying+clumpy

#o <-c( 2 , 4 , 1,  5 , 6,  7 , 3 , 8,  9, 10) #most outlying
#o <- c( 4, 10,  2,  9 , 1 , 7 , 8 , 6 , 5,  3) #most striated
#o <-c( 5,  8,  4, 10,  6 , 7,  9,  3 , 2 , 1) #most sparse
#o <-  c(4, 10 , 2 , 9 , 1,  7 , 8,  6 , 5 ,3) #most"Sparse"+ Striated"
#o <- c( 3 , 7  6 ,10 , 4,  1  ,8 , 2,  9 , 5)  # #most "Outlying"+ "Clumpy"   +"Sparse" 
#o <- c( 6 , 7  ,5  ,4,  1,  3 , 2,  9 ,10 , 8) # "Monotonic","Convex"
#o <- c(5 , 9,  2,  7 , 3 , 1 , 8,  6, 10,  4) # "clumpy"

#____________________________________________

# find second hamiltonian from decomposition using selected scagnostics
oo <- find_path(-sc1, weighted_hpaths,path1=o)  
dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=oo[2,], pcp.col=cols, lwd=1.4,
         main=pc_scag_title("Second Hamiltonian",sel_scag),bar.col = scag_cols[sel_scag],bar.axes=TRUE,bar.ylim=c(0,.6))
         


         
         
        

# ----scagnostics legend --------------
dev.new(width=1.5,height=3)
par(mar=c(0,4.5,0,.5))

barplot(rep(1,9),col=scag_cols,horiz=TRUE,space=0, axes=FALSE,names.arg=scags,las=2,cex.names=.8)



dev.new(width=1,height=2)
par(mar=c(0,3,0,.5))
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))

sel_scags <- c("Outlying","Clumpy","Striated","Sparse")

barplot(rep(1,length(sel_scags)),col=scag_cols[sel_scags],
horiz=TRUE,space=0, axes=FALSE,names.arg=sel_scags,las=2,cex.names=.8)

# ----Eulerian on thresholded graph --------------

sel_scag <- c("Outlying")

sc1 <- select_scagnostics(sc,sel_scag)

g <- mk_complete_graph(edge2dist(-sc1))
gn <- dn_graph(g,-.2) # removes edges with low scagnostic values, ie edges above -.2
o <- as.numeric(eulerian(gn))


dev.new(width=4.5,height=3)
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))

guided_pcp(data,sc1, path=o, pcp.col=cols, lwd=1.4,
         main=pc_scag_title("NN Eulerian",sel_scag),bar.col = scag_cols[sel_scag],bar.axes=TRUE,bar.ylim=c(0,0.6))


# ----Eulerian on graph with scagnostics outliers--------------


q <- apply(sc,2,quantile, .9)
s.out <- apply(sc,1,function(x) any(x >=q))

sc1 <- t(sc)
sc1[sc1 <= q] <- 0 # values less than q are set to 0

sgrid <- as.matrix(scagnosticsGrid(t(sc)))
sweight <- colSums(sc1)
gs <- ftM2graphNEL(sgrid[s.out,],-sweight[s.out],edgemode="undirected")

sc1 <- t(sc1)
# gs has edges for pairs of variables with a scagnostic index above 90th percentile. The edge weight is a sum of scag indices above q.

o <- as.numeric(eulerian(gs))


dev.new(width=6.5,height=3)
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))


guided_pcp(data,sc1, path=o, pcp.col=cols, lwd=1.4,
         main="Large Scagnostics",bar.col = scag_cols,bar.axes=TRUE,bar.ylim=c(0,2))
