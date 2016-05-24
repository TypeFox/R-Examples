library(PairViz)


#----------
regression_graph <- function(y,x)
  {
  	g <- mk_hypercube_graph(colnames(x)) 
  	nodeDataDefaults(g, "residuals") <- 0
  	nodeDataDefaults(g, "sse") <- 0
    edgeDataDefaults(g,"weight") <- 1

  
    preds <- nodes(g)
    res <-sapply(preds,function(z) {
      p <- strsplit(z,character(0))[[1]]
	  if (p[1]=="0") 
	     r <- y - mean(y)
	  else {d <- as.data.frame(x[,p])
	  r <- lm(y ~., data=d)$residuals}
	  
	  return(r)})
	  
	  
    colnames(res) <- preds
    for (i in 1:ncol(res)){
      r <- res[,i]
      nodeData(g,nodes(g)[i],"residuals"	) <- list(r)
      nodeData(g,nodes(g)[i],"sse"	) <- sum(r*r)
      }
    for (n in nodes(g)) {
    	m <- edges(g,n)[[1]]
    	ew <- unlist(nodeData(g,m,"sse")) - nodeData(g,n,"sse")[[1]]
    	edgeData(g,n,m,"weight") <- abs(ew)
    	}  	

    return(g)
  	}
  	
#----------



library(alr3)
data(sleep1)

data <- sleep1



# logging the brain and body weights
data[,c(4,5)] <- log(data[,c(4,5)])
colnames(data) <- c("SW","PS" ,"TS" ,"Bd", "Br","L","GP","P" ,"SE" , "D"  )

#data <- na.omit(data)

data <- data[,c(1,2,4,5,6)]

data <- na.omit(data)
y <- data[,4]
x <- data[,-4]
colnames(x) <-c("A","B","C","D")


g <- regression_graph(y,x)


res <- sapply(nodeData(g,attr="residuals"),identity) # matrix of residuals

# g is  eulerian, since it is Q4
eulerian(g,weighted=FALSE)

# This version uses weights, and picks as a start a node connected to the loweset weight edge, which happens to be node AD
# To do a backwards selection version, specify start as "ABCD"


o <- eulerian(g,start="ABCD")

	
ew <- NULL	
for (i in 2:length(o)) ew <- c(ew,
  nodeData(g,o[i],"sse")[[1]] - nodeData(g,o[i-1],"sse")[[1]])


dev.new(width=8,height=3)


par(tcl = -.2, cex.axis=.4,mgp=c(3,.3,0),cex.main=.8)
pcp0 <- function(...){pcp(...)
	abline(h=0,col="grey70",lwd=2)}

cols <- desaturate_color(rainbow(10,alpha=0.7))
cols <- cols[cut(rank(res[,"ABCD"]),10,labels=FALSE)]

guided_pcp(res,path=match(o,nodes(g)),pathw=ew,bar.col="#FDCDAC",pcp.scale=FALSE,pcpfn=pcp0,lwd=2,pcp.col=cols,main="Sleep data: Model residuals.",pcp.mar=c(1,1,1.5,1))

guided_pcp(res,path=match(o,nodes(g)),pcp.scale=FALSE,pathw=ew,bar.col="grey80",pcpfn=pcp0,lwd=1,pcp.col="grey40",main="Sleep data: Model residuals.",pcp.mar=c(1,1,1.5,1))



cols <- rep("grey70",nrow(data))
cols[rownames(data) =="Human"] <- "red"
cols[rownames(data) =="Asian_elephant"  ] <- "black"
cols[rownames(data) =="Big_brown_bat"    ] <- "purple"
cols[rownames(data) =="Little_brown_bat"      ] <- "purple"
cols[rownames(data) =="Echidna"       ] <- "blue"
cols[rownames(data) =="Lesser_short-tailed_shrew"] <- "cyan"
cols[rownames(data) =="Ground_squirrel"  ] <- "magenta"

  

#---------------
#remove intercept-only node
g1 <- g
g1 <- removeNode("0" ,g1)
#Graph is not even- so construction of eulerian requires extra edges
o <- eulerian(g1,start="ABCD")
#In this case extra edges are beween A-D and C-B

ew <- NULL	
for (i in 2:length(o)) ew <- c(ew,
  nodeData(g1,o[i],"sse")[[1]] - nodeData(g1,o[i-1],"sse")[[1]])

	


dev.new(width=8,height=3)

par(tcl = -.2, cex.axis=.4,mgp=c(3,.3,0))

guided_pcp(res,path=match(o,nodes(g1)),pcp.scale=FALSE,pathw=ew,bar.col="#FDCDAC",pcpfn=pcp0,lwd=2,pcp.col=cols,main="Sleep data: Model residuals.",pcp.mar=c(1,1,2,1))



#---------------

#remove nodes that do not include c

g1 <- g
for (n in nodes(g1)) {
	if (! ("C" %in% strsplit(n,"")[[1]]))
	g1 <- removeNode(n,g1)
	}

o <- eulerian(g1,start="ABCD")

ew <- NULL	
for (i in 2:length(o)) ew <- c(ew,
  nodeData(g1,o[i],"sse")[[1]] - nodeData(g1,o[i-1],"sse")[[1]])

ecols <- rep("#FDCDAC",length(ew))


for (i in 2:length(o)) {
	a <- o[i]
	b <- o[i-1]
	if (is.na(match(b,edges(g1,a)[[1]]))) ecols[i-1] <- "grey70"
	}
ecols[10] <- "grey70"
dev.new(width=5,height=3)

par(tcl = -.2, cex.axis=.4,mgp=c(3,.3,0))

guided_pcp(res,path=match(o,nodes(g)),pcp.scale=FALSE,pathw=ew,bar.col=ecols,pcpfn=pcp0,lwd=2,pcp.col=cols,main="Sleep data: Model residuals.",pcp.mar=c(1,1,2,1))




