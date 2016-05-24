####### Code for constructing figures of Section 3.3 (Hurley and Oldford, 2008) ###################


require(PairViz)
mt9 <- as.matrix(mtcars)[c(7,28,27,31,30,1,12:14),1:7]
rownames(mt9) <- 1:9


col4 <- c("red","blue","green","purple")

colh <- "orange"
cole <- "magenta"




#col4 <- rep("grey50",4)
#colh <- "grey50"
#cole <- "grey50"

o <- hpaths(7)
# for consistency with the paper execute the following line
o <- o -1 ; o[o==0] <- 7 # changes the `join' from 1 to 7

o <- rbind(1:7,o)




mains <- c("Dataset order H0", "Order H1","Order H2" ,"Order H3")
for (i in 1:nrow(o)) {
	dev.new(width=2,height=2)
	
	stars(mt9[,o[i,]], len = 0.8,col.stars=rep(col4[i], nrow(mt9)), 
	 key.loc = NULL, mar=c(.6,0,0,0),   cex=.6,radius = FALSE)
	}
	
#------------------------------------
	
	
dev.new(width=2.5,height=2.5)
stars(mt9[,as.vector(t(o[-1,]))], len = 0.8,col.stars=rep(colh, nrow(mt9)), 
	 key.loc = NULL,  mar=c(.6,0,0,0),cex=.7,radius = FALSE)



#------------------------------------

m <- cor(mt9)
o <- eulerian(-m)
o <- o[-length(o)] # tour is closed, last element is the same as first and is uneeded for star glyps


dev.new(width=2.5,height=2.5)

stars(mt9[,o], len = 0.8,col.stars=rep(cole, nrow(mt9)), 
	 key.loc = NULL, mar=c(.6,0,0,0),   cex=.7,radius = FALSE)
	 
	 

#------------------------------------
dev.new(width=4,height=2.5) 
par(mar=c(2,5,1,1))

h <- hclust(dist(scale(mt9)), "average")
plot(h,main="",xlab="",cex=.7,cex.lab=1,axes=FALSE)
axis(2,cex.axis=.7)

stars(mt9[h$order,o], col.stars=rep(cole, nrow(mt9)),
locations=cbind((1:9),c(.5,.5,1,1,-.5,-.5,-.5,1,1)-.7),
      main = "",radius = FALSE,mar=c(0,3,0,0),labels=NULL,len=.5,add=TRUE)
      
#------------------------------------

dev.new(width=4,height=3) 

loc <- prcomp(mt9, scale=TRUE,retx=TRUE)$x

stars(mt9[,o], col.stars=rep(cole, nrow(mt9)),
locations=loc[,1:2]*3, mar=c(.6,0,0,0), cex=.7,
    #  main = "Eulerian",
      radius = FALSE)
