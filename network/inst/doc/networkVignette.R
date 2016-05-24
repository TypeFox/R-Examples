### R code from vignette source 'networkVignette.Rnw'

###################################################
### code chunk number 1: networkVignette.Rnw:151-153
###################################################
library(network)
set.seed(1702)


###################################################
### code chunk number 2: networkVignette.Rnw:171-173
###################################################
data("flo")
data("emon")


###################################################
### code chunk number 3: networkVignette.Rnw:184-186
###################################################
net <- network.initialize(5)
net


###################################################
### code chunk number 4: networkVignette.Rnw:213-216
###################################################
nmat <- matrix(rbinom(25, 1, 0.5), nr = 5, nc = 5)
net <- network(nmat, loops = TRUE)
net


###################################################
### code chunk number 5: networkVignette.Rnw:218-219
###################################################
summary(net)


###################################################
### code chunk number 6: networkVignette.Rnw:221-222
###################################################
all(nmat == net[,])


###################################################
### code chunk number 7: networkVignette.Rnw:234-236
###################################################
net <- as.network(nmat, loops = TRUE)
all(nmat == net[,])


###################################################
### code chunk number 8: networkVignette.Rnw:242-244
###################################################
nflo <- network(flo, directed = FALSE)
nflo


###################################################
### code chunk number 9: networkVignette.Rnw:248-253
###################################################
nflo[9,]
nflo[9,1]
nflo[9,4]
is.adjacent(nflo, 9, 1)
is.adjacent(nflo, 9, 4)


###################################################
### code chunk number 10: networkVignette.Rnw:260-268
###################################################
network.size(nflo)                 #Number of vertices
network.edgecount(nflo)            #Number of edges
network.density(nflo)              #Network density
has.loops(nflo)                    #Can nflo have loops?
is.bipartite(nflo)                 #Is nflo coded as bipartite?
is.directed(nflo)                  #Is nflo directed?
is.hyper(nflo)                     #Is nflo hypergraphic?
is.multiplex(nflo)                 #Are multiplex edges allowed?


###################################################
### code chunk number 11: networkVignette.Rnw:274-278
###################################################
as.sociomatrix(nflo)
all(nflo[,]==as.sociomatrix(nflo))
all(as.matrix(nflo)==as.sociomatrix(nflo))
as.matrix(nflo,matrix.type="edgelist")


###################################################
### code chunk number 12: networkVignette.Rnw:287-305
###################################################
#Add edges to an empty network
net <- network.initialize(5,loops=TRUE)
net[nmat>0] <- 1                       #One way to add edges
all(nmat==net[,])                      #Should be TRUE
net[,] <- 0                            #Remove the edges
net[,] <- nmat                         #Not quite kosher, but _will_ work....
all(nmat==net[,])                      #Should still be TRUE
net[,] <- 0                            #Remove the edges
for(i in 1:5)                          #Add the hard way!
  for(j in 1:5)
    if(nmat[i,j])
      net[i,j] <- 1
all(nmat==net[,])                      #Should STILL be TRUE
net[,] <- 0                            #Remove the edges
add.edges(net,row(nmat)[nmat>0],col(nmat)[nmat>0])
all(nmat==net[,])                      #When will it all end??
net[,] <- as.numeric(nmat[,])
all(nmat==net[,])                      #When will it all end??


###################################################
### code chunk number 13: networkVignette.Rnw:309-317
###################################################
#Add edges (redux)
net<-network.initialize(5)                  #Create empty graph
add.edge(net,2,3)                           #Create 2->3 edge
net[,]                                      #Trust, but verify
add.edges(net,c(3,5),c(4,4))                #3 and 5 send ties to 4
net[,]                                      #Again, verify edges
net[,2]<-1                                  #Everyone sends ties to 2
net[,]                                      #Note that loops are not created!


###################################################
### code chunk number 14: networkVignette.Rnw:323-328
###################################################
#Deleting vertices
delete.vertices(net,4)                                  #Remove vertex 4
net[,]                                                  #It's gone!
add.vertices(net,2)                                     #Add two new vertices
net[,]                                                  #Both are isolates


###################################################
### code chunk number 15: networkVignette.Rnw:334-338
###################################################
#Retrieving edges
get.edges(net,1)                                  #Out-edges sent by vertex 1
get.edges(net,2,neighborhood="in")                #In-edges to vertex 2
get.edges(net,1,alter=2)                          #Out-edges from 1 to 2


###################################################
### code chunk number 16: networkVignette.Rnw:343-347
###################################################
#Retrieving edge IDs
get.edgeIDs(net,1)                        #Same as above, but gets ID numbers
get.edgeIDs(net,2,neighborhood="in") 
get.edgeIDs(net,1,alter=2)     


###################################################
### code chunk number 17: networkVignette.Rnw:351-354
###################################################
#Vertex neighborhoods
get.neighborhood(net,1)                                    #1's out-neighbors
get.neighborhood(net,2,type="in")                          #2's in-neighbors


###################################################
### code chunk number 18: networkVignette.Rnw:358-364
###################################################
#Deleting edges
net[2,3]<-0                                            #This deletes the 2->3
                                                       #edge
net[2,3]==0                                            #Should be TRUE
delete.edges(net,get.edgeIDs(net,2,neighborhood="in")) #Remove all->2
net[,]                        


###################################################
### code chunk number 19: networkVignette.Rnw:376-379
###################################################
net <- network.initialize(5)
set.network.attribute(net, "boo", 1:10)
net %n% "hoo" <- letters[1:7]


###################################################
### code chunk number 20: networkVignette.Rnw:382-388
###################################################
#List attributes
list.network.attributes(net)

#Retrieve attributes
get.network.attribute(net,"boo")
net %n% "hoo"


###################################################
### code chunk number 21: networkVignette.Rnw:392-395
###################################################
#Delete attributes
delete.network.attribute(net,"boo")
list.network.attributes(net)


###################################################
### code chunk number 22: networkVignette.Rnw:403-417
###################################################
#Add vertex attributes
set.vertex.attribute(net,"boo",1:5)              #Create a numeric attribute
net %v% "hoo" <- letters[1:5]                    #Now, a character attribute

#Listing attributes
list.vertex.attributes(net)                      #List all vertex attributes

#Retrieving attributes
get.vertex.attribute(net,"boo")                  #Retrieve 'em
net %v% "hoo"

#Deleting attributes
delete.vertex.attribute(net,"boo")               #Remove one
list.vertex.attributes(net)                      #Check to see that it's gone


###################################################
### code chunk number 23: networkVignette.Rnw:426-447
###################################################
#Create a network with some edges
net <- network(nmat)

#Add attributes
set.edge.attribute(net,"boo",sum(nmat):1)
set.edge.value(net,"hoo",matrix(1:25,5,5)) #Note: only sets for extant edges!
net %e% "woo" <- matrix(rnorm(25),5,5)     #Ditto
net[,,names.eval="zoo"] <- nmat*6           #Ditto if add.edges!=TRUE

#List attributes
list.edge.attributes(net)

#Retrieving attributes
get.edge.attribute(get.edges(net,1),"boo") #Get the attribute for 1's out-edges
get.edge.value(net,"hoo")
net %e% "woo"
as.sociomatrix(net,"zoo")

#Deleting attributes
delete.edge.attribute(net,"boo")
list.edge.attributes(net)


###################################################
### code chunk number 24: networkVignette.Rnw:462-477
###################################################
#Extract location information
MtSHloc<-emon$MtStHelens%v%"Location"

#Build an incidence matrix based on Local/Non-local/Both placement
MtSHimat<-cbind(MtSHloc%in%c("L","B"),MtSHloc%in%c("NL","B"))

#Convert incidence matrix to a hypergraph
MtSHbyloc<-network(MtSHimat,matrix="incidence",hyper=TRUE,directed=FALSE,
  loops=TRUE)

#Set vertex names, for convenience
MtSHbyloc%v%"vertex.names"<-emon$MtStHelens%v%"vertex.names"

#Examine the result
MtSHbyloc


###################################################
### code chunk number 25: networkVignette.Rnw:489-491
###################################################
plot(nflo, displaylabels = TRUE, boxed.labels = FALSE)
plot(nflo, displaylabels = TRUE, mode = "circle")


###################################################
### code chunk number 26: networkVignette.Rnw:502-507
###################################################
op<-par(no.readonly=TRUE) # cache the plot params
par(mfcol=c(1,2),mar=c(1,1,1,1),cex=0.5) # adjust margins and text size to fit two panels
plot(nflo, displaylabels = TRUE,boxed.labels = TRUE)
plot(nflo, displaylabels = TRUE, mode = "circle")
par(op) # reset the plot params


###################################################
### code chunk number 27: networkVignette.Rnw:513-514
###################################################
plot(emon$MtSi)


###################################################
### code chunk number 28: networkVignette.Rnw:521-522
###################################################
plot(emon$MtSi)


###################################################
### code chunk number 29: networkVignette.Rnw:532-541
###################################################
library(sna)
network.layout.degree <- function(d, layout.par){
      id <- degree(d, cmode = "indegree")
      od <- degree(d, cmode = "outdegree")
      cbind(id, od)
    }
plot(emon$MtStHelens, mode = "degree", displaylabels = TRUE, 
    boxed.labels = FALSE, suppress.axes = FALSE, label.cex = 0.5,
    xlab = "Indegree", ylab = "Outdegree", label.col = 3)


###################################################
### code chunk number 30: networkVignette.Rnw:548-551
###################################################
plot(emon$MtStHelens, mode = "degree", displaylabels = TRUE, 
    boxed.labels = FALSE, suppress.axes = FALSE, label.cex = 0.5,
    xlab = "Indegree", ylab = "Outdegree", label.col = 3)


###################################################
### code chunk number 31: networkVignette.Rnw:559-564
###################################################
plot(MtSHbyloc, displaylabels = TRUE, label = 
    c(network.vertex.names(MtSHbyloc), "Local", "Non-Local"), 
    boxed.labels = FALSE, label.cex = rep(c(0.5, 1), times = c(27, 2)),
    label.col = rep(c(3, 4), times = c(27, 2)), vertex.col = rep(c(2, 5),
    times = c(27, 2)))


###################################################
### code chunk number 32: networkVignette.Rnw:573-578
###################################################
plot(MtSHbyloc, displaylabels = TRUE, label = 
    c(network.vertex.names(MtSHbyloc), "Local", "Non-Local"), 
    boxed.labels = FALSE, label.cex = rep(c(0.5, 1), times = c(27, 2)),
    label.col = rep(c(3, 4), times = c(27, 2)), vertex.col = rep(c(2, 5),
    times = c(27, 2)))


###################################################
### code chunk number 33: networkVignette.Rnw:718-730
###################################################
rnbernexp <- function(n, nv, p = 0.5, onset.hazard = 1, 
    termination.hazard = 1){
      nets <- list()
      for(i in 1:n)
        nets[[i]] <- .Call("rnbernexp_R", network.initialize(nv, 
          directed = FALSE), p, onset.hazard, termination.hazard, 
          PACKAGE = "networkapi.example")
      if(i > 1)
        nets
      else
        nets[[1]]
    }


