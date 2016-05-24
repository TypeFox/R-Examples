library(testthat)
assignInNamespace("cedta.override", c(data.table:::cedta.override,"AFM"),"data.table")
assignInNamespace("cedta.override", c(data.table:::cedta.override,"AFMImageAnalyser"),"data.table")

source("./configuration.R")


test_that("generateReport", {
  #   #exportDirectory="C:/Users/mbeauvai/Documents/AccPlatform/person_display_R"
  #   exportDirectory=tempdir()
  #   
  #   AFMImage<-importFromNanoscope("C:\\Users\\i0272418\\Documents\\personnel\\AccPlatform\\person_display_R\\20140521.001_flatten.txt")
  #   AFMImage<-extractAFMImage(AFMImage,15,15,48)
  #   edgeList=c()
  #   
  #   AFMImage@data$h<-AFMImage@data$h*10
  #   AFMImage@data$h[AFMImage@data$h<39.9]<-0
  #   displayIn3D(AFMImage)
  #   #   
  #   #   print(AFMImage@samplesperline)
  #   #   print(AFMImage@lines)
  #   #   print(length(AFMImage@data$h))
  #   #   AFMImage@data$h<-as.numeric(AFMImage@data$h)
  #   #   print(getVertexId(AFMImage,1,1))
  #   #   print(getVertexId(AFMImage,512,512))
  #   # 
  #   #   x=7
  #   #   y=512
  #   #   currentVertexId<-getVertexId(AFMImage,x,y) 
  #   #   print(currentVertexId)
  #   #   print("--")
  #   #   print(getSurroundingVertexesList(AFMImage,x,y))
  #   #   print("--")
  #   #   print(getSurroundingVertexesList(AFMImage,x+1,y))
  #   #   print("--")
  #   #   x=8
  #   #   y=512
  #   #   currentVertexId<-getVertexId(AFMImage,x,y) 
  #   #   print(currentVertexId)
  #   #   print(getSurroundingVertexesList(AFMImage,x,y))
  #   #   print("--")
  #   #   
  #   #   print(AFMImage@data$h[currentVertexId])
  #   #   if (existsEdge(AFMImage, currentVertexId)) edgeList<-c(edgeList, getSurroundingVertexesList(AFMImage,x,y))
  #   
  #   library(igraph)
  #   
  #   edgeList=data.table()
  #   counter<-0
  #   for (x in seq(1: AFMImage@samplesperline)) {
  #     for (y in seq(1: AFMImage@lines)) {
  #       currentVertexId<-getVertexId(AFMImage,x,y)
  #       if (existsEdge(AFMImage, currentVertexId)) {
  #         edgeList<-rbind(edgeList, getSurroundingVertexesList(AFMImage,x,y))
  #         counter<-counter+1
  #       }
  #     }
  #   }
  # 
  #   
  #   print(paste("Vertexes= ", counter))
  #   el=as.matrix(edgeList)
  #   #♀el <- matrix(edgeList, nc = 2, byrow = TRUE) 
  #   
  #   g<-graph_from_edgelist(el[,1:2], directed=FALSE)
  # #   g
  # #   E(g)$weight=as.numeric(el[,3])
  # #   plot(g,layout=layout.fruchterman.reingold,edge.width=E(g)$weight/2)
  # #   
  #   
  #   
  #   
  #   AFMImage_lines=512
  #   AFMImage_sampleperline=512
  #   
  #   
  # define vertices coordinates from names
  
  # display on grid layout
  
  #################################### 
  
  library(igraph)
  
  #AFMImage<-importFromNanoscope("P:\\MAT_DOCUMENTS\\AccPlatform\\person_display_R\\jessem\\bones\\20140521.001_flatten.txt")
  
  AFMImage<-importFromNanoscope("C:\\Users\\i0272418\\Documents\\personnel\\AccPlatform\\person_display_R\\20140521.001_flatten.txt")
  AFMImage<-importFromNanoscope("C:\\Users\\i0272418\\Documents\\personnel\\AccPlatform\\person_display_R\\20140521.009_flatten.txt")
  AFMImage<-extractAFMImage(AFMImage,0,0,128)
  #AFMImage<-extractAFMImage(AFMImage,15,15,256)
#   AFMImage@data$h<-AFMImage@data$h*10
#   AFMImage@data$h[AFMImage@data$h<39.9]<-0

  
  AFMImageNetworksAnalysis = new("AFMImageNetworksAnalysis")
  AFMImageNetworksAnalysis@heightNetworksslider=10
  AFMImageNetworksAnalysis@filterNetworkssliderMin=89.9
  AFMImageNetworksAnalysis@filterNetworkssliderMax=200

  newAFMImage<-copy(AFMImage)
  heights<-newAFMImage@data$h*AFMImageNetworksAnalysis@heightNetworksslider
  heights<-heights+abs(min(heights))
  heights[heights<AFMImageNetworksAnalysis@filterNetworkssliderMin]<-0
  heights[heights>AFMImageNetworksAnalysis@filterNetworkssliderMax]<-0
  heights[heights!=0]<-1
  newAFMImage@data$h<-heights
  displayIn3D(newAFMImage)
  
  #distance=getInverseDistanceMatrix(newAFMImage, 1, 1)
  
  #print(system.time(getInverseDistanceMatrix(newAFMImage)))
  
  
  topologyAFMImage=getTopologyAFMImage(newAFMImage)
  topologyAFMImage@data$h<-topologyAFMImage@data$h*1000
  displayIn3D(topologyAFMImage)

  topologyAFMImage2<-extractAFMImage(topologyAFMImage,0,0,32)
  displayIn3D(topologyAFMImage2)
  
  
  boxplot(topologyAFMImage2@data$h)
  
  localMean<-mean(topologyAFMImage2@data$h)
  heights<-topologyAFMImage2@data$h
  heights[heights<=localMean]<-0
  heights[heights>localMean]<-10
  topologyAFMImage2@data$h<-heights
  displayIn3D(topologyAFMImage2)
  
  
  # 
  # 
  # AFMImageNetworksAnalysis = new("AFMImageNetworksAnalysis")
  # 
  # 
  # 
  # AFMImageNetworksAnalysis@heightNetworksslider=10
  # AFMImageNetworksAnalysis@filterNetworkssliderMin=89.9
  # AFMImageNetworksAnalysis@filterNetworkssliderMax=200
  # 
  # newAFMImage<-copy(AFMImage)
  # heights<-newAFMImage@data$h*AFMImageNetworksAnalysis@heightNetworksslider
  # heights<-heights+abs(min(heights))
  # heights[heights<AFMImageNetworksAnalysis@filterNetworkssliderMin]<-0
  # heights[heights>AFMImageNetworksAnalysis@filterNetworkssliderMax]<-0
  # newAFMImage@data$h<-heights
  # 
  # 
  # 
  # 
  # 
  # getAllDistances<-getInverseDistanceMatrix(newAFMImage, 1, 1)*matrix(newAFMImage@data$h, ncol=newAFMImage@samplesperline)
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # AFMImageNetworksAnalysis<-calculateNetworks(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis, AFMImage= newAFMImage)
  #   
  # 
  # gridIgraphPlot(AFMImage, AFMImageNetworksAnalysis@originalGraph)
  # gridIgraphPlot(AFMImage, AFMImageNetworksAnalysis@skeletonGraph)
  # 
  # 
  # edgeList<-V(g)$name
  # 
  # uniqueVerticesList<-unique(edgeList)
  # uniqueVerticesList
  # # degree de chaque noeud
  # edgeDegreeList<-degree(g, v=uniqueVerticesList, mode = c("all"), loops = FALSE, normalized = FALSE)
  # edgeDegreeList
  # 
  # # liste ordonn'e9e croissante des noeuds en fonction du degree
  # 
  # 
  # 
  # allVertices<-data.table(vertexId=uniqueVerticesList, degree=edgeDegreeList)
  # # get-list of adjacent vertices with degree > 2 (can't remove if degree < 2)
  # allVertices<-allVertices[order(degree)]
  # listOfCandidateVertices<-allVertices[degree==4]
  # 
  # 
  # g3<-delete_vertices(g3, listOfCandidateVertices$vertexId)
  # 
  # gridIgraphPlot(AFMImage, g)
  # gridIgraphPlot(AFMImage, g3)
  # 
  # 
  # ###########################
  # # network analysis
  # 
  # degr<-degree(g)/2
  # hist(degr)
  # length(degr)
  # print(paste("number of ", length(degr[degr<3])))
  # print(paste("number of ", length(degr[degr>=3],"(should be more simplified)")))
  # 
  # decompose(g)
  #   
  # 
  #  hist(degree(g)/2, bin=5)
  #  plot(degree.distribution(g))
  # 
  #  diameter(g)
  #  transitivity(g)
  #  average.path.length(g)
  #  degree(g)
  #  degree.distribution(g)
  # 
  #  
  #  degree<-degree(g)
  #  hist(degree)
  #  betweenness<-betweenness(g)
  #  hist(betweenness)
  #  closeness<-closeness(g)
  #  hist(closeness)
  #  coreness<-coreness(g)
  #  hist(coreness)
  #  
  # ################
  # 
  # 
  # 
  # degreeAFMImage<-copy(AFMImage)
  # betweennessAFMImage<-copy(AFMImage)
  # closenessAFMImage<-copy(AFMImage)
  # corenessAFMImage<-copy(AFMImage)
  # 
  # coordinatesVector<-getNetworkGridLayout(AFMImage, V(g)$name)
  # coordinatesVector$degree<-degree(g)
  # coordinatesVector$betweenness<-betweenness(g)
  # coordinatesVector$closeness<-closeness(g)
  # coordinatesVector$coreness<-coreness(g)
  # 
  # coordinatesVectorDT<-data.table(coordinatesVector)
  # coordinatesVectorDT<-coordinatesVectorDT[order(y,x)]
  # 
  # data=data.table()
  # nextE<-1
  # nextC<-coordinatesVectorDT[nextE]
  # sizeC<-nrow(coordinatesVectorDT)
  # for(y in seq(1,AFMImage@lines)){
  #   for(x in seq(1,AFMImage@samplesperline)){
  #     if (((x==nextC$x)&&(y==nextC$y))&&(nextE<sizeC)) {
  #       data=rbind(data, 
  #                  data.table(x=x,y=y,
  #                             degree=nextC$degree, 
  #                             betweenness=nextC$betweenness, 
  #                             closeness=nextC$closeness,
  #                             coreness=nextC$coreness))
  #       nextE<-nextE+1
  #       nextC<-coordinatesVectorDT[nextE]
  #     }else{
  #       data=rbind(data, data.table(x=x,y=y,degree=0,betweenness=0, closeness=0, coreness=0))
  #     }
  #   }
  # }
  # 
  # 
  # degreeAFMImage@data$x<-data$x
  # degreeAFMImage@data$y<-data$y
  # degreeAFMImage@data$h<-data$degree
  # displayIn3D(degreeAFMImage)
  # 
  # 
  # betweennessAFMImage@data$x<-data$x
  # betweennessAFMImage@data$y<-data$y
  # 
  # hist(data$betweenness)
  # max(data$betweenness)
  # 
  # betweennessAFMImage@data$h<-data$betweenness*16/max(data$betweenness)
  # betweennessAFMImage@data$h[is.na(betweennessAFMImage@data$h)]<-0
  # displayIn3D(betweennessAFMImage)
  # 
  # closenessAFMImage@data$x<-data$x
  # closenessAFMImage@data$y<-data$y
  # max(data$closeness)
  # closenessAFMImage@data$h<-data$closeness*16/max(data$closeness)
  # displayIn3D(closenessAFMImage)
  # 
  # 
  # corenessAFMImage@data$x<-data$x
  # corenessAFMImage@data$y<-data$y
  # max(data$coreness)
  # corenessAFMImage@data$h<-data$coreness*16/max(data$coreness)
  # displayIn3D(corenessAFMImage)
  # 
  # ############################
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # allVertices<-data.table(vertexId=uniqueVerticesList, degree=edgeDegreeList)
  # # get-list of adjacent vertices with degree > 2 (can't remove if degree < 2)
  # allVertices<-allVertices[order(degree)]
  # 
  # allVertices[vertexId=="262161"]
  # allVertices[vertexId=="524305"]
  # 
  # 
  # adjacent_vertices(g, v=c("262161"), mode = c("all"))
  # 
  # edgeDegreeList$name
  # 
  # "524304"
  # [1] "262161" "524305"
  # 
  # 
  # 
  # 
  # av<-adjacent_vertices(g, "7340033")[["7340033"]]
  # unique(av$name)[1]
  # 
  # 
  # unlist(av)$name
  # 
  # print(av[[]])
  # 
  # names(av[]$name)
  # 
  # 
  # 
  # graph_from_edgelist
  # 
  # get.vertex.attribute(av)
  # 
  # #   
  # #   
  # #   
  # #   
  # #   #   
  # #   #   subEdgeList<-edgeList[1:1000]
  # #   #   
  # #   #   sg<-induced_subgraph(g, 1:250)
  # #   #   sg
  # #   #   plot.igraph(sg)
  # #   #   plot(sg, layout= layout_on_sphere, vertex.color="green")
  # #   #   
  # #   #   tkid <- tkplot(sg) #tkid is the id of the tkplot that will open
  # #   #   l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
  # #   #   plot(sg, layout=l)
  # #   #   
  # #   #   bt<-betweenness(sg)
  # #   #   bt[bt!=0]
  # #   #   
  # #   
  # #   # betweeness centrality network
  # # #   library(AFM)
  # # #   library(data.table)
  # # #   
  # # #   bt<-betweenness(g)
  # # #   bt[bt!=0]
  # # #   go<-bt
  # # #   nm<-go
  # # #   nm<-c(nm, rep(0, length(AFMImage@data$h)-length(nm)))
  # # #   nm[nm>0]<-250
  # # #   nm[is.infinite(nm)]<-500
  # # #   cAFMImageh<-copy(AFMImage@data$h)
  # # #   nm[is.na(nm)]<-(cAFMImageh[is.na(nm)])/10
  # # #   
  # # #   
  # # #   cl<-closeness(g)
  # # #   cl[cl!=0]
  # # #   cl<-cl*250/max(cl)
  # # #   go<-cl
  # # #   nm<-go
  # # #   nm<-c(nm, rep(0, length(AFMImage@data$h)-length(nm)))
  # # #   nm[nm>0]<-250
  # # #   nm[is.infinite(nm)]<-500
  # # #   cAFMImageh<-copy(AFMImage@data$h)
  # # #   nm[is.na(nm)]<-(cAFMImageh[is.na(nm)])/10
  # # #   
  # # #   
  # # #   dg<-degree(g)
  # # #   dg[dg!=0]
  # # #   dg<-dg*64/max(dg)
  # # #   go<-dg
  # # #   nm<-go
  # # #   nm<-c(nm, rep(0, length(AFMImage@data$h)-length(nm)))
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   # create a 128 pixels by 128 pixels AFM image
  # # #   Lines=AFMImage@lines
  # # #   Samplesperline=AFMImage@samplesperline
  # # #   fullfilename="BetweennessNetworkAFMImage"
  # # #   # the size of scan is 128 nm
  # # #   ScanSize=AFMImage@scansize
  # # #   # the heights is a normal distribution in nanometers
  # # #   
  # # #   
  # # #   
  # # #   scanby<-ScanSize/Samplesperline
  # # #   endScan<-ScanSize*(1-1/Samplesperline)
  # # #   BetweennessNetworkAFMImage<-AFMImage(
  # # #     data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
  # # #                       y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
  # # #                       h = nm),
  # # #     samplesperline = Samplesperline, lines = Lines, 
  # # #     vscansize = ScanSize, hscansize = ScanSize, scansize = ScanSize, 
  # # #     fullfilename = fullfilename )
  # # #   
  # # #   displayIn3D(BetweennessNetworkAFMImage)
  # # #   
  # # #   ####################################################################################
  # # #   AFMImage<-importFromNanoscope("C:\\Users\\nobody\\Desktop\\20140521.001_flatten.txt")
  # # #   AFMImage<-extractAFMImage(AFMImage,15,15,48)
  # # #   edgeList=c()
  # # #   
  # # #   AFMImage@data$h<-AFMImage@data$h*10
  # # #   AFMImage@data$h[AFMImage@data$h<29.9]<-0
  # # #   displayIn3D(AFMImage)
  # # #   BetweennessNetworkAFMImage<-AFMImage
  # # #   
  # # # #   
  # # # #   hello<-function(AFMImage, counter)    {
  # # # #     
  # # # #     for (x in seq(1: AFMImage@samplesperline)) {
  # # # #       for (y in seq(1: AFMImage@lines)) {
  # # # #         if(isAdjacentToBetterVertex(AFMImage,x,y)) {
  # # # #           counter<-counter+1
  # # # #           vertexId<-getVertexId(AFMImage,x,y) 
  # # # #           AFMImage@data$h[vertexId]<-0
  # # # #         }
  # # # #       }
  # # # #     }
  # # # #     return(counter)
  # # # #     
  # # # #   }
  # # #   
  # # #   ####### several times
  # # #   mycoun<-0
  # # #   while(mycoun<10) {
  # # # 
  # # #     edgeList=c()
  # # #     counter<-0
  # # #     for (x in seq(1: BetweennessNetworkAFMImage@samplesperline)) {
  # # #       for (y in seq(1: BetweennessNetworkAFMImage@lines)) {
  # # #         currentVertexId<-getVertexId(BetweennessNetworkAFMImage,x,y)
  # # #         if (existsEdge(BetweennessNetworkAFMImage, currentVertexId)) {
  # # #           edgeList<-c(edgeList, getSurroundingVertexesList(BetweennessNetworkAFMImage,x,y))
  # # #           counter<-counter+1
  # # #         }
  # # #       }
  # # #     }
  # # #     print(paste("Vertexes= ", counter))
  # # #     el <- matrix(edgeList, nc = 2, byrow = TRUE) 
  # # #     
  # # #     g<-graph_from_edgelist(el, directed=FALSE)
  # # #     
  # # #     vertexId<-2208
  # # #     ladv<-adjacent_vertices(g, vertexId)
  # # #     
  # # #     listOfVertices<-unlist(ladv[[1]][seq(1,length(ladv[[1]]),by=2)])
  # # #     listOfAllVertices<-unique(edgeList)
  # # #     
  # # #     
  # # #     for(vertice in listOfAllVertices) {
  # # #       
  # # #       ladv<-adjacent_vertices(g, vertice)
  # # #       newlistOfVertices<-ladv[[1]][seq(1,length(ladv[[1]]),by=2)]
  # # #       if (compareEqual(listOfAllVertices,newlistOfVertices, round=function(x) { x })) print("found")
  # # #     }
  # # #     
  # # #     
  # # #    # plot.igraph(g)
  # # #     
  # # #     dg<-degree(g)
  # # #     dg[dg!=0]
  # # #     print(max(dg))
  # # #     dg<-dg*64/max(dg)
  # # #     go<-dg
  # # #     nm<-go
  # # #     nm<-c(nm, rep(0, length(BetweennessNetworkAFMImage@data$h)-length(nm)))
  # # #     
  # # #     # create a 128 pixels by 128 pixels AFM image
  # # #     Lines=AFMImage@lines
  # # #     Samplesperline=BetweennessNetworkAFMImage@samplesperline
  # # #     fullfilename="BetweennessNetworkAFMImage"
  # # #     # the size of scan is 128 nm
  # # #     ScanSize=BetweennessNetworkAFMImage@scansize
  # # #     # the heights is a normal distribution in nanometers
  # # #     
  # # #     
  # # #     
  # # #     scanby<-ScanSize/Samplesperline
  # # #     endScan<-ScanSize*(1-1/Samplesperline)
  # # #     BetweennessNetworkAFMImage<-AFMImage(
  # # #       data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
  # # #                         y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
  # # #                         h = nm),
  # # #       samplesperline = Samplesperline, lines = Lines, 
  # # #       vscansize = ScanSize, hscansize = ScanSize, scansize = ScanSize, 
  # # #       fullfilename = fullfilename )
  # # #     
  # # #     displayIn3D(BetweennessNetworkAFMImage)
  # # # 
  # # #     newc<-0
  # # #     for (x in seq(1: BetweennessNetworkAFMImage@samplesperline)) {
  # # #       for (y in seq(1: BetweennessNetworkAFMImage@lines)) {
  # # #         if(isAdjacentToBetterVertex(BetweennessNetworkAFMImage,x,y)) {
  # # #           newc<-newc+1
  # # #           vertexId<-getVertexId(BetweennessNetworkAFMImage,x,y) 
  # # #           BetweennessNetworkAFMImage@data$h[vertexId]<-0
  # # #         }
  # # #       }
  # # #     }
  # # # 
  # # #     print(newc)
  # # #     displayIn3D(BetweennessNetworkAFMImage)
  # # # 
  # # #     mycoun<-mycoun+1
  # # #   }
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   
  # # #   # max-cliques
  # # #   library(AFM)
  # # #   library(data.table)
  # # #   
  # # #   # create a 128 pixels by 128 pixels AFM image
  # # #   Lines=512
  # # #   Samplesperline=512
  # # #   fullfilename="BetweennessNetworkAFMImage"
  # # #   # the size of scan is 128 nm
  # # #   ScanSize=512
  # # #   # the heights is a normal distribution in nanometers
  # # #   mc<-max_cliques(g)
  # # #   nm<-mc
  # # #   
  # # #   scanby<-ScanSize/Samplesperline
  # # #   endScan<-ScanSize*(1-1/Samplesperline)
  # # #   BetweennessNetworkAFMImage<-AFMImage(
  # # #     data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
  # # #                       y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
  # # #                       h = nm),
  # # #     samplesperline = Samplesperline, lines = Lines, 
  # # #     vscansize = ScanSize, hscansize = ScanSize, scansize = ScanSize, 
  # # #     fullfilename = fullfilename )
  # # #   
  # # #   BetweennessNetworkAFMImage@data$h[is.na(BetweennessNetworkAFMImage@data$h)]<-0
  # # #   BetweennessNetworkAFMImage@data$h[is.infinite(BetweennessNetworkAFMImage@data$h)]<-5000
  # # #   BetweennessNetworkAFMImage@data$h[BetweennessNetworkAFMImage@data$h>0]<-20
  # # #   
  # # #   displayIn3D(AFMImage)
  # # #   displayIn3D(BetweennessNetworkAFMImage)
  # # #   
  # # #   
  # # #   min(BetweennessNetworkAFMImage@data$h)
  # # #   
  # # #   max(BetweennessNetworkAFMImage@data$h)
  # # #   
  # # #   
  # # #   rglplot(sg, layout=layout_on_grid(g, dim = 3))
  # # #   
  # # #   layout_on_grid(g, width = 0, height = 0, dim = 2)
  # # #   rglplot(g, layout=layout_on_grid(g, dim = 3))
  # # #   
  # # #   plot.igraph(g)
  # # #   #crash plot(g, layout=layout_with_kk, vertex.color="green")
  # # #   mc<-max_cliques(g, min=5)
  # # #   bt<-betweenness(g)
  # # #   g[2]
  # # #   
  # #   
})
