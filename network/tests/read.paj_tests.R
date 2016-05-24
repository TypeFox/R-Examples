# test for reading pajek formatted files
require(testthat)
require(network)


# test for case of verticse, but no edges/arcs
tmptest<-tempfile()
cat("*Vertices          2
        1    1231062
        2    1231095
*Arcs
*Edges
",file=tmptest)
noEdges<-read.paj(tmptest)
expect_equal(network.size(noEdges),2)
expect_equal(network.edgecount(noEdges),0)


# check arcs vs edges parsing

# arcs only

tmptest<-tempfile()
cat("*Vertices          3
1   'A' 
2   'B'
3   'C' 
*Arcs
1 2 1
1 3 1
",file=tmptest)
arcsOnly<-read.paj(tmptest)
expect_true(is.directed(arcsOnly))
expect_equal(network.edgecount(arcsOnly),2)

# edges only
tmptest<-tempfile()
cat('*Vertices      9
1 "1"    0.3034    0.7561
2 "2"    0.4565    0.6039
3 "3"    0.4887    0.8188
4 "4"    0.5687    0.4184
5 "5"    0.3574    0.4180
6 "6"    0.7347    0.2678
7 "7"    0.9589    0.3105
8 "8"    0.8833    0.1269
9 "9"    0.7034    0.0411
*Arcs
*Edges
1      2       1
1      3       1
2      3       1
2      4       1
2      5       1
4      5       1
4      6       1
6      7       1
6      8       1
6      9       1
7      8       1
8      9       1 
',file=tmptest)
edgesOnly<-read.paj(tmptest)
expect_false(is.directed(edgesOnly))
expect_equal(network.edgecount(edgesOnly),12)


# both arcs and edges
# network will be directed, each *edges record will create one arc in each direction
tmptest<-tempfile()
cat("*Vertices          4
1   'A' 
2   'B' 
3   'C'
4   'D'
*Arcs
1 2 1
1 3 1
*Edges
3 4 1
",file=tmptest)
arcsNEdges<-read.paj(tmptest)
expect_true(is.directed(arcsNEdges))
expect_equal(network.edgecount(arcsNEdges),4)
as.matrix(arcsNEdges)


# ----- error testing
tmptest<-tempfile()
cat("*Vertices          2
1   'A' 
2   'B' 
*Arcs
1 
",file=tmptest)
expect_error(read.paj(tmptest),regexp = 'does not appear to have the required')

tmptest<-tempfile()
cat("*Vertices          2
1   'A' 
2   'B' 
*Arcs
1 A 1
",file=tmptest)
expect_error(read.paj(tmptest),regexp = 'contains non-numeric or NA values')

tmptest<-tempfile()
cat("*Vertices          2
1   'A' 
2   'B' 
*Arcs
1 2.5 1
",file=tmptest)
expect_error(read.paj(tmptest),regexp = 'contains non-integer values')


# check vertex graphic attribute fill-in
tmptest<-tempfile()
cat("*Vertices          4
1   'A' 0 0 0 box
2   'B' 0 0 0
3   'C' 0 0 0 
4   'D' 0 0 0 ellipse
*Arcs
1 2 1
1 3 1
",file=tmptest)
fillIn<-read.paj(tmptest)
expect_equal(fillIn%v%'shape',c('box','box','box','ellipse'))


# test stuff in file comments
########## but multirelational ############ only ~200  nodes 
#GulfLDays.net 
#GulfLMonths.net
#GulfLDow.net 
#gulfAllDays.net     #GulfADays.zip
#gulfAllMonths.net   #GulfAMonths.zip
#LevantDays.net 
#LevantMonths.net
#BalkanDays.net 
#BalkanMonths.net 

#arcs and edges both present   search for " #these have both arc and edge lines " or "URL has a net file"
#Graph drawing competition page (GD)
#C95,C95,B98,A99,C99,A99m


#things to do:
#handle ragged array .net files like "CSphd.net"     DONE!!
#handel two mode networks                            DONE!!
#handle mix of edges and arcs                        DONE!!
#handle multirelational pajek files

#issue with read.table and number.cols and fill...SanJuanSur_deathmessage.net has one row with 8 all the rest (including the first 5 have 5)



# this file has character encoding issues
scotland<-tempfile('scotland',fileext='.zip')
download.file('http://vlado.fmf.uni-lj.si/pub/networks/data/esna/scotland.zip',scotland)
scotpaj<-tempfile('Scotland',fileext='.paj')
cat(readLines(unz(scotland,'Scotland.paj')),sep='\n',file = scotpaj)
scotproj<-read.paj(scotpaj)

# produces two element list, containing networks and partitions
expect_equal(names(scotproj),c("networks","partitions"))
expect_equal(network.size(scotproj[[1]][[1]]),244)
expect_equal(names(scotproj$partitions),c("Affiliation.partition.of.N1.[108,136]","Industrial_categories.clu"))


A95net<-read.paj("http://vlado.fmf.uni-lj.si/pub/networks/data/GD/gd95/A95.net")
expect_equal(network.size(A95net),36)
expect_equal(network.vertex.names(A95net)[1:5],c("MUX","INSTRUCTION BUFFER (4 x 16)", "RCV","DRV","ROM REG"))

# test reading a .paj project file
bkfratProj<-read.paj('http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/bkfrat.paj')

# should have two networks
expect_equal(sapply(bkfratProj,class),c('network','network'),check.attributes=FALSE)
# .. with wierd names
expect_equal(names(bkfratProj),c('UciNet\\BKFRAT.DAT : BKFRAB','UciNet\\BKFRAT.DAT : BKFRAC'))
# and 58 vertices
expect_equal(sapply(bkfratProj,network.size),c(58,58),check.attributes=FALSE)
expect_equal(sapply(bkfratProj,network.edgecount),c(1934,3306),check.attributes=FALSE)

#check edge values and attribute naming
expect_equal((bkfratProj[[1]]%e%"UciNet\\BKFRAT.DAT : BKFRAB")[1900:1934],c(1, 1, 1, 5, 2, 4, 2, 1, 3, 1, 3, 1, 2, 5, 1, 1, 1, 2, 1, 2, 2, 1, 6, 2, 1, 2, 2, 1, 1, 1, 1, 3, 3, 1, 1))

# check vert attrs
expect_equal(list.vertex.attributes(bkfratProj[[1]]),c('na','vertex.names','x','y','z'))

# check network attrs
expect_equal(bkfratProj[[1]]%n%'title',"UciNet\\BKFRAT.DAT : BKFRAB")
expect_equal(bkfratProj[[2]]%n%'title',"UciNet\\BKFRAT.DAT : BKFRAC")

# check loop flagging

tmptest<-tempfile()
cat("*Vertices          2
1   'A' 
2   'B' 
*Arcs
1 1 1
",file=tmptest)
loopTest<-read.paj(tmptest,verbose=TRUE)
expect_true(has.loops(loopTest))

# check edge.name argument

tmptest<-tempfile()
cat("*Vertices          2
1   'A' 
2   'B' 
*Arcs
1 1 1
",file=tmptest)
loopTest<-read.paj(tmptest,verbose=TRUE,edge.name='weight')
expect_equal(list.edge.attributes(loopTest),c('na','weight'))

# the rest of these will take longer, so including in opttest block so won't run on CRAN
require(statnet.common)
opttest(testvar = "ENABLE_statnet_TESTS",{


#  ----- examples from http://vlado.fmf.uni-lj.si/pub/networks/doc/ECPR/08/ECPR01.pdf ---

GraphSet<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/GraphSet.net')
expect_true(is.directed(GraphSet))
expect_equal(network.edgecount(GraphSet),27)
# network contains some repeated edges
expect_true(is.multiplex(GraphSet))
expect_equal(network.vertex.names(GraphSet),letters[1:12])

Tina<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/TinaSet.net')

# arcslist
GraphList<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/GraphList.net')  
# http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/TinaList.net  # arcslist

# matrix
GraphMat <-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/GraphMat.net') 
expect_equal(network.vertex.names(GraphMat),letters[1:12])
# check that edge attribute created and parsed correctly
expect_equal(as.matrix(GraphMat,attrname='GraphMat')[3,7],2)

# partition
TinaPaj<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Tina.paj')
expect_equal(class(TinaPaj$partitions),'data.frame')
expect_equal( TinaPaj$partitions[,1],c(2,1,2,2,2,2,2,2,3,3,3),use.names=FALSE)
expect_true(is.network(TinaPaj$networks$Tina))

# --- crude timing info --
# by default timing info should be added as attribute
timetest<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Time.net')
expect_equal(timetest%e%'pajekTiming',c("[7]","[6-8]"))
expect_equal(timetest%v%'pajekTiming',c("[5-10,12-14]", "[1-3,7]", "[4-*]"))
expect_true(setequal(list.vertex.attributes(timetest),c('na','pajekTiming','vertex.names'))) # no x or y
expect_true(setequal(list.edge.attributes(timetest),c('na','pajekTiming','Time')))

# test converting to networkDynamic format
timetestNd<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Time.net',time.format='networkDynamic')
expect_equal(class(timetestNd),c('networkDynamic','network'))
# check that activiy matrices are built as expected
expect_equal(get.vertex.attribute(timetestNd,'active',unlist=FALSE),list(structure(c(5, 12, 11, 15), .Dim = c(2L, 2L)), structure(c(1, 7, 4, 8), .Dim = c(2L, 2L)), structure(c(4, Inf), .Dim = 1:2)))
expect_equal(get.edge.attribute(timetestNd,'active',unlist=FALSE),list(structure(c(7, 8), .Dim = 1:2), structure(c(6, 9), .Dim = 1:2)))

# read a *big* one http://vlado.fmf.uni-lj.si/pub/networks/data/CRA/Days.zip
# 1.3 Mb, 13k vertices, 256K lines. 
#  days<-tempfile('days',fileext='.zip')
#  download.file('http://vlado.fmf.uni-lj.si/pub/networks/data/CRA/Days.zip',days)
#  terrorTerms<-read.paj(unz(days,'Days.net'),verbose=TRUE,time.format='networkDynamic',edge.name='count')



# multiple networks
sampson<-read.paj('http://vlado.fmf.uni-lj.si/pub/networks/pajek/data/Sampson.net')  
lapply(sampson,class)  # for some reason it is a formula?
expect_equal(names(sampson$networks),c("SAMPLK1", "SAMPLK2", "SAMPLK3",  "SAMPDLK",  "SAMPES","SAMPDES","SAMPIN","SAMPNIN","SAMPPR","SAMNPR"))

# multiple networks in arcslist format
# sampsonL<-read.paj('http://vlado.fmf.uni-lj.si/pub/networks/pajek/data/SampsonL.net') 

# two-mode
sandi<-read.paj('http://vlado.fmf.uni-lj.si/pub/networks/data/2mode/sandi/sandi.net')  
expect_true(is.bipartite(sandi))
expect_equal(sandi%n%'bipartite',314)
Davis<-read.paj('http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Davis.paj')         # two-mode
expect_equal(Davis$networks[[1]]%n%'bipartite',18)

# lots of edge and vertex attributes
A96<-read.paj('http://vlado.fmf.uni-lj.si/pub/networks/data/GD/gd96/A96.net')
expect_equal(network.size(A96),1096)
expect_equal(list.vertex.attributes(A96),c("bw","fos","na","shape","vertex.names", "x","x_fact","y","y_fact"))   # note no z attribute
expect_equal(head(A96%v%'shape'),c("box","ellipse", "ellipse", "ellipse", "ellipse", "ellipse"))
# check edge attribute parsing
expect_equal(list.edge.attributes(A96),c("A96", "fos", "l" ,  "lr",  "na",  "s",   "w"  ))
# l is the only one with unique values
expect_equal(head(A96%e%'l'),c("a", "s","n","r","s","t"))

})  # end of non-cran tests



# temporal versions http://vlado.fmf.uni-lj.si/pub/networks/data/KEDS/KEDS.htm

# temporal events data (not supported)
# http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Time.tim
# http://vlado.fmf.uni-lj.si/vlado/podstat/AO/net/Friends.tim
