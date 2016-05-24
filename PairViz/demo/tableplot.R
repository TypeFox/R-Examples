library(PairViz) 	
		
tabc <- apply(HairEyeColor, c(1, 2), sum)


dev.new()
par(mar=c(3,3,1,1))
par(cex=.6,mgp=c(2, -.5, 0))
table_plot(sqrt(tabc),sqrt(tabc))
# this table plot has cells with widths and heights proportional to the square root of cell counts.

tabp <- prop.table(tabc,2)

table_plot(apply(tabc,2,sum),tabp) # make cell widths proportional to margin totals, heights to conditional prob

cols <- 2:5
table_plot(apply(tabc,2,sum),tabp, yjust="bottom",col=cols,yruler=c("left","right")) # add colours, rulers and  bottom-justify

# The result is similar to the mosaic, without the mosaic effect of equalizing gaps. In the table version the rectangles line up across rows, so comparing heights, ie. conditional probs is easier.

o <- hpaths(1:4)[2,]
table_plot(apply(tabc,2,sum)[o],tabp[,o], yjust="bottom",col=cols,yruler=c("left","right"))
# Permutes the columns so all pairs of columns can be compared. In the second permutation can easily see that   p(black|blue eyes)> p(black|green eyes)


#mosaicplot(t(tabc)[,nrow(tabc):1],col=rev(cols),main="")
# mosaic- good for seeing deviations from independence. hard to compare conditional probs,
# except for those in the bottom and top rows. 


# 
# This data was reported in "A handbook of statistical analyses using SAS" Geoff Der and Brian Everitt page 307.
# A large number of people in the UK were surveyed and asked which of 13 characteristics they would associate
# with the nationals of the UK's partner countries in the Eurpoean Community. Entries are
# percents of respondants,


x <- c(37,29,21,19,10,10,8,8,6,6,5,2,1,7,14,8,9,27,7,3,7,3,23,12,1,3,
30,12,19,10,20,7,12,6,6,13,10,1,2,9,14,4,6,27,12,2,13,26,16,29,6,25,
1,7,1,16,30,3,10,9,5,11,22,2,27,
5,4,2,2,15,2,0,13,24,1,28,4,6,
4,48,1,12,3,9,2,11,41,1,38,8,8)

x <- matrix(x,ncol=7)
colnames(x) <- c("france","spain","italy","uk","ireland","holland","germany")
rownames(x) <- c("stylish","arrogant","sexy","devious","easygoing","greedy","cowardly","boring","efficient","lazy","hardworking","clever","courageous")



# for column colours
#cols <- matrix(rainbow(ncol(x)),nrow(x),ncol(x),byrow=TRUE)

# for row colours
cols <- matrix(rainbow(nrow(x)),nrow(x),ncol(x))

# First, get an ordering of countries that uses clustering to place similar countries adjacently, and do the same to group characteristics

library(gclus)
c.ord <- order.hclust(-dist(t(x)),method="average")
r.ord <- order.hclust(-dist(x),method="average")

dev.new(width=4,height=4)
par(mar=c(3,3.5,1,1))
par(cex=.6,mgp=c(2, -.5, 0))

table_plot(sqrt(x)[r.ord,c.ord],sqrt(x)[r.ord,c.ord],col=cols,yruler="center")
# Ireland's ratings are quite similar to those of the UK, less so to its other neighbour Spain

# We use weighted_hpaths to give the next "lowest-weight" hamiltonian amoung those in the hpaths decomposition, where here "lowest-weight" is the sum of distances  between adjacent countries.
h <- weighted_hpaths(dist(t(x)),c.ord,cycle=FALSE)
o <- h[2,]
dev.new(width=4,height=4)
par(mar=c(3,3.5,1,1))
par(cex=.6,mgp=c(2, -.5, 0))

table_plot(sqrt(x[r.ord,o]),sqrt(x[r.ord,o]),col=cols,yruler="center")

# Note with cycle=FALSE hpaths obtains a hamiltonian path decomposition, which is not exact (has a few duplicate edges). These duplicate edges appear here- Ireland neighbours UK in both displays shown.
