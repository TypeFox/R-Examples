
SHAH.plot =function(out.BUUHWE_2D, Evol=TRUE){
BUUHWE_2D.plot(out.BUUHWE_2D, Evol)
}

BUUHWE_2D.plot <- function(out.BUUHWE_2D, Evol=TRUE) {

if(Evol){
# ===========================INITIALISATION ===========================

#dimension
m =out.BUUHWE_2D$d[1]  
n =out.BUUHWE_2D$d[2]

# graph parameters
oldpar = par(no.readonly = TRUE)
par(mfrow=c(m,n), oma=c(0,0,0,0), mar=c(0,0,0,0))

# number of edges
noe <- 2 * m * n - m - n    


# NEW EDGES
edges <- matrix(0, 2*(n*m-1), 2)

	# Columns: From, to.
grid=matrix(1:(n*m),nrow=m)
tgrid=t(grid)


for (k in 1:(m*n-1)){
    edges[k,] =c(grid[k],grid[k+1])
    edges[(m*n-1)+k,] =c(tgrid[k],tgrid[k+1]) 
	}

RM_V = seq(m,m*n-1, by=m)
RM_H = seq(n,m*n-1, by=n)	+(m*n-1)
RM=c(RM_V,RM_H)
edges =edges[-RM,] 

## modify order of the edges so as to coincide with Buuhwe_2D_Step

lengthH = (n*m-1)-length(RM_H)
selectHinput= edges[(nrow(edges)-lengthH+1):nrow(edges),1]
reorderingH = order(selectHinput)
lengthV = (n*m-1)-length(RM_V)
edges = edges[c((1:lengthV), lengthV+reorderingH ) ,]

	
 
#---------- nodes -----------------------------------------------------
nodeIn = matrix(NA,ncol=2, nrow= nrow(edges))
nodeOut = matrix(NA,ncol=2, nrow= nrow(edges))

for(i in 1:nrow(edges)){
nodeIn[i,] = c(row(grid)[grid==edges[i,1]],col(grid)[grid==edges[i,1]])
nodeOut[i,] = c(row(grid)[grid==edges[i,2]],col(grid)[grid==edges[i,2]])
}


### identifying pixel reduction

edges_select= t(out.BUUHWE_2D$decomp.hist[1,,])

# Pixel reduction at step 1
pixelIn = c(row(grid)[grid==edges_select[1,1]],col(grid)[grid==edges_select[1,1]])
pixelOut =  c(row(grid)[grid==edges_select[1,2]],col(grid)[grid==edges_select[1,2]])


#### draw 

plot(1:3,1:3, type='n', xaxt='n',yaxt='n',bty='n',
    xlab='',ylab='', xlim=c(0,n+1),ylim=c(0,m+1))
    
for(i in 1:nrow(edges)){
segments(nodeIn[i,2],nodeIn[i,1],nodeOut[i,2],nodeOut[i,1])
points(nodeIn[i,2],nodeIn[i,1], pch=21,cex =4, bg='white' )
points(nodeOut[i,2],nodeOut[i,1],pch=21,cex =4, bg='white' )
}

points(pixelIn[2],pixelIn[1], pch=21,cex =4, bg='green' )
points(pixelOut[2],pixelOut[1],pch=21,cex =4, bg='green' )

for(i in 1:nrow(edges)){
text(nodeIn[i,2],nodeIn[i,1], labels= edges[i,1],offset=0)
text(nodeOut[i,2],nodeOut[i,1], labels= edges[i,2],offset=0)
}


#===================== FURTHER STEPS =====================================

#edges_select= t(out.BUUHWT_2D$decomp.hist[1,,])

for (k in 1:nrow(edges_select)){

edges[edges == edges_select[k,2]] <- edges_select[k,1]

# Pixel reduction at step k+1
if (k< nrow(edges_select)){
pixelIn = c(row(grid)[grid==edges_select[k+1,1]],col(grid)[grid==edges_select[k+1,1]])
pixelOut =  c(row(grid)[grid==edges_select[k+1,2]],col(grid)[grid==edges_select[k+1,2]])
}


nodeIn = matrix(NA,ncol=2, nrow= nrow(edges))
nodeOut = matrix(NA,ncol=2, nrow= nrow(edges))
for(i in 1:nrow(edges)){
nodeIn[i,] = c(row(grid)[grid==edges[i,1]],col(grid)[grid==edges[i,1]])
nodeOut[i,] = c(row(grid)[grid==edges[i,2]],col(grid)[grid==edges[i,2]])
}

plot(1:3,1:3, type='n', xaxt='n',yaxt='n',bty='n',xlab='',ylab='', ylim=c(0,m+1),xlim=c(0,n+1))

for(i in 1:nrow(edges)){
segments(nodeIn[i,2],nodeIn[i,1],nodeOut[i,2],nodeOut[i,1])
points(nodeIn[i,2],nodeIn[i,1], pch=21,cex =4, bg='white' )
points(nodeOut[i,2],nodeOut[i,1],pch=21,cex =4, bg='white' )
}

if (k< nrow(edges_select)){
points(pixelIn[2],pixelIn[1], pch=21,cex =4, bg='green' )
points(pixelOut[2],pixelOut[1],pch=21,cex =4, bg='green' )
}
for(i in 1:nrow(edges)){
text(nodeIn[i,2],nodeIn[i,1], labels= edges[i,1],offset=0)
text(nodeOut[i,2],nodeOut[i,1], labels= edges[i,2],offset=0)
}



}


par(oldpar)}
else{ #================================================================
#======================================================================
# ===========================INITIALISATION ===========================

#dimension
m =out.BUUHWE_2D$d[1]  
n =out.BUUHWE_2D$d[2]

# graph
oldpar = par(no.readonly = TRUE)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(0,0,0,0))

# number of edges 
noe <- 2 * m * n - m - n    


# NEW EDGES
edges <- matrix(0, 2*(n*m-1), 2)

	# Columns: From, to.
grid=matrix(1:(n*m),nrow=m)
tgrid=t(grid)


for (k in 1:(m*n-1)){
    edges[k,] =c(grid[k],grid[k+1])
    edges[(m*n-1)+k,] =c(tgrid[k],tgrid[k+1]) 
	}

RM_V = seq(m,m*n-1, by=m)
RM_H = seq(n,m*n-1, by=n)	+(m*n-1)
RM=c(RM_V,RM_H)
edges =edges[-RM,] 	

## modify order of the edges so as to coincide with Buuhwe_2D_Step

lengthH = (n*m-1)-length(RM_H)
selectHinput= edges[(nrow(edges)-lengthH+1):nrow(edges),1]
reorderingH = order(selectHinput)
lengthV = (n*m-1)-length(RM_V)
edges = edges[c((1:lengthV), lengthV+reorderingH ) ,]



#---------- nodes -----------------------------------------------------
nodeIn = matrix(NA,ncol=2, nrow= nrow(edges))
nodeOut = matrix(NA,ncol=2, nrow= nrow(edges))

for(i in 1:nrow(edges)){
nodeIn[i,] = c(row(grid)[grid==edges[i,1]],col(grid)[grid==edges[i,1]])
nodeOut[i,] = c(row(grid)[grid==edges[i,2]],col(grid)[grid==edges[i,2]])
}


### identifying pixel reduction

edges_select= t(out.BUUHWE_2D$decomp.hist[1,,])

# Pixel reduction at step 1
pixelIn = c(row(grid)[grid==edges_select[1,1]],col(grid)[grid==edges_select[1,1]])
pixelOut =  c(row(grid)[grid==edges_select[1,2]],col(grid)[grid==edges_select[1,2]])

#-----------------------------------------------------------------------------
#### Draw 
#-----------------------------------------------------------------------------

palette(topo.colors(nrow(edges_select)))     # six color rainbow

Bigwidth = abs(out.BUUHWE_2D$decomp.hist[3,1,])>abs(quantile(out.BUUHWE_2D$decomp.hist[3,1,],0.8))
width =rep(1,length(Bigwidth))
width[Bigwidth] = rep(5, length(width[Bigwidth]))

plot(1:3,1:3, type='n', xaxt='n',yaxt='n',bty='n',
    xlab='',ylab='', xlim=c(0,n+1),ylim=c(0,m+1))
    

segments(pixelIn[2],pixelIn[1],pixelOut[2],pixelOut[1],col=1,lwd=width[1] )




#===================== FURTHER STEPS =====================================

#edges_select= t(out.BUUHWT_2D$decomp.hist[1,,])

for (k in 1:nrow(edges_select)){

edges[edges == edges_select[k,2]] <- edges_select[k,1]

# Pixel reduction at step k+1
if (k< nrow(edges_select)){
pixelIn = c(row(grid)[grid==edges_select[k+1,1]],col(grid)[grid==edges_select[k+1,1]])
pixelOut =  c(row(grid)[grid==edges_select[k+1,2]],col(grid)[grid==edges_select[k+1,2]])
}


nodeIn = matrix(NA,ncol=2, nrow= nrow(edges))
nodeOut = matrix(NA,ncol=2, nrow= nrow(edges))
for(i in 1:nrow(edges)){
nodeIn[i,] = c(row(grid)[grid==edges[i,1]],col(grid)[grid==edges[i,1]])
nodeOut[i,] = c(row(grid)[grid==edges[i,2]],col(grid)[grid==edges[i,2]])
}

if (k< nrow(edges_select)){
segments(pixelIn[2],pixelIn[1],pixelOut[2],pixelOut[1],col=(k+1), lwd=width[k+1] )

#print(width[k+1])

}

palette("default") 

}


par(oldpar)}









}

  
 

				