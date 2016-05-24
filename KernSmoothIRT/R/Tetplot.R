Tetplot <-
function(OBJ,items, main, ...){


		if(missing(main)){main<--1}

		tethraedron <- function(listamatricedist,punti.colore="red",labels,points=TRUE,colors=colors,titl,OBJ=OBJ,...){

			
		
			a <- matrix(c(-sqrt(2)/sqrt(3),-sqrt(2)/3,-1/3,0,(2*sqrt(2))/3,-1/3,0,0,1),3,3,byrow=TRUE)
			puntini <- NULL
			for(j in 1:length(listamatricedist)){
				listamatricedist[[j]] -> lista
				punti <- matrix(0,nrow(lista),(ncol(lista)-1))
				for(i in 1:nrow(lista)) {
					b<-lista[i,2:4]
					p<-solve(a,b)
					punti[i,]<-p
				}

				puntini[[j]] <- punti 
			}

			segments3d(c(0,-sqrt(3)/(2*sqrt(2))), c(0,3/(2*sqrt(2))), c(0,0),size=1.2)
			segments3d(c(-sqrt(3)/(2*sqrt(2)),-sqrt(3)/sqrt(2)), c(3/(2*sqrt(2)),0), c(0,0),size=1.2)
			segments3d(c(-sqrt(3)/sqrt(2),0), c(0,0), c(0,0),size=1.2)
			segments3d(c(0,-sqrt(3)/(2*sqrt(2))), c(0,1/(2*sqrt(2))), c(0,1),size=1.2)
			segments3d(c(-sqrt(3)/(2*sqrt(2)),-sqrt(3)/sqrt(2)), c(1/(2*sqrt(2)),0), c(1,0),size=1.2)
			segments3d(c(-sqrt(3)/sqrt(2),-sqrt(3)/(2*sqrt(2))), c(0,3/(2*sqrt(2))), c(0,0),size=1.2)
			segments3d(c(-sqrt(3)/(2*sqrt(2)),-sqrt(3)/(2*sqrt(2))), c(3/(2*sqrt(2)),1/(2*sqrt(2))), c(0,1),size=1.2)
		
			for(j in 1:length(puntini)){
				
				spheres3d(puntini[[j]],radius=.01,col=colors,...)

				alarm()



			}	
			
			text3d(0,0,0,labels[1])
			text3d(-sqrt(3)/(2*sqrt(2)),1,0,labels[3])
			text3d(-sqrt(3)/(2*sqrt(2)),1/(2*sqrt(2)),1,labels[4])
			text3d(-sqrt(3)/(sqrt(2)),0,0,labels[2])
			text3d(-.4,.4,1.2,titl,font=2)

			
	
		}
	
	plotit<-function(x,OBJ,main,...){
		


		Torank<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
		
		order<-rank(apply(Torank[,-c(1:3)],1,max))
		order<-max(order)+1-order	
	
	
		highest4<-sapply(order,function(yyy)ifelse(yyy>4,0,1))
		subset<-which(highest4==max(highest4))
		
		Normed<-Torank[subset,-c(1:3)]

		ToPlot<-t(Normed)
		tonamecols<-as.character(Torank[subset,2])
		tonamecols[which(tonamecols==-1)]='NA'
		
		
		colnames(ToPlot)<-tonamecols
		if(length(tonamecols)<4){return("Item has too few options")}
		
		ToPlot<-ToPlot/apply(ToPlot,1,sum)
		
		
		pts<-length(OBJ$evalpoints)
		onethird<-ceiling(pts/3)

		colors<-c(rep("red",onethird),rep("green",onethird),rep("blue",onethird))
	
		if(main==-1){main<-paste("Item : ",OBJ$itemlabels[[x]])}


		
		open3d()
		tethraedron(listamatricedist=list(ToPlot),labels=tonamecols,colors=colors,titl=main,OBJ=OBJ,...)
		
	
	}
	
	nada<-sapply(items,plotit,OBJ=OBJ,main,...)



	cat("Trait Levels")
	cat("\n")
	cat("Red: Low")
	cat("\n")
	cat("Green: Medium")
	cat("\n")
	cat("Blue: High")
	

}

