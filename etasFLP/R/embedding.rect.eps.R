embedding.rect.eps <-
function(x,y,cycle=FALSE,eps=0){
			eps.x	=diff(range(x))*eps/2
			eps.y	=diff(range(y))*eps/2
			rx	=c(min(x)-eps.x,max(x)+eps.x)
			ry	=c(min(y)-eps.y,max(y)+eps.y)
   
			rect	=c(rx[2],rx[1],rx[1],rx[2],rx[2],ry[2],ry[2],ry[1],ry[1],ry[2])
			rect	=matrix(rect,5,2)
			if(!cycle)rect=rect[1:4,]
			return(rect)
							}
