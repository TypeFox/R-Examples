paraplot <- function(para1, para2=NULL, para3=NULL, para4=NULL,
                      tickx=NULL, jitter=NULL, pch = 1:4, col=1:4, 
                      legpos=NULL, legtext, annotate=FALSE)
{
  K <- 4
  if(is.null(para4)) K <- K-1
  if(is.null(para3)) K <- K-1  
  if(is.null(para2)) K <- K-1  
 
  par(mfrow=c(1,1),mar=c(3,2,1,1))
  
  if(K==1){
    p <- nrow(para1)  
    index <- 1:p
    l <- min(para1[,2])
    u <- max(para1[,3])
    plot(para1[,1], index, xlim=c(l-abs(l)*0.1,u), ylim=c(1,p+1), axes=F, 
         pch=pch[1], col=col[1], xlab="", ylab="")
    segments(para1[,2],index, para1[,3], index, col=col[1])
    axis(2, at=NULL, labels=FALSE, tick = FALSE)
    if(is.null(tickx)) axis(1,round(seq(l,u,length.out=8),2))
    else axis(1,tickx)
    if(annotate){
      if(is.null(rownames(para1))) stop("Assign row names to the parameter object")
      else text(rep(l-abs(l)*0.15,p), 1:p, rownames(para1), cex=1.2)
    }
    legend(legpos[1], legpos[2], legtext, pch=pch[1], col=pch[1], lty=1, cex=1.2)  
  }
  else if(K>=2){
    mesg <- "Assign row names to the parameter object. \n Use the same row name if two parameter \n objects share the same parameter."

    if(K==2) {
      if(any(c(is.null(rownames(para1)),is.null(rownames(para2))))) stop(mesg)
      l <- min(c(min(para1),min(para2)))
      u <- max(c(max(para1),max(para2))) 
      
      para1 <- data.frame(para1); 
      para1 <- cbind(rownames(para1), para1); 
      colnames(para1) <- c("id","e1","l1","u1")
      
      para2 <- data.frame(para2);
      para2 <- cbind(rownames(para2), para2); 
      colnames(para2) <- c("id","e2","l2","u2") 

      para<- merge(para1, para2, by="id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))       
    }
    if(K==3) {
      if(any(c(is.null(rownames(para1)),is.null(rownames(para2)), 
               is.null(rownames(para3))))) stop(mesg)
      l <- min(c(min(para1),min(para2),min(para3)))
      u <- max(c(max(para1),max(para2),max(para3)))  
      
      para1<- data.frame(para1); id<- rownames(para1); 
      para1<- cbind(id, para1); colnames(para1) <- c("id","e1","l1","u1")
      
      para2<- data.frame(para2); id<- rownames(para2); 
      para2<- cbind(id, para2); colnames(para2) <- c("id","e2","l2","u2")
      
      para3<- data.frame(para3); id<- rownames(para3); 
      para3<- cbind(id, para3); colnames(para3) <- c("id","e3","l3","u3")  

      para<- merge(para1, para2, by = "id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))
      para<- merge(para, para3, by = "id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))
    }    
    if(K==4) {
      if(any(c(is.null(rownames(para1)),is.null(rownames(para2)), 
               is.null(rownames(para3)),is.null(rownames(para4))))) stop(mesg)
      l <- min(c(min(para1),min(para2),min(para3),min(para4)))
      u <- max(c(max(para1),max(para2),max(para3),max(para4)))  
      
      para1<- data.frame(para1); id<- rownames(para1); 
      para1<- cbind(id, para1); colnames(para1) <- c("id","e1","l1","u1")
      
      para2<- data.frame(para2); id<- rownames(para2); 
      para2<- cbind(id, para2); colnames(para2) <- c("id","e2","l2","u2") 
      
      para3<- data.frame(para3); id<- rownames(para3); 
      para3<- cbind(id, para3); colnames(para3) <- c("id","e3","l3","u3") 
      
      para4<- data.frame(para4); id<- rownames(para4); 
      para4<- cbind(id, para4); colnames(para4) <- c("id","e4","l4","u4")  
      
      para<- merge(para1, para2, by = "id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))  
      para<- merge(para, para3, by = "id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))
      para<- merge(para, para4, by = "id", all.x=TRUE, 
                   all.y=TRUE, suffixes = c(" "," "))  
    }  
    p<- nrow(para)
    index <- 1:p
        
    for(i in 1:p){
      for(j in 1:ncol(para)){
        if(is.na(para[i,j])) para[i,j]<- -99999
      }}
    if(is.null(jitter)) jitter <- p*0.01
    plot(para[,2], index, xlim=c(l-abs(l)*0.3,u), 
         ylim=c(1,p+max(1,jitter*K)), 
         col=col[1], pch=pch[1], axes=F, xlab="", ylab="")
    segments(para[,3], index, para[,4], index, col=col[1])      
    for(i in 2:K){
      points(para[,1+(i-1)*3+1], index+jitter*i, col=col[i], pch=pch[i])
      segments(para[,1+(i-1)*3+2], index+jitter*i, 
               para[,1+(i-1)*3+3], index+jitter*i, col=col[i])
    }
    axis(2, at=NULL, labels=FALSE, tick=FALSE)
    if(is.null(tickx)) axis(1,round(seq(l,u,length.out=8),2))
    else axis(1,tickx)
    if(annotate) text(rep(l-abs(l)*0.15,p), 1:p, para$id, cex=1.2)
    legend(legpos[1],legpos[2],legtext,
           pch=pch[1:K],col=pch[1:K],lty=rep(1,K), cex=1.2)  
  }
}



          

