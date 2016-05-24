peaks.on.som <- function(result, solution,
clusterPalette=rainbow(solution), cell.size = 0.9, mfrow=c(2,ceiling(n.errors/2)), new.order=1:solution){
    
    halfsize <- cell.size/2
    cords<-result$som$code.sum[,1:2]+1
    palette(clusterPalette)

    n.errors <- dim(result$measures.synthetic.peaks)[1]
    cluster.rer <- change.order.clusters(result$cluster.assignment[[solution]], new.order)


    par(mfrow = mfrow ,mar=c(1,1,5,1)+0.1, xpd=NA)
    for(i in 1:n.errors){
        plot(c(1,result$som.dim[2]),c(1,result$som.dim[1]), type="n",xaxt="n", yaxt="n",xlab="", ylab="")
        mtext(3,text=result$error.names[i], line=1, adj=0.1, cex=1)

        #turn arround coordinate to be consistent with shards
        rect(cords[,2]-halfsize,cords[,1]-halfsize,cords[,2]+halfsize,cords[,1]+halfsize,col=cluster.rer$cluster, border="transparent")
        unique.pos <- unique(result$som.pos.errors[i,,1:2])
        labels<-rep(NA, NROW(unique.pos))
        for(j in 1:NROW(unique.pos)){
           labels[j] <- paste(which(unique.pos[j,1]+result$som.dim[1]*unique.pos[j,2]==result$som.pos.errors[i,,1]+result$som.dim[1]*result$som.pos.errors[i,,2]), collapse=",")
        }
        #turn arround coordinate to be consistent with shards
        text(labels=labels,y=unique.pos[,1]+1, x=unique.pos[,2]+1)
        text(labels="*", x=result$som.pos.no.error[,2]+1 ,y=result$som.pos.no.error[,1]+1)
    }
}
