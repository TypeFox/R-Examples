plotProfile <-
function(data,rank){

# data: dataForStabPlot object
# rank: rank of the gene set
# gene.perm: TRUE for plotting gene perm data


    rv <- data$rank.vector
    rnk <- match(rank,rv)
    geneset <- data$gene.sets.stb[rnk]
    cat(geneset)
    ind <- match(data$gene.sets.stb,colnames(data$gene.sets))
    ord <- order(ind)
    rank.vector <- data$rank.vector[ord]
    rank <- match(rank, rank.vector)
    pos.data <- (data$positive.prof.up)[[rank]]
    perm.data <- matrix(0,data$perm.number,length(pos.data))
    for(i in 1:data$perm.number){
           perm.data[i,]= data$perm.data[[i]][[rank]]}
    
    
    quantiles <-  apply(perm.data,2,quantile,probs=c(0,0.10,0.25,0.50,0.75,0.90,1))
    total.data <- rbind(pos.data,quantiles)
    df <- data.frame(Real.data= total.data[1,],Q1 = total.data[2,],Q2=total.data[3,], Q3 = total.data[4,], Q4=total.data[5,],Q5=total.data[6,],Q6=total.data[7,],Q7=total.data[8,])
    plot(df$Real.data,lwd=3,type="l",xlab="Ordered gene list",ylab="Running scores",font=2,font.lab=2,family="sans",ylim=c(min(as.numeric(perm.data),as.numeric(pos.data)),max(as.numeric(perm.data),as.numeric(pos.data))),xlim=c(1,length(pos.data)+ (length(pos.data))/5))
    lines(df$Q1,col="red",lty=2,lwd=3)
    lines(df$Q2,col="blue",lty=2,lwd=3)
    lines(df$Q3,col="green",lty=2,lwd=3)
    lines(df$Q4,col="purple",lty=2,lwd=3)
    lines(df$Q5,col="magenta",lty=2,lwd=3)
    lines(df$Q6,col="pink",lty=2,lwd=3)
    lines(df$Q7,col="grey",lty=2,lwd=3)
    legend("topright",c("True data","0.00","0.10","0.25","0.50","0.75","0.90","1.00"),lty=c(1,2,2,2,2,2,2,2),col=c("black","red","blue","green","purple","magenta","pink","grey"),bty="n",lwd=3,cex=.70)}
