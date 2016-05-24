som.cluster <- function(result, maxc=result$maxc, solutions=2:maxc, new.order=NA){

for(centers in solutions){

        #Für jeden Datensatz den Cluster bestimmen, um Plot machen zu können
        if(any(is.na(new.order))){
            regions.cmeans <- som.regions.for.values(result=result, solution=centers)
        } else {
            if(length(solutions)!=1) stop("Plot with new order only possilbe for single solution")
            regions.cmeans <- som.regions.for.values(result=result, solution=centers, new.order=new.order)
        }

        shardsplot(result$som, data.or=cbind(result$measures.uniform[!result$na.rows,], regions.cmeans)  ,xaxt="n", yaxt="n")
        opar <- par(xpd = NA)
            legend(((result$som$ydim-1)/2)+1, result$som$xdim+0.2, pt.bg = rainbow(max(regions.cmeans)), xjust = 0.5, yjust = 0, legend = LETTERS[1:max(regions.cmeans)], pch = 21, ncol=5, col="black", cex=0.7, pt.cex=1.5)
        par(opar)

}

}
