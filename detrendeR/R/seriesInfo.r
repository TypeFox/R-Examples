seriesInfo = function (rwl, print=T)
{
    yr.range = function(x) {
        yr.vec = as.numeric(names(x))
        mask = !is.na(x)
        range(yr.vec[mask])
    }

info.fun = function (x) {
out=c(rep(NA,3))
        out[1] = yr.range(x)[1]
      out[2] = yr.range(x)[2]
      out[3] = out[2]-out[1]+1

return(out)
}
out=t(apply(rwl,2, info.fun))

        colnames(out)<-c("First","Last","Span")



if (print){
cat(rep("=", 34),"\n",sep="")
WriteMatrix(out, na="", sep="|",ID=T, ID.name="Seq",col.width=6, row.name="Series   ")
cat(rep("=", 34),"\n",sep="")
} else {
out<-cbind(rownames(out), out)
out<-cbind(1:nrow(out), out)
out <- rbind(c("ID", "Series", "First", "Last", "Span"), out)
format(out[,2], justify="left", width=10)->out[,2]
for(i in 1:5) format(out[,i], justify="right")->out[,i]


return(out)
}
}
