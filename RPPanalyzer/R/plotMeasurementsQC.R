`plotMeasurementsQC` <-
function(x,file="QC_plots.pdf",arrays2rm=c("protein")) {

         x <- remove.arrays(x,param="target",arrays2rm=arrays2rm)
         
         inc.run <- levels(as.factor(x[[3]]["incubation_run",]))

         mmnt.lines <- which (x[[4]][,"sample_type"]=="measurement")
     pdf(file=file)
 for (j in seq(along=inc.run)){

         ## define columns of a single incubation run
          run.cols <- which((x[[3]]["incubation_run",])==inc.run[j])
          ## select data of a single incubation run
            temp.dat <- x[[1]][,run.cols]
            temp.arrays <- x[[3]] [,run.cols]
            blk.col <- grep("lank$",temp.arrays["target",])

     for (i in 1:ncol(temp.dat)){

          par(lwd=2)

          plot(temp.dat[mmnt.lines,i],temp.dat[mmnt.lines,blk.col],col="darkgreen",pch=5
               , main=c(paste("target: ",temp.arrays["target",i],"  Ab: ",temp.arrays["AB_ID",i])
                        ,paste("incubation run: ",inc.run[j]))
               ,xlab="target signal",ylab="blank signal"
               ,ylim=c(0,max(c(temp.dat[mmnt.lines,i],temp.dat[mmnt.lines,blk.col])))
               ,xlim=c(0,max(c(temp.dat[mmnt.lines,i],temp.dat[mmnt.lines,blk.col]))))

          abline(0,1,col="red")
          }
  }
   dev.off()

}

