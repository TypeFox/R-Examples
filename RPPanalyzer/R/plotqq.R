`plotqq` <-
function(x,fileName="qqplot_and_line.pdf") {

pdf(file=fileName)

      if( length(x)==4){
      mmnt.lines <- which (x[[4]][,"sample_type"]=="measurement")
     for (i in 1:ncol(x[[1]])){
     
          par(lwd=2)

          qqnorm(x[[1]][mmnt.lines ,i],col="green",pch=5
                ,main=c(paste("target: ",x[[3]]["target",i])
                        ,paste("Antibody: ",x[[3]]["AB_ID",i])
                        ,paste("incubation run: ",x[[3]]["incubation_run",i])))
               
          qqline(x[[1]][mmnt.lines ,i],col="red")
         }
       }
       
       else{
          mmnt.lines <- which (x[[3]][,"sample_type"]=="measurement")
     for (i in 1:ncol(x[[1]])){

          par(lwd=2)

          qqnorm(x[[1]][mmnt.lines ,i],col="green",pch=5
                ,main=c(paste("target: ",x[[2]]["target",i])
                        ,paste("Antibody: ",x[[2]]["AB_ID",i])
                        ,paste("incubation run: ",x[[2]]["incubation_run",i])))

          qqline(x[[1]][mmnt.lines ,i],col="red")
         }
       
       }
dev.off()

}

