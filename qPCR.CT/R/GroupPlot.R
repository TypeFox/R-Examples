GroupPlot <-
function(filename, 
                  con_group_name,tr_group_name, 
                  con_gene_name, repeatnum, ... ){
                    
rawdata<-read.csv(filename)
unknown<-   subset(rawdata, Type != "NTC")
suppressWarnings(aggdata <-aggregate.data.frame(unknown, by=list(unknown$Fluor),FUN=mean, na.rm=TRUE))
#print(aggdata[1])
groups<-paste(aggdata[,1])
lengr <- length(groups)
temp<- data.frame()
i<- 1
for(i in 1: lengr)
{
   if(groups[i] != con_gene_name ){
     
          tr_gene_name <- groups[i] # assign group names
          temp <- Calc.Expr(filename=filename, con_group_name, tr_group_name, 
          con_gene_name = con_gene_name,tr_gene_name = tr_gene_name, repeatnum = repeatnum)
          plot.CT(temp)
   }
}
}

