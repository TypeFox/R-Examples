kullnagar.boot<-function(data, i, ...)
{
         data$Observed<-data$Observed[i]

         kullnagar.stat(data, ...)

}
