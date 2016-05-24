"fun.disc.estimation" <-
function(x,nint){
# Create the categories:
categ<-pretty.su(x,nint)
len<-length(categ)
freq<- tabulate(cut(x,categ,include.lowest=TRUE))

# Caculate the mean and variance using categorised midpoint values.

midpoint<-(categ[-len]+categ[-1])/2

mean.e<-sum(midpoint*freq/sum(freq))

var.e<-sum((midpoint-mean.e)^2*freq/sum(freq))*(sum(freq)/(sum(freq)-1))

return(c(mean.e,var.e))
}

