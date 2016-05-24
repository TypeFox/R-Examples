plotbary.slide<-function(tt)
{
d<-dim(tt$center)[1]
coordi<-1
plotbary(tt,paletti=seq(1:1000),coordi=coordi)
loc<-locator(1)
while (loc$y>=0){
    if (coordi==d) coordi<-1 else coordi<-coordi+1         
    plotbary(tt,paletti=seq(1:1000),coordi=coordi)          
    loc<-locator(1)
}

}

