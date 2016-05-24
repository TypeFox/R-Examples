overly2 <-
function(veg,Plot.no,y,sint) { 
# sint: sampling interval
# The next 2 lines set defaults if linewidths are null
#  defwidth<- is.null(l.widths)
#  if(defwidth == TRUE) l.widths<- seq(0.5,3.0,0.5)
# The same now with colors
#  defcol<- is.null(colors)
#  if(defcol == TRUE) colors<- rgb(0,0,0,seq(100,230,15),maxColorValue=255)
 vegtypes<- names(veg)
 veg<- veg^y
 veg<- veg[order(Plot.no),]                        # rearranging the series 
 ser<- as.integer(table(Plot.no))                  # number of releves per time series
 lev<- levels(as.factor(Plot.no))                  # original names of plots, used for legends
 nt<- length(ser)                                  # no. of time series
 dser<- rep(0,nt*nt)                               # distance matrix of time series
 dim(dser)<-c(nt,nt)
 mser<- rep(0,nt*nt*2)                             # labels of closest observations
 dim(mser)<-c(nt,nt,2)
 kids<-rep(0,nt)
 jj<-rep(0,2)
 drel<- as.matrix(dist(veg))                       # drel is the full distance matrix (all releves)
 nspec <- length(veg[1,])                          # number of species                   
# distance matrix of series, processing col by col (i), while row is (j)
# notation for the subset of the full distance matrix is drel[i1:i2,j1:j2]
 i1<- 0 ; i2<- 0
 for (i in 1:nt){
     i1<- i2+1  
     i2<- i2+ser[i]
     j1<- 0 ; j2<- 0
     for (j in 1:nt) {
        j1<- j2+1
        j2<- j2+ser[j]
        if (i != j) {
           m<- min(drel[i1:i2,j1:j2])
           jj<- which(drel[i1:i2,j1:j2] == m)
           dser[j,i]<- m                           # dser: distance matrix of time series
           mser[j,i,2]<- ceiling(jj[1]/ser[i])
           mser[j,i,1]<- jj[1]-((ser[i]*(mser[j,i,2]-1)))
        }
     }
 }
 out <- pco(dser,k=2)                             # This is PCOA
 tree<- spantree(dser)                                                  # minimum spanning tree
# shifting time series into proper position
 kids[2:nt]<- tree$kid    ; kids[1]<-kids[2]
 parents<- c(1:nt) 
 merged<- rep(0,nt) ; merged[1]<- 1
 shift<- rep(0,nt)
 shift[parents[2]]<- mser[kids[2],parents[2],1]-mser[kids[2],parents[2],2]+shift[kids[2]]
 merged[2]<- 1
 while(sum(merged) < nt) {
   for (i in 2:nt) {
      if (merged[kids[i]]==1 & merged[parents[i]]==0) {
      shift[parents[i]]<- mser[kids[i],parents[i],2]-mser[kids[i],parents[i],1]+shift[kids[i]]
      merged[parents[i]]<- 1
      }
   }
 }
 shift<- shift-min(shift)
 range<-max(shift+ser)
 align.series<- rep(0,nt*range*nspec)
 count.series<- rep(0,nt*range)
 dim(align.series)<- c(nt,range,nspec)
 dim(count.series)<- c(nt,range)
# back-transform veg
 veg<- veg^(1/y)
# linex1, linex2 and ltext are for the second plot
 linex1<- rep(0,nt)
 linex2<- rep(0,nt)
 ltext<- rep("x",nt)
 is<- 1
 for(i in 1:nt) {
    align.series[i,(shift[i]+1):(shift[i]+ser[i]),1:nspec]<- veg[is:(is+ser[i]-1),1:nspec]
    is<-is+ser[i]
    count.series[i,(shift[i]+1):(shift[i]+ser[i])]<- rep(1,ser[i])
    linex1[i]<- shift[i]+1
    linex2[i]<- shift[i]+ser[i]
    ltext[i]<- lev[i]
#    text((shift[i]+ser[i]),i,lev[i],pos=4,cex=0.6)
 }
# abline(v=0,lwd=1.0,col="gray")
 C<- colSums(count.series)
 S<- colSums(align.series)
 M<- S/C
 timescal<- seq(0,(range-1)*sint,sint)
# The output list
 M<- as.data.frame(M)
 colnames(M) <- vegtypes
 rownames(M) <- as.character(seq(1,range,1))
overly<- list(plot.labels=levels(Plot.no),n.tsteps=range,tsteps=timescal,tser.data=M,ord.scores=out$points,d.mat=dser,tree=tree,vegraw=veg,linex1=linex1,linex2=linex2,ltext=ltext,sint=sint,vegtypes=vegtypes)
 }
