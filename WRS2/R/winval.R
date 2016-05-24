winval<-function(x,tr=.2){
#
#  Winsorize the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
#  This function is used by several other functions that come with this book.
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
winval<-ifelse(x<=xbot,xbot,x)
winval<-ifelse(winval>=xtop,xtop,winval)
winval
}