trim.me<-function(Y,trim=.2){
n.i<-length(Y)
Y1<-sort(Y)
bot<-floor(n.i*trim)+1
top<-n.i-bot+1
Y2<-Y1[bot:top]
Y2
}