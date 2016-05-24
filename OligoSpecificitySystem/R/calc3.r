"calc3"<-function(){
a1<-c(1,0,0,1,1,0,1)
a2<-c(0,1,0,1,0,1,1)
a3<-c(0,0,1,0,1,1,1)
a4<-c(1,1,0,1,1,1,1)
a5<-c(1,0,1,1,1,1,1)
a6<-c(0,1,1,1,1,1,1)
a7<-c(1,1,1,1,1,1,1)
gauche<-rbind(a1,a2,a3,a4,a5,a6,a7)
droite<-c(length(txt1),length(txt2),length(txt3),length(levels(factor(c(txt1,txt2)))),length(levels(factor(c(txt1,txt3)))),length(levels(factor(c(txt2,txt3)))),length(levels(factor(c(txt1,txt2,txt3)))))
#droite<-c(bd1,bd2,bd3,4,5,4,1)
zz<-solve(gauche,droite)
zz1<-matrix(nr=1,nc=7)
colnames(zz1)<-c("a","b","c","d","e","f","g")
zz1[1,]<-zz
zz1<<-zz1 }