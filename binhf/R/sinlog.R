"sinlog" <-
function(t){
t1<-t-0.4
y<-sin(5*pi*t1)+sin(6*pi*(t1-1/60))+sin(7*pi*(t1-1/35))+sin(8*pi*(t1-3/80))+sin(9*pi*(t1-2/45))+sin(10*pi*(t1-1/20))+sin(11*pi*(t1-3/55))
y<-y/7
y<-(0.7+y+cos((6/5)*pi*(t-.5)))/3
}

