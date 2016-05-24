see.adddel<-function(){

action<-tclVar("add") 
pts <- matrix(c(12, 56,
                20, 94,
                33, 98,
                32, 120,
                61, 180,
                75, 160,
                98, 223), ncol=2, byrow=TRUE)
x<-pts[,1]
y<-pts[,2]
dev.new(height=5,width=5*1.4)
par(mar=c(5,4,1,9),xpd=FALSE,bg="gray90")
plot(x,y,xlab=expression(italic(x)),ylab=expression(italic(y)),type="n")
pu<-par("usr")
rect(pu[1],pu[3],pu[2],pu[4],col="white")
points(x,y,pch=21,cex=1.2,bg="skyblue2")
l<-lm(y~x);co<-coef(l);c1<-summary(l)$coefficients;p<-c1[2,4]
abline(co[1],co[2])  
par(xpd=TRUE)
legend(pu[2]+.5,pu[4],legend=c(paste(" Y int = ",round(co[1],4)),paste(" Slope = ",round(co[2],4)),"",paste(" r = ",round(cor(y,x),4)),paste(" r\u00b2 = ",round(cor(y,x)^2,4)),"",paste(" P-value = ",round(p,4))),box.col="gray90",cex=.8,bg="white")
box()

add.points<-function(){
loc<-locator(1)
x<-c(x,loc$x)
y<-c(y,loc$y)
cbind(x,y) 
} 

delete.points<-function(){
ans<-identify(x,y, n=1, plot=FALSE)
x<-x[-ans]
y<-y[-ans]
cbind(x,y)
}
 
repeat{
    refresh<-function(){
        
        action<-tclvalue(action)
        if(action=="add")m<-add.points()
        if(action=="delete")m<-delete.points()
        tkdestroy(tt)
        x<<-m[,1];y<<-m[,2]
        dev.hold()
        par(xpd=FALSE)
        plot(x,y,xlab=expression(italic(x)),ylab=expression(italic(y)),type="n")
        pu<-par("usr")
        rect(pu[1],pu[3],pu[2],pu[4],col="white")
        points(x,y,pch=21,cex=1.2,bg="skyblue2")
        l<-lm(y~x);co<-coef(l);c1<-summary(l)$coefficients;p<-c1[2,4]
        abline(co[1],co[2])  
        par(xpd=TRUE)
        legend(pu[2]+.5,pu[4],legend=c(paste(" Y int = ",round(co[1],4)),paste(" Slope = ",round(co[2],4)),"",paste(" r = ",round(cor(y,x),4)),paste(" r\u00b2 = ",round(cor(y,x)^2,4)),"",paste(" P-value = ",round(p,4))),box.col="gray90",cex=.8,bg="white")
        box()
        dev.flush()
        }

tclServiceMode(TRUE)   
tt <- tktoplevel()
tkwm.geometry(tt, "+50+4")
tkwm.title(tt, "Demonstration of least squares regression -- Add/delete points")
tkpack(tklabel(tt,text="Adding/deleting points in\nsimple linear regression"))
tkpack(tklabel(tt,text=""))
tkpack(tklabel(tt, text = "  Action: "), side = "top")
        for ( i in c("add", "delete")){                           
            tmp <- tkradiobutton(tt, text=i, variable=action, value=i)
            tkpack(tmp, anchor = "w")}
tkpack(tkbutton(tt, text = "Exit", command = function()tkdestroy(tt))) 
refresh()
}
}

                