see.power<-function(alpha=NULL,sigma=NULL,n=NULL,effect=NULL,test="lower",xlim=c(-3,3),strict=FALSE){
    upper.titlel<-bquote(paste("Distribution assuming ",H[0],": ",mu >= 0))
    upper.titleu<-bquote(paste("Distribution assuming ",H[0],": ",mu <= 0))
    upper.titleb<-bquote(paste("Distribution assuming ",H[0],": ",mu," = 0"))
    lower.titlel<-bquote(paste("Distribution assuming ",H[A],": ",mu," < 0"))
    lower.titleu<-bquote(paste("Distribution assuming ",H[A],": ",mu," > 0"))
    lower.titleb<-bquote(paste("Distribution assuming ",H[A],": ",mu != 0))

    effect=abs(effect)
    layout(matrix(c(1,rep(2,2),rep(3,2)), 5, 1, byrow = TRUE))
    par(mar=c(4, 4, 2, 1))
    
    

if(test == "lower"){
    dev.hold()
    powerp<-power.z.test(alpha=alpha,sigma=sigma,effect=effect,power=NULL,n=n,test="one.tail")$power
    power1<-round(powerp,3)
    qalpha<-qnorm(alpha,mean=0,sd=sigma/sqrt(n))
    qpower<-qnorm(powerp,mean=-1*effect,sd=sigma/sqrt(n))
    plot(seq(0,1),seq(0,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(.5,0.8,bquote(paste(alpha," = ",.(alpha),",     1 - ",beta," = ", .(power1),",    ",italic(n)," = ",.(n),",    ",sigma," = ",.(sigma))),cex=2)
    text(.5,0.14,bquote(paste("effect size = ",.(effect),",     'tail' of test = '", .(test),"'")),cex=2)


    shade.norm(x=qalpha,sigma=sigma/sqrt(n),mu=0,tail="lower",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightyellow2",main=upper.titlel,xlim=xlim)
    abline(v=qalpha,col=1)
    legend("topright",pch=22,pt.bg="lightyellow2",pt.cex=2,legend=expression(alpha))
    text(qalpha,0.5*dnorm(0,mean=0,sd=sigma/sqrt(n)),"rejection region \u2190", adj = 1)
    shade.norm(x=qpower,sigma=sigma/sqrt(n),mu=-1*effect,tail="lower",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightcoral" ,main=lower.titlel,xlim=xlim)
    abline(v=qpower,col=1)
    legend("topright",pch=22,pt.bg="lightcoral",pt.cex=2,legend=expression(paste("1-",beta)))
    text(qpower,0.5*dnorm(effect,mean=effect,sd=sigma/sqrt(n)),"rejection region \u2190", adj = 1)
    dev.flush()
    }

if(test == "upper"){
    dev.hold()
    powerp<-power.z.test(alpha=alpha,sigma=sigma,effect=effect,power=NULL,n=n,test="one.tail")$power
    power1<-round(powerp,3)
    qalpha<-qnorm(1-alpha,mean=0,sd=sigma/sqrt(n))
    qpower<-qnorm(1-powerp,mean=effect,sd=sigma/sqrt(n))
    plot(seq(0,1),seq(0,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(.5,0.8,bquote(paste(alpha," = ",.(alpha),",     1 - ",beta," = ", .(power1),",    ",italic(n)," = ",.(n),",    ",sigma," = ",.(sigma))),cex=2)
    text(.5,0.14,bquote(paste("effect size = ",.(effect),",     'tail' of test = '", .(test),"'")),cex=2)

    shade.norm(x=qalpha,sigma=sigma/sqrt(n),mu=0,tail="upper",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightyellow2",main=upper.titleu,xlim=xlim)
    abline(v=qalpha,col=1)
    legend("topright",pch=22,pt.bg="lightyellow2",pt.cex=2,legend=expression(alpha))
    text(qalpha,0.5*dnorm(0,mean=0,sd=sigma/sqrt(n)),"\u2192 rejection region", adj = 0)
    
    shade.norm(x=qpower,sigma=sigma/sqrt(n),mu=effect,tail="upper",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightcoral" ,main=lower.titleu,xlim=xlim)
    abline(v=qpower,col=1)
    legend("topright",pch=22,pt.bg="lightcoral",pt.cex=2,legend=expression(paste("1-",beta)))
    text(qpower,0.5*dnorm(effect,mean=effect,sd=sigma/sqrt(n)),"\u2192 rejection region", adj = 0)
    dev.flush()
    }
if(test == "two"){
    dev.hold()
    powerp<-power.z.test(alpha=alpha,sigma=sigma,effect=effect,power=NULL,n=n,test="two.tail")$power
    if(strict == TRUE)powerp<-power.z.test(alpha=alpha,sigma=sigma,effect=effect,power=NULL,n=n,test="two.tail",strict=TRUE)$power
    power1<-round(powerp,3)
    qalpha<-qnorm(1-(alpha/2),mean=0,sd=sigma/sqrt(n))
    qpower<-qnorm(1-power.z.test(alpha=alpha,sigma=sigma,effect=effect,power=NULL,n=n,test="two.tail")$power,mean=effect,sd=sigma/sqrt(n))
    plot(seq(0,1),seq(0,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    text(.5,0.8,bquote(paste(alpha," = ",.(alpha),",     1 - ",beta," = ", .(power1),",    ",italic(n)," = ",.(n),",    ",sigma," = ",.(sigma))),cex=2)
    text(.5,0.14,bquote(paste("effect size = ",.(effect),",     'tail' of test = '", .(test),"'")),cex=2)

    shade.norm(x=qalpha,sigma=sigma/sqrt(n),mu=0,tail="two",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightyellow2",main=upper.titleb,xlim=xlim)
    abline(v=qalpha,col=1)
    abline(v=-qalpha,col=1)
    legend("topright",pch=22,pt.bg="lightyellow2",pt.cex=2,legend=expression(alpha))
    text(qalpha,0.5*dnorm(0,mean=0,sd=sigma/sqrt(n)),"\u2192 rejection region", adj = 0)
    text(-qalpha,0.5*dnorm(0,mean=0,sd=sigma/sqrt(n)),"rejection region \u2190", adj = 1)
    
    shade.norm(x=qpower,sigma=sigma/sqrt(n),mu=effect,tail="two",show.p=FALSE,show.d=FALSE,show.dist=TRUE,shade.col="lightcoral" ,main=lower.titleb,xlim=xlim)
    abline(v=qpower,col=1)
    abline(v=-qpower,col=1)
    legend("topright",pch=22,pt.bg="lightcoral",pt.cex=2,legend=expression(paste("1-",beta)))
    text(qpower,0.5*dnorm(effect,mean=effect,sd=sigma/sqrt(n)),"\u2192 rejection region", adj = 0)
    text(-qpower,0.5*dnorm(0,mean=0,sd=sigma/sqrt(n)),"rejection region \u2190", adj = 1)
    dev.flush()
    }
}


 
see.power.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    alpha <- 0.05
    sigma <- 1.5
    effect <-1
    n <- 5
    test <-tclVar("lower")
    strict <- tclVar("FALSE")
    
    assign("alpha", tclVar(alpha),envir= slider.env)
    assign("sigma", tclVar(sigma),envir= slider.env)
    assign("effect", tclVar(effect),envir= slider.env)
    assign("n", tclVar(n),envir= slider.env)
       
    xmin <- -3
    assign("xmin", tclVar(xmin),envir= slider.env)
    xmax <- 3
    assign("xmax", tclVar(xmax),envir= slider.env)
           
   norm.refresh <- function(...) {
        alpha <- as.numeric(evalq(tclvalue(alpha),envir= slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma),envir= slider.env))
        n <- as.numeric(evalq(tclvalue(n),envir= slider.env))
        effect <- as.numeric(evalq(tclvalue(effect),envir= slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin),envir= slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax),envir= slider.env))
        xixa <- c(xmin, xmax)
        test <- tclvalue(test)
        strict <- tclvalue(strict)
        see.power(alpha=alpha,sigma=sigma,n=n,effect=effect,xlim=xixa,test=test, strict=strict)
    }
                     
    m <- tktoplevel()
    tkwm.title(m, "Visualizing Power")
    tkpack(tklabel(m,text="      Visualizing Power      "))
    tkwm.geometry(m, "+0+0")
    
    tkpack(tklabel(m,text="            "))
    
    tkpack(fr <- tkframe(m), side = "top", anchor = "w")
    tkpack(fr1 <- tkframe(fr), side = "left",anchor = "w")
    tkpack(fr2 <- tkframe(fr), side = "right", anchor = "e")
    
     tkpack(tklabel(fr1, text = "  Test: "), side = "left")
        for ( i in c("lower", "upper", "two")){
            tmp <- tkradiobutton(fr1, text=i, variable=test, value=i)
            tkpack(tmp,anchor="w")
            }
    
    tkpack(tklabel(fr2, text = "               Strict: "), side = "left")
        for ( i in c("TRUE", "FALSE")){
            tmp <- tkradiobutton(fr2, text=i, variable=strict, value=i)
            tkpack(tmp,anchor="w")
            }
        
   
              
            
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03b1',font=c("Helvetica","10","italic"), width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.01, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc,envir= slider.env)
    evalq(tkconfigure(sc, variable = alpha),envir= slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03c3',font=c("Helvetica","10","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.5, 
        to = 3, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc,envir= slider.env)
    evalq(tkconfigure(sc, variable = sigma),envir= slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n",font=c("Helvetica","10","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc,envir= slider.env)
    evalq(tkconfigure(sc, variable = n),envir= slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Effect size", font=c("Helvetica","10"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 3, orient = "horiz", resolution = .1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc,envir= slider.env)
    evalq(tkconfigure(sc, variable = effect),envir= slider.env)
    
    
    
    tkpack(tklabel(m,text="            "))  
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e,envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmin),envir= slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e,envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmax),envir= slider.env)
    
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right") 
}



