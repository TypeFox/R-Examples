investigate <- function() {

  # We investigate here the examples and the exercices from the article cited in the references section above

    above.env <- environment()
    current.env <- new.env()

    local({
    graphics.off()
    
    exists.investigate <- exists("tt.investigate", envir = above.env)
    if (exists.investigate) {eval(parse(text="tkdestroy(tt.investigate)"), envir = above.env);eval(parse(text="rm(tt.investigate,pos=1)"), envir = above.env)}
    
    for (i in 1:9) {
        
        x1 <- paste("tt",format(i),".1",sep="")
        x2 <- paste("tt",format(i),".2",sep="")
        exists.x1 <- exists(x1, envir = above.env)
        exists.x2 <- exists(x2, envir = above.env)
        
        if (exists.x1) {eval(parse(text=paste("tkdestroy(",x1,")",sep="")), envir = above.env);eval(parse(text=paste("rm(",as.character(x1),",pos=1)",sep="")), envir = above.env)}
        if (exists.x2) {eval(parse(text=paste("tkdestroy(",x2,")",sep="")), envir = above.env);eval(parse(text=paste("rm(",as.character(x2),",pos=1)",sep="")), envir = above.env)}
        
    }
    
    
    plotright <- function(...) {
        
        val <- as.numeric(tkcurselection(tl))+1
        if (length(val) == 0) val <- 0
        
        if (val == 0) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let Z be a uniform U[0,1] random variable."))
            text(0.5,0.8,expression("Define " ~ X[n]=="1"["[m.2^(-k);(m+1).2^(-k))"](Z)))
            text(0.5,0.7,expression("where" ~ n==2^k+m ~ "for " ~ k<=1 ~ "and with " ~ 0<=m ~ ""<2^k ~ "."))
            text(0.5,0.6,expression("Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~0 ~ "? Does" ~ X[n] ~ " " ~ frac(P,"     ")~">" ~0 ~ "?"))
        }
        
        if (val == 1) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let Z be a uniform U[0,1] random variable."))
            text(0.5,0.8,expression("Define " ~ X[n]=="1"["[m.2^(-k);(m+1).2^(-k))"](Z)))
            text(0.5,0.7,expression("where" ~ n==2^k+m ~ "for " ~ k<=1 ~ "and with " ~ 0<=m ~ ""<2^k ~ "."))
            text(0.5,0.6,expression("Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~0 ~ "? Does" ~ X[n] ~ " " ~ frac(P,"     ")~">" ~0 ~ "?"))
            
        }
        
        
        if (val == 2) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ X[1] ~ ",... , " ~ X[n] ~ "be i.i.d. N(0,1) random variables and " ~ X==X[1] ~ "."))
            text(0.5,0.8,expression("Does" ~ X[n] ~ " " ~ frac(L,"     ")~">" ~X ~ "? Does" ~ X[n] ~ " " ~ frac(P,"     ")~">" ~X ~ "?"))
        }
        
        if (val ==3) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ X[1] ~ ",... , " ~ X[n] ~ "be independent random variables"))
            text(0.5,0.8,expression("such that P[" ~ X[n]==sqrt(n) ~ "]=1/n and P[" ~ X[n]==0 ~ "]=1-1/n."))
            text(0.5,0.7,expression("Does" ~ X[n] ~ " " ~ frac(2,"     ")~">" ~0 ~ "? Does" ~ X[n] ~ " " ~ frac(P,"     ")~">" ~0 ~ "?"))
        }
        
        
        if (val == 4) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let Z be U[0,1] and let" ~ X[n]==2^n ~ "1"["[0,1/n)"](Z)))
            text(0.5,0.7,expression("Does" ~ X[n] ~ " " ~ frac("r","     ")~">" ~0 ~ "? Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~0 ~ "?"))
        }
        
        
        if (val ==5) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ Y[1] ~ ",... , " ~ Y[n] ~ "be independent random variables with mean 0"))
            text(0.5,0.8,expression("and variance 1. Define" ~ X[1]==X[2] ~ "=1 and"))
            text(0.5,0.68,expression(X[n]==frac(sum(Y[i], i==1, n),"(2n log log n)"^{1/2}) ~ "," ~ n>=3))
            text(0.5,0.5,expression("Does" ~ X[n] ~ " " ~ frac(2,"     ")~">" ~0 ~ "? Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~0 ~ "?"))
            
        }
        
        
        if (val == 6) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let Z be Unif [0,1]."))
            text(0.5,0.8,expression("Let " ~ Y[1] ~ ",... , " ~ Y[n] ~ "be i.i.d. Unif{0,1,...,9} and let " ~ X[n]==sum(frac(Y[i],10^i), i==1, n) ~ "."))
            text(0.5,0.7,expression("It can be proved that " ~ X[n] ~ " " ~ frac("a.s.","     ") ~ ">" ~ X==sum(frac(Y[i],10^i), i==1, infinity)))
            text(0.5,0.6,expression("which follows a Unif [0,1] distribution."))
            text(0.5,0.5,expression("Does" ~ X[n] ~ " " ~ frac(L,"     ")~">" ~Z ~ "? Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~Z ~ "?"))
            
            
        }
        
        if (val == 7) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ Y[1] ~ ",... , " ~ Y[n] ~ "be i.i.d. N(0,1) random variables and let " ~ X[n]==frac(1,n) ~ sum(Y[i], i==1, n) ~ "."))
            text(0.5,0.8,expression("Does" ~ X[n] ~ " " ~ frac(P,"     ")~">" ~0 ~ "?"))
            text(0.5,0.7,expression("Does" ~ X[n] ~ " " ~ frac("a.s.","     ")~">" ~0 ~ "?"))
            
        }
        
        
        
        if (val == 8) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ X[1] ~ ",... , " ~ X[n] ~ "be independent random variables"))
            text(0.5,0.8,expression("such that P[" ~ X[n]==n^0.4 ~ "]=1/n and P[" ~ X[n]==0 ~ "]=1-1/n."))
            text(0.5,0.7,expression("Does" ~ X[n] ~ " " ~ frac(r,"     ")~">" ~0 ~ "?, r=1, 2, 3."))
            
            
        }
        
        if (val == 9) {
            
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            box()
            text(0.5,0.9,expression("Let " ~ Z[1] ~ ",... , " ~ Z[n] ~ "be " ~ chi[1]^{2} ~ " independent random variables."))
            text(0.5,0.7,expression("Let " ~ X[n]==frac(1,sqrt(n))*bgroup("[",frac(sum(Z[i], i==1, n)-n,sqrt(2)),"]") ~ "."))
            text(0.5,0.5,expression("Does" ~ X[n] ~ " " ~ frac(L,"     ")~">" ~"N(0,1)" ~ "?"))
            
            
        }
        
        
    }
    
    
    listboxfunc <- function(...) {
        
        tkrreplot(imgright)
        
    }
    
    assign("tt.investigate", tktoplevel(), above.env)
    eval(parse(text="tkgrid(tklabel(tt.investigate,text=\"What do you want to investigate?\"))"), envir = current.env)
    exercises <- c("Exercise 1","Exercise 2","Exercise 3","Exercise 4","Exercise 5","Exercise 6","Example 1","Example 2","Example 3","Example 4")
    tl <- eval(parse(text="tklistbox(tt.investigate,height=9,selectmode=\"single\",background=\"white\")"), envir = current.env)
    if (.Platform$OS.type == "unix") {imgright <- eval(parse(text="tkrplot(tt.investigate,plotright,hscale=1.2,vscale=1)"), envir = current.env)}
    if (.Platform$OS.type == "windows") {imgright <- eval(parse(text="tkrplot(tt.investigate,plotright,hscale=2,vscale=2)"), envir = current.env)}
    tkgrid(tl,imgright)
    tkbind(tl, "<ButtonRelease-1>",listboxfunc)
    assign("val", 0, above.env)
    
    for (i in 1:9)
        {
            tkinsert(tl,"end",exercises[i])
        }
    tkselection.set(tl,0)  # Default value.  Indexing starts at zero.
    
    
    OnOK <- function(...)
        {
            
            graphics.off()
            
            exerciseChoice <- exercises[as.numeric(tkcurselection(tl))+1]
            
            for (i in 1:9) {
                
                x1 <- paste("tt",format(i),".1",sep="")
                x2 <- paste("tt",format(i),".2",sep="")
                exists.x1 <- exists(x1, envir = above.env)
                exists.x2 <- exists(x2, envir = above.env)
                
                if (exists.x1) {eval(parse(text=paste("tkdestroy(",x1,")",sep="")), envir = above.env);eval(parse(text=paste("rm(",as.character(x1),",pos=1)",sep="")), envir = above.env)}
                if (exists.x2) {eval(parse(text=paste("tkdestroy(",x2,")",sep="")), envir = above.env);eval(parse(text=paste("rm(",as.character(x2),",pos=1)",sep="")), envir = above.env)}
                
            }
            
            if (exerciseChoice=="Exercise 1") {
                
########### Exercise 1 ###########

                pnotasgen <- function(n){
                    Z<-runif(1)
                    k<-floor(log2(1:n))
                    m<-1:n-2^k
                    res<-(m*2^(-k)<= Z & Z<(m+1)*2^(-k))
                    return(as.integer(res))
                }

                assign("tt1.1", check.convergence(nmax=2000,M=500,genXn=pnotasgen,mode="as"), above.env)
                
            }
            
            if (exerciseChoice=="Exercise 2") {
                
########### Exercise 2 ###########
                
                lnotpgen <- function(n){x<-rnorm(n);x-x[1]}

                assign("tt2.2", check.convergence(nmax=2000,M=500,genXn=lnotpgen,mode="p"), above.env)
                
            }
            
            if (exerciseChoice=="Exercise 3") {
                
########### Exercise 3 ###########
                
                pnotrgen <- function(n){rbinom(n,1,1/(1:n))*sqrt(1:n)}
                
                                        #       if (.Platform$OS.type=="unix") X11()
                                        #       if (.Platform$OS.type=="windows") windows()
                                        #       if (.Platform$OS.type=="mac") quartz()
                dev.new()
                
                check.convergence(nmax=1000,M=10000,genXn=pnotrgen,mode="r",r=2,ylim=c(0,5))
                legend("topleft",legend=expression(hat(e)[n~bold(',')~'2']),lty=1)
                assign("tt3.1", check.convergence(nmax=2000,M=500,genXn=pnotrgen,mode="p"), above.env)
                
                
            }
            
            if (exerciseChoice=="Exercise 4") {
                
########### Exercise 4 ###########
                
                asnotrgen <- function(n){x<-2^(1:n);res<-(runif(1)<1/(1:n))*x;res[is.infinite(x)]<-0;return(res)}

                                        #       if (.Platform$OS.type=="unix") X11()
                                        #       if (.Platform$OS.type=="windows") windows()
                                        #       if (.Platform$OS.type=="mac") quartz()
                dev.new()

                check.convergence(nmax=10,M=500,genXn=asnotrgen,mode="r",r=2)
                legend("topleft",legend=expression(hat(e)[n~bold(',')~'2']),lty=1)
                assign("tt4.1", check.convergence(nmax=2000,M=500,genXn=asnotrgen,mode="as"), above.env)

            }

            if (exerciseChoice=="Exercise 5") {

########### Exercise 5 ###########

                rnotasgen <- function(n){if (n == 1) res <- 1 else if (n == 2) res <- c(1,1) else res<-c(1,1,cumsum(rnorm(n))[-(1:2)]/sqrt(2*(3:n)*log(log(3:n))))}

                                        #       if (.Platform$OS.type=="unix") X11()
                                        #       if (.Platform$OS.type=="windows") windows()
                                        #       if (.Platform$OS.type=="mac") quartz()
                dev.new()

                check.convergence(nmax=2000,M=500,genXn=rnotasgen,mode="r",r=2,col="red")
                points(3:1000,1/(2*log(log(3:1000))),type="l",col="blue")
                legend("topright",legend=c(expression(hat(e)[n~bold(',')~'2']),expression(e[n~bold(',')~'2'])),col=c("red","blue"),lty=1)
                assign("tt5.1", check.convergence(nmax=2000,M=500,genXn=rnotasgen,mode="as"), above.env)

            }

            if (exerciseChoice=="Exercise 6") {

########### Exercise 6 ###########

                gen6.1 <- function(n){res<-cumsum(floor(10*runif(n))/(10^(1:n)))}
                gen6.2 <- function(n){res<-cumsum(floor(10*runif(n))/(10^(1:n)))-runif(1)}

                                        #       if (.Platform$OS.type=="unix") X11()
                                        #       if (.Platform$OS.type=="windows") windows()
                                        #       if (.Platform$OS.type=="mac") quartz()
                dev.new()

                assign("tt6.1", check.convergence(nmax=20,M=5000,genXn=gen6.1,mode="L",density=FALSE,densfunc=dunif,probfunc=punif,tinf=-0.1,tsup=1.1), above.env)
                assign("tt6.2", check.convergence(nmax=2000,M=500,genXn=gen6.2,mode="as"), above.env)

            }


            if (exerciseChoice=="Example 1") {

########### Example 1 ###########

                moyrand <- function(n,...){cumsum(rnorm(n,...))/(1:n)}

                assign("tt7.1", check.convergence(nmax=2000,M=500,genXn=moyrand,mode="p"), above.env)

            }


            if (exerciseChoice=="Example 2") {

########### Example 2 ###########

                myrbinom <- function(n,alpha){rbinom(n,1,1/(1:n))*((1:n)**alpha)}

                                        #       if (.Platform$OS.type=="unix") X11()
                                        #       if (.Platform$OS.type=="windows") windows()
                                        #       if (.Platform$OS.type=="mac") quartz()
                dev.new()

                check.convergence(nmax=2000,M=500,genXn=myrbinom,argsXn=list(alpha=0.5),mode="r",r=3,plotfunc=plot,col="green",ylim=c(0,300))
                check.convergence(nmax=2000,M=500,genXn=myrbinom,argsXn=list(alpha=0.5),mode="r",r=2,plotfunc=points,col="blue",ylim=c(0,300))
                check.convergence(nmax=2000,M=500,genXn=myrbinom,argsXn=list(alpha=0.5),mode="r",r=1,plotfunc=points,col="red",ylim=c(0,300))
                text(100,-5,'r=1',xpd=TRUE,col='red')
                text(200,5,'r=2',xpd=TRUE,col='blue')
                text(300,50,'r=3',xpd=TRUE,col='green')
                legend("topleft",legend=c(expression(hat(e)[n~bold(',')~'1']),expression(hat(e)[n~bold(',')~'2']),expression(hat(e)[n~bold(',')~'3'])),col=c("red","blue","green"),lty=1)

            }

            if (exerciseChoice=="Example 3") {

########### Example 3 ###########

                
                                        #      if (.Platform$OS.type=="unix") {X11(width=2,height=2);plot.new();title("Please wait ...");X11()}
                                        #      if (.Platform$OS.type=="windows") {windows(width=2,height=2);plot.new();title("Please wait ...");windows()}
                                        #      if (.Platform$OS.type=="mac") {quartz(width=2,height=2);plot.new();title("Please wait ...");quartz()}
                dev.new(width=2,height=2);plot.new();title("Please wait ...");dev.new()

                if (exists("bringToTop", envir = above.env)) bringToTop(dev.prev())
                rand1 <- function(n){(cumsum(rchisq(n,df=1))-(1:n))/sqrt(2*(1:n))}
                assign("tt8.1", check.convergence(nmax=200,M=5000,genXn=rand1,mode="L",density=FALSE,densfunc=dnorm,probfunc=pnorm,tinf=-4,tsup=4), above.env)
                dev.off() 
            }
        }

    OK.but <- eval(parse(text="tkbutton(tt.investigate,text=\"   OK   \",command=OnOK)"), envir = current.env)
    tkgrid(OK.but)

},envir=current.env)
    
}
