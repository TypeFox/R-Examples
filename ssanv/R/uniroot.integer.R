"uniroot.integer" <-
function (f, interval, lower = min(interval), upper = max(interval), step.power=6,step.up=TRUE,pos.side=FALSE, print.steps=FALSE,
maxiter = 1000, ...) {
    ## iter counts how many times f is evaluated 
    iter<-0
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
        stop("lower < upper  is not fulfilled")
    if (lower==-Inf && step.up==TRUE) stop("lower cannot be -Inf when step.up=TRUE")
    if (upper==Inf && step.up==FALSE) stop("upper cannot be Inf when step.up=FALSE")
    if (step.up){
        f.old<-f(lower,...)
        iter<-iter+1
        sign<-1
        xold<-lower }
    else{
        f.old<-f(upper,...)
        iter<-iter+1
        sign<- -1 
        xold<-upper
        }
     if (print.steps){ print(paste("x=",xold," f(x)=",f.old)) }
     ever.switched<-FALSE
     tried.extreme<-FALSE
     while (step.power>-1){
        # Jun-22-2015: fix problems when f(i)=0 exactly for some i
        # break out when f.old=0, since xold will be the root 
        if (f.old==0) break()
        if (iter>=maxiter) stop("reached maxiter without a solution")
        xnew<- xold + sign*2^step.power
        if ((step.up & xnew< upper) || (!step.up & xnew> lower) ){  
		f.new<-f(xnew,...) 
            iter<-iter+1
		if (print.steps){ print(paste("x=",xnew," f(x)=",f.new)) }
		}
        else{ 
		#### Since stepped beyond extreme, move x back to xold
      	  xnew<- xold
		#### define f.new for the `if' statements following
            f.new<-f.old
		#### Decrease the step size if you step beyond the extreme
           step.power<-step.power-1 
		#### Only run the f(extreme) once, and test that you have opposite ends at both extremes
		if (tried.extreme==FALSE){
            	if (step.up){ f.extreme<-f(upper,...);iter<-iter+1; x.extreme<-upper }
			else{ f.extreme<-f(lower,...);iter<-iter+1; x.extreme<-lower }
               	tried.extreme<-TRUE 
                  xswitch<-x.extreme
                  f.switch<-f.extreme
			if (print.steps){ print(paste("x=",x.extreme," f(x)=",f.extreme)) }
                  # Jun-22-2015: fix problems when f(i)=0 exactly for some i
                  # break out when f=0, since then x will be the root
                  if (f.extreme==0){
                      # set xold to root and f.old=f(root) 
                      xold<-x.extreme
                      f.old<-f.extreme
                      break()
                  }
               	if (f.old*f.extreme>=0){stop("f() at extremes not of opposite sign")}
			}
		 }

        if (f.old*f.new <0){ 
                 sign<- sign*(-1)
                 ever.switched<-TRUE
		     xswitch<-xold
		     f.switch<-f.old
         }
        if (ever.switched){ 
			#### Only decrease the step size if you have already either switched directions
			#### (implying that the root is before the switch)
			#### Previously, we had decreased the step size if we had stepped beyond the extreme
                 step.power<-step.power-1 
                 if (step.power==-1){ break()} }


        xold<- xnew 
        f.old<-f.new

      }
    # Jun-22-2015: fix problems when f(i)=0 exactly for some i
    # for breaks out of while loop for f.old==0, make sure 
    # get the right root
    if (f.old==0){
       root<-xold
       f.root<-f.old
    } else if (f.new==0){
       root<-xnew
       f.root<-f.new
    } else if (f.switch==0){
       root<-xswitch
       f.root<-f.switch
    } else if (pos.side){  
        root<-ifelse(f.new>0,xnew,xswitch) 
        f.root<-ifelse(f.new>0,f.new,f.switch) 
   } else { 
        root<-ifelse(f.new<0,xnew,xswitch)
        f.root<-ifelse(f.new<0,f.new,f.switch) 
   }  
   # list(xold=xold,f.old=f.old,xswitch=xswitch,f.switch=f.switch,xnew=xnew,f.new=f.new,root=root)
   list(iter=iter,f.root=f.root,root=root)
}

