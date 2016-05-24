kern.fun.default <-
function(x,t,h,type_data=c("discrete","continuous"),
				   ker=c("bino","triang","dirDU","BE","GA","LN","RIG"),a0=0,a1=1,a=1,c=2,...)
              {

    if (missing(type_data))  stop("argument 'type_data' is omitted")
    if ((type_data=="discrete") & (ker=="GA"||ker=="LN"||ker=="BE" ||ker=="RIG")) 
       stop(" Not appropriate kernel for type_data")
   if ((type_data=="continuous") & (ker=="bino"||ker=="triang"||ker=="dirDU")) 
      stop(" Not appropriate kernel for 'type_data'")
   if ((type_data=="discrete") & missing(ker)) ker<-"bino"
   if ((type_data=="continuous") & missing(ker)) ker<-"GA"
     kx <- kef(x,t,h,type_data,ker,a0,a1,a,c)
     structure(list(kernel = ker,x=x,t=t,kx=kx),class="kern.fun")
}
