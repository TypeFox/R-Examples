reg.fun.default <-
function(Vec,y,type_data=c("discrete","continuous"),ker=c("bino","triang","dirDU","BE","GA","LN","RIG"),
			h,
			x=NULL,
			a0=0,
			a1=1,
			a=1,
			c=2,
			...)
    	{
if (missing(type_data))  stop("argument 'type_data' is omitted")
    if ((type_data=="discrete") & (ker=="GA"||ker=="LN"||ker=="BE" ||ker=="RIG")) 
       stop(" Not appropriate kernel for type_data")
   if ((type_data=="continuous") & (ker=="bino"||ker=="triang"||ker=="dirDU")) 
      stop(" Not appropriate kernel for 'type_data'")
   if ((type_data=="discrete") & missing(ker)) ker<-"bino"
   if ((type_data=="continuous") & missing(ker)) ker<-"GA"

  if(is.null(x)){x=Vec}
    n <- length(x)
    m=matrix(0,n,length(Vec))
    m2=matrix(0,n,length(Vec))

    for(i in 1:n){
     for(j in 1:length(Vec)){
         m[i,]= kef(x[i],Vec,h,type_data,ker,a0,a1,a,c)
         m2[i,]= m[i,]*y
	        }
    }
    res<-apply(m,1,sum)
    res2<-apply(m2,1,sum)
    result<-res2/res
    moyY<-mean(y)
    R<-sum(((result-moyY)^2))/sum(((y-moyY)^2)) # Coefficient de détermination R^2
		
	structure(list(data=Vec,y=y,n=length(Vec),kernel=ker,h=h, 
                  eval.points=x, m_n = result,Coef_det=R),class="reg.fun")		

}
