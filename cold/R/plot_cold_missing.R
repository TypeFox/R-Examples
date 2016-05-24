
setMethod(f="plot",signature=c(x="cold",y="missing"), 
definition=function(x,which=c(1:3),xlab=NULL,ylab=NULL,main=NULL,
subSET,ident=FALSE,caption=c("Individual mean profiles"),cex.caption=1)

{


if (!is.numeric(which) || any(which < 1) || any(which > 3)) 
       stop("'which' must be in 1:3")

	if(missing(which)) 
#        stop("choose 'which' in 1:3")
         {which=1}


show <- rep(FALSE, 3)
    show[which] <- TRUE
 


 x2<-x@Fitted 

 depend<- x@call$dependence




##################################################
  factor<-x@call$aggregate
  subset1<-x@call$subSET

  if (length(factor)==0 && length(subset1)==0)
	{ trace.label=" "}

  else if(length(factor)==0 && length(subset1)!=0)
  	{trace.label=deparse(substitute(subset1))}

  else if (length(factor)!=0)
  	{trace.label=deparse(substitute(factor))}


#################################################

   time<-x@Time 
   x3<-unique(x@Time)
   n.time<-length(x3)

   x4<-unique(x@model.matrix)

   n.id<-x@n.cases
   data<-x@subset.data


  data$id<-x@data.id  

  data$y<-as.vector(t(x@y.matrix))  



  id.new<-unique(data$id)

  f<-x@f.value
  level.f<-unique(f)
  nlevel.f<-length(unique(f))

  levg<-as.vector(level.f)
  colg<-seq(1:nlevel.f)

  n<-dim(x4)[1]

  ncurves<-nlevel.f

 	aa<-matrix(x@Fitted,n.id,n.time,byrow=TRUE)
 	aanew<-numeric(n.id)
	x3new<-numeric(n.id)

	pos.id<-seq(1:n.id)
	id.numb<-unique(data$id)
	afit<-data.frame(id.numb,pos.id,aa)

        xlab2<-"Time"
 
if (any(show[1:2])){
         if(length(ylab)==0) ylab="Number"
         else ylab<-ylab
         if(length(main)==0) main=" "
         else main<-main
                   }

if (any(show[1:2])){
         if(length(xlab)==0) xlab="Time"
         else xlab<-xlab
         if(length(main)==0) main=" "
         else main<-main
                   }
 


	if(show[1]) {


##################################
##      Parametric fit          ##
##################################
	         


 		coef<-x@coefficients
 		npar <-length(coef)
 		omega<-as.double(coef[npar])
 
 		n.col<-ncol(x4)
 		coef1<-as.double(x@coefficients[1:n.col])
 		eta<-x4%*%coef1
         

 		n.time<-length(x3)
                l<-exp(eta)

		ylim<-range(l,na.rm=TRUE)
		ylim<-extendrange(l,f=0.08)
 
 		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab,ylab=ylab,main=" ",ylim=ylim)
		
                mtext(main,side=3,0.25, cex=1)

		n1<-1
		n2<-n.time
		
		p<-round(ylim[2])
		p1<-ylim[2]-0.08
                


		for(j1 in 1:ncurves)
				{
                n2<-j1*(n.time)
		lines(x3,exp(eta[n1:n2]),col=j1,lty=1)
          	n1<-n2+1
		                }



		if(ncurves>1){legend(x3[2],p1,levg,
 		lty=c(1,1),col=colg,bty="n",cex=0.75)}
                 
		if(ncurves==1){legend(x3[1],p,paste(trace.label),bty="n",cex=0.75)}

			
	
                     }



####################################################

##		Individual mean profiles          ##

####################################################

  
##	Individual mean profiles

    	if(show[2])  {
 
	   if (depend=="ind"||depend=="AR1") 
       		stop("dependence must be indR or AR1R")
		


		l1<-aa
		ylim<-range(l1,na.rm=TRUE)
		ylim<-extendrange(l1,f=0.08)

		

		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab,ylab=ylab,main=" ",ylim=ylim)
 		mtext(main,side=3,1, cex=0.8)
		mtext(caption[1],3,0.25,cex=0.8)
                

            if(missing(subSET)){
 		for(i in 1:n.id) {
 		lines(x3,aa[i,],col=1,lty=i)
                if(i<=n.time){aanew[i]<-aa[i,i]
                x3new[i]<-x3[i]}
		else{aanew[i]<-aa[i,n.time]
                x3new[i]<-x3[n.time]}
                     } 

             	aanew<-matrix(aanew,n.id,1)
             	aanew1<-data.frame(id.new,x3new,aanew)
                if(ident==TRUE)
		{
		text(aanew1[,2],aanew1[,3],labels=aanew1[,1],cex=0.8)
		}

                           }


	   if(!missing(subSET)){
		

		id1 <- eval(substitute(subSET), data)
	        data<-subset(data, id1)
	
                n.subj<-length(unique(data$id))
             
                id.new1<-unique(data$id)
                aanew2<-numeric(n.subj)
		x3new1<-numeric(n.subj)

                for(j in 1:n.subj) {

                id.numbj<-id.new1[j]
                i<-afit[afit$id.numb==id.numbj,2]
         
                
		lines(x3,aa[i,],col=1,lty=i)
		if(i<=n.time){aanew2[j]<-aa[i,i]
                x3new1[j]<-x3[i]}
		else{aanew2[j]<-aa[i,n.time]
                x3new1[j]<-x3[n.time]}
                                   }     

            	aanew2<-matrix(aanew2,n.subj,1)
             	aanew3<-data.frame(id.new1,x3new1,aanew2)
		if(ident==TRUE){ 
		text(aanew3[,2],aanew3[,3],labels=aanew3[,1],cex=0.8)} 

                               }
         	
                 }

	if(show[3])  {

		if (depend=="ind"||depend=="AR1") 
       		stop("dependence must be indR or AR1R")

		if(missing(subSET))
		stop("number of id must be given")


		l1<-data$y
		ylim<-range(l1,na.rm=TRUE)
		ylim<-extendrange(l1,f=0.08)


		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab,ylab=ylab,main=" ",ylim=ylim)

 		mtext(main,side=3,2, cex=0.8)
		mtext(caption[1],3,0.25,cex=0.8)



		id1 <- eval(substitute(subSET), data)
	        data<-subset(data, id1)
	
                n.subj<-length(unique(data$id))
             
                id.new1<-unique(data$id)
                aanew2<-numeric(n.subj)
		x3new1<-numeric(n.subj)
                ynew<-data$y

                for(j in 1:n.subj) {

                id.numbj<-id.new1[j]
                i<-afit[afit$id.numb==id.numbj,2]
                p2<-round(ylim[2])

		p3<-p2-1

                points(x3,ynew,pch=19,col="grey50",cex=1.3)

		points(x3,aa[i,],pch=19,col="grey70")
		if(i<=n.time){aanew2[j]<-aa[i,i]
                x3new1[j]<-x3[i]}
		else{aanew2[j]<-aa[i,n.time]
                x3new1[j]<-x3[n.time]}
                                   }     

            	aanew2<-matrix(aanew2,n.subj,1)
            	aanew3<-data.frame(id.new1,x3new1,aanew2)
                
		legend(x3[1],p2,c("Observed","Fitted"),pch=c(19,19),col=c("grey50","grey70"),cex=0.8,bty="n")


                

			}

}

 )



