setMethod(f="plot",signature=c(x="bild",y="missing"), 
definition=function(x,which=c(1:5),ylab=NULL,main=NULL,ask=prod(par("mfcol"))<length(which)&&dev.interactive(),
subSET,add.unadjusted=FALSE,ident=FALSE,caption=c("Residuals vs Fitted","Residuals vs Time",
"ACF residuals","PACF residuals","Individual mean profiles"),cex.caption=1)
{
if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
       stop("'which' must be in 1:6")


show <- rep(FALSE, 6)
    show[which] <- TRUE
 

 x1<-x@Fitted.av
 x2<-x@Fitted #new
 r<-x@s.residuals
 r1<-x@residuals	#all residuals
 
  

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

 time<-x@Time #new
 x3<-unique(x@Time)
 n.time<-length(x3)

 x4<-unique(x@model.matrix)
 y.av<-x@y.av
 n.id<-x@n.cases
 data<-x@subset.data

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

 	ylab1<-"Standardized residuals"
	xlab1<-"Fitted values"
        xlab2<-"Time"
 
if (any(show[1:6])){
         if(length(ylab)==0) ylab="probability"
         else ylab<-ylab
         if(length(main)==0) main=" "
         else main<-main
        }
 
  if(length(which)==5)
 	{ 

            	one.fig<-prod(par("mfcol"))==1
             
	     	if(ask){
             	oask<-devAskNewPage(TRUE)
            	on.exit(devAskNewPage(oask))
            	 }

##	 Fitted vs residuals
 
		ylim<-range(r,na.rm=TRUE)
             	ylim<-extendrange(r=ylim,f=0.08)
          	plot(x1,r,xlab=xlab1,ylab=ylab1,ylim=ylim)
 		abline(h=0,lty=3,col="gray")
             	mtext(caption[1],3,0.25,cex=cex.caption)
 
##	 Time vs residuals
 
		ylim<-range(r,na.rm=TRUE)
             	ylim<-extendrange(r=ylim,f=0.08)
 		plot(x3,r,pch=" ",xlab=xlab2,ylab=ylab1,ylim=ylim)
 		for(i in 1:n.time) {
 		lines(c(x3[i],x3[i]),c(0,r[i]))
 		}
 		abline(h=0,lty=3,col="gray")
        	mtext(caption[2],3,0.25,cex=cex.caption)
 
##	ACF 
 		acf(r, main=" ")
		mtext(caption[3],side=3,0.25, cex=1)

##	PACF 
 		pacf(r, main=" ")
 		mtext(caption[4],side=3,0.25, cex=1)

##	Parametric fit & Adjusted parametric fit & unadjusted parametric fit
	 
 		coef<-x@coefficients
 		npar <-length(coef)
 		omega<-as.double(coef[npar])
 
 		n.col<-ncol(x4)
 		coef1<-as.double(x@coefficients[1:n.col])
 		eta<-x4%*%coef1
         
 		c2<-(16*sqrt(3)/(15*pi))^2
 		eta.c<-eta*(c2*exp(omega)+1)^(-0.5)
 
 		n.time<-length(x3)
 
 		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab2,ylab=ylab,main=" ")
                mtext(main,side=3,0.25, cex=1)

		n1<-1
		n2<-n.time
		
                
		if(depend=="ind"||depend=="MC1"||depend=="MC2"){#NOVO
                

		for(j1 in 1:ncurves){
                n2<-j1*(n.time)
		lines(x3,plogis(eta[n1:n2]),col=j1,lty=1)#
          	n1<-n2+1
		}
		legend(x3[1],1,c("parametric fit"),bty="n",cex=0.75)
		legend(x3[1],0.95,paste(trace.label),bty="n",cex=0.75)#


		if(ncurves>1){legend(x3[1],0.90,levg,
 		lty=c(1,1),col=colg,bty="n",cex=0.75)}

    								 }#NOVO

		if(depend=="indR"||depend=="MC1R"||depend=="MC2R"){#Novo1
       		
		if(add.unadjusted==FALSE){

		n1<-1
		n2<-n.time
		for(j1 in 1:ncurves){                
		n2<-j1*(n.time)
		lines(x3,plogis(eta.c[n1:n2]),col=j1,lty=1)
          	n1<-n2+1
		
		                    }
 		
		y.time<-max(x3)-0.4*(max(x3)-min(x3))
		legend(y.time,1,c("adjusted parametric fit"),bty="n",cex=0.75)
		legend(y.time,0.95,paste(trace.label),bty="n",cex=0.75)

		if(ncurves>1){

		legend(y.time,0.90,levg,lty=c(1,1),col=colg,bty="n",cex=0.75)		
		
                }

                else{
#		legend(x3[1],1,c("adjusted parametric fit"),lty=2,bty="n",cex=0.75)
#		legend(x3[1],0.95,paste(trace.label),bty="n",cex=0.75)##############
                  }##########
	                                   }

		
           
 		if(add.unadjusted==TRUE){

		
		n1<-1
		n2<-n.time
		y.time<-max(x3)-0.4*(max(x3)-min(x3))
		for(j1 in 1:ncurves){                
		n2<-j1*(n.time)
		lines(x3,plogis(eta[n1:n2]),col=j1,lty=2)
		lines(x3,plogis(eta.c[n1:n2]),col=j1,lty=1)
          	n1<-n2+1	
		                   }
                legend(x3[1],1,c("unadjusted parametric fit"),bty="n",cex=0.75)
		legend(y.time,1,c("adjusted parametric fit"),bty="n",cex=0.75)
		legend(x3[1],0.95,paste(trace.label),bty="n",cex=0.75)#
		legend(y.time,0.95,paste(trace.label),bty="n",cex=0.75)#

		if(ncurves>1){legend(x3[1],0.90,levg,
 		lty=c(2,2),col=colg,bty="n",cex=0.75)
		legend(y.time,0.90,levg,
		lty=c(1,1),col=colg,bty="n",cex=0.75)}

                else{legend(x3[1],0.90,c("    "),lty=2,col=colg,bty="n",cex=0.75)
		legend(y.time,0.90,c("    "),lty=1,col=colg,bty="n",cex=0.75)}

					}
                     }#NOVO1


        }

######################################################
 
##	 Fitted vs residuals
 
  else if(show[1]){
                
                ylim<-range(r,na.rm=TRUE)
             	ylim<-extendrange(r=ylim,f=0.08)
     		plot(x1,r,xlab=xlab1,ylab=ylab1,ylim=ylim)
     		abline(h=0,lty=3,col="gray")
             	mtext(caption[1],3,0.25,cex=cex.caption)

    		}
##	 Time vs residuals
 
  else if(show[2])
                {
		ylim<-range(r,na.rm=TRUE)
             	ylim<-extendrange(r=ylim,f=0.08)
 		plot(x3,r,pch=" ",xlab=xlab2,ylab=ylab1,ylim=ylim)
 		for(i in 1:n.time) {
 		lines(c(x3[i],x3[i]),c(0,r[i]))
 		                   }
 		abline(h=0,lty=3,col="gray")
         	mtext(caption[2],3,0.25,cex=cex.caption)

                }

##	ACF 
 

  else if(show[3]){

	         acf(r, main=" ")
 		 mtext(caption[3],side=3,0.25, cex=1)

                 }
    
##	PACF
 
  else if(show[4]){
        
                 pacf(r, main=" ")
 		 mtext(caption[4],side=3,0.25, cex=1)

                 }

##	Parametric fit & Adjusted parametric fit & unadjusted parametric fit

# ##############################################################
  else if(show[5]){

         	coef<-x@coefficients
 		npar <-length(coef)
 		omega<-as.double(coef[npar])

 		n.col<-ncol(x4)
 		coef1<-as.double(x@coefficients[1:n.col])
 		eta<-x4%*%coef1
        
 		c2<-(16*sqrt(3)/(15*pi))^2
 		eta.c<-eta*(c2*exp(omega)+1)^(-0.5)
 
 		n.time<-length(x3)
 
 		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab2,ylab=ylab,main=" ")
 		mtext(main,side=3,0.25, cex=1)

		n1<-1
		n2<-n.time


                if(depend=="ind"||depend=="MC1"||depend=="MC2"){#Novo
 
		for(j1 in 1:ncurves){                
		n2<-j1*(n.time)
		lines(x3,plogis(eta[n1:n2]),col=j1,lty=1)
          	n1<-n2+1	
		                   }
                                             
		legend(x3[1],1,c("parametric fit"),bty="n",cex=0.75)
		legend(x3[1],0.95,paste(trace.label),bty="n",cex=0.75)#
		if(ncurves>1){legend(x3[1],0.90,levg,
 		lty=c(1,1),col=colg,bty="n",cex=0.75)}

    								 }#NOVO

 		


                if(depend=="indR"||depend=="MC1R"||depend=="MC2R"){#Novo1
       		
		if(add.unadjusted==FALSE){

		n1<-1
		n2<-n.time
		for(j1 in 1:ncurves){
                n2<-j1*(n.time)
		lines(x3,plogis(eta.c[n1:n2]),col=j1,lty=1)
          	n1<-n2+1
				     }
 		
		y.time<-max(x3)-0.4*(max(x3)-min(x3))
		
		legend(y.time,1,c("adjusted parametric fit"),bty="n",cex=0.75)
		legend(y.time,0.95,paste(trace.label),bty="n",cex=0.75)#
		if(ncurves>1){
		legend(y.time,0.90,levg,lty=c(1,1),col=colg,bty="n",cex=0.75)

                }


                else{
#		legend(x3[1],1,c("unadjusted parametric fit"),lty=1,bty="n",cex=0.75)
#		legend(x3[1],1,c("adjusted parametric fit"),lty=2,bty="n",cex=0.75)
                     }
					 }


                if(add.unadjusted==TRUE){

		
		n1<-1
		n2<-n.time
		y.time<-max(x3)-0.4*(max(x3)-min(x3))
		for(j1 in 1:ncurves){                
		n2<-j1*(n.time)
		lines(x3,plogis(eta[n1:n2]),col=j1,lty=2)
		lines(x3,plogis(eta.c[n1:n2]),col=j1,lty=1)
          	n1<-n2+1	
		                   }
                legend(x3[1],1,c("unadjusted parametric fit"),bty="n",cex=0.75)
		legend(y.time,1,c("adjusted parametric fit"),bty="n",cex=0.75)
		legend(x3[1],0.95,paste(trace.label),bty="n",cex=0.75)#
		legend(y.time,0.95,paste(trace.label),bty="n",cex=0.75)#

		if(ncurves>1){legend(x3[1],0.90,levg,
 		lty=c(2,2),col=colg,bty="n",cex=0.75)
		legend(y.time,0.90,levg,
		lty=c(1,1),col=colg,bty="n",cex=0.75)}
                 

 		else{legend(x3[1],0.90,c("    "),lty=2,col=colg,bty="n",cex=0.75)
		legend(y.time,0.90,c("    "),lty=1,col=colg,bty="n",cex=0.75)}


					}
                     }#NOVO1
                 }#FIM 5
####################################################

  
##	Individual mean profiles

	if(show[6])  {
 
if (depend=="ind"||depend=="MC1"||depend=="MC2") 
       stop("dependence must be indR, MC1R or MC2R")

		plot(x3,c(0,rep(1,n.time-1)), type="n", xlab=xlab2,ylab=ylab,main=" ")
 		mtext(main,side=3,1, cex=0.9)
		mtext(caption[5],3,0.25,cex=0.9)
                

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
                i<-afit[afit$id.numb==id.numbj,2]#ident position
         
                
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

}#Aqui

 )




