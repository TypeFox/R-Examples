AS.test <-function(grid=10,x1,x2,supp=NULL){

#NOTE: this function should be called *after* choosing statistic subsets

## does Joyce and Marjoram's sufficient statistics odds-ratio computation
##
## grid is the grid of theta values on which to interpolate/construct posterior
## or number of points.
## x1 & x2 accepted thetas, (based on "subset1" (Sk-1) 
## and "subset2". "subset2" should be one more than "subset1" (Sk)).

x1<-as.matrix(x1)
x2<-as.matrix(x2)

Nkminus1<-length(x1)	#in our case, these two should be the same
Nk<-length(x2)		#

rstar<-Ti<-Ti2<-matrix(0,1,grid)
s<-0:grid		#makes sure that the number of bins stays at grid.

r1<-range(x1)
r2<-range(x2)
rL<-min(c(r1[1],r2[1]))
rR<-max(c(r1[2],r2[2]))
s<-seq(from=min(c(supp[1],rL)),to=max(c(supp[2],rR)),length=(grid+1))

tmp1<-hist(x1,breaks=s,plot=F)
tmp2<-hist(x2,breaks=s,plot=F)


Nkminus1counts<-tmp1$counts
Nkcounts<-tmp2$counts

expected<-matrix(Nkminus1counts*Nk/Nkminus1,nrow=1)

sd<-sqrt(expected*(Nkminus1-Nkminus1counts)/Nkminus1)

Ti<-expected+4*sd		
Ti2<-expected-4*sd		

violate1<-(Nkcounts>Ti)
violate2<-(Nkcounts<Ti2)

extreme<-(any(violate1)+any(violate2))>0	#Sk adds information 

#print(extreme)

return(extreme)

}

