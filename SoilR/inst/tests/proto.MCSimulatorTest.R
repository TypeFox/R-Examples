#
# vim:set ff=unix expandtab ts=2 sw=2:
test.MC=function(){
	require("parallel")
	tcn="time"				#name of time column
	attcn="averageTransitTime"		#name ot the columns
	extractColumn=function(df,colname){df[,colname]}
	reduce2singledf=function(l){
		colavg=function(colname){
			mat=mcmapply(extractColumn,l,MoreArgs=list(colname))
			rs=function(i){mean(mat[i,])}
			newcol=mcmapply(rs,1:nrow(mat))
			return(newcol)
		}
		colnames=setdiff(names(l[[1]]),tcn)
		res=cbind(l[[1]][,tcn],mapply(colavg,colnames))
		res=as.data.frame(res)
		names(res)<-names(l[[1]])
		return(res)
	}
	npcn="np"
	f=function(pseudoarg){
		np0=1e2					#number of particles initially in the pool
		ir=function(t){2e3}			#inputrate
		st=0					#Starttime
		k=1e-1					#decay constant
		nt=1e2					#number of timesteps 
		sts=5e-1				#timestepsize	
		etcn="entryTime"			#name of entryTime column
		#####################################################################
		resultline=function(t,np,aTT){
			rd=data.frame(t,np,aTT)
			names(rd)<-c(tcn,npcn,attcn)
			return(rd)
		}
		#####################################################################
		newParticles=function(t,n){#computes the amount of Input particles for this timestep and sets the entry time property
			# in this case constant
			npd=data.frame(entryTime=rep(t,n),color=rep("green",n))
			names(npd)<-c(etcn,"color")
			return(npd)
		}
		#####################################################################
		pd=newParticles(st,np0)
		results=resultline(st,np0,0)
		t=st
		for (i in 1:nt){
			pd=rbind(pd,newParticles(t,ir(t)*sts)) #add inputrate x timestepsize new particles with entry time t
			np=nrow(pd)
			pl <- k*sts 		#Probability for a single particle to leave during the present timestep 
						#in this case constant but in general time dependent and therefor inside the loop
			r <- runif(np,0,1.0)	#random value for each particle between 0 and 1
			ll=r<pl			#logical vector true for leaving particles
			ls=!ll			#logical vector  true for stayin particles
			lpd <- as.data.frame(pd[ll,])	#dataframe of leaving particles
			t=t+sts
			mtt=mean(t-lpd[,etcn])
			#print(nrow(lpd))
			#print(mtt)
			pd <- as.data.frame(pd[ls,])	#dataframe of particles still in the pool
			results=rbind(results,resultline(t,np,mtt))
		}	
		return(results)
	}
	r1=f(1)
	np=detectCores()
	dfl=mclapply(1:np,f,mc.cores=np)
	cr=reduce2singledf(dfl)
	plot(r1[,tcn],r1[,npcn])
	lines(cr[,tcn],cr[,npcn])
	plot(r1[,tcn],r1[,attcn])
	lines(cr[,tcn],cr[,attcn])
}

