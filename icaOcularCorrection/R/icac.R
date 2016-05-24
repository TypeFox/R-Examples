icac <- function(x, # a data frame or icac object
         channel,
         noise.sig, # noise against which to threshold, can be a vector
		 trial.cn="Trial",
		 include=TRUE, # whether to include noise channels or not in ICA
         threshold=0.4,
         n.comp=length(channel),
         ica.method="R", # or "C"
		 correct=TRUE,
		 ica.only=FALSE,
         proctime = TRUE,
		 seed=NULL, # or an interger
         verbosity=5,
		 ...
){
    #       initialize proc.time calculation       #
    if(proctime){
        if(verbosity>3)cat("cleaning up ...\n")
		sink(file=file.path(tempdir(),"tmp.txt"))
        gc(verbose=FALSE,reset=TRUE)
		sink(file=NULL)
    	ptm1<-proc.time()
    }

    # some other preliminaries
	if(!is.null(seed)){
    	set.seed(seed)
	}
    require(fastICA,quietly=TRUE)
    correction="by.trial"
	channels<-NULL
	col.means<-NULL
	X<-NULL

	options(warn=-1)

	# perform ICA
	if(class(x)[1]=="data.frame"){
		channels<-channel
		
		if(!include){
			for(i in noise.sig){
				if(i%in%channels){
					channels<-channels[-which(channels==i)]
				}
			}
		}
		X<-x[,channels] # select channels only
		if(verbosity>0)cat("fastICA ...\n")

        # getting channel means
        col.means<-colMeans(X)
        names(col.means)<-channels
  
		# perform ICA decomposition (with pre-whitening and centering)
		a<-fastICA(X,n.comp=n.comp,method=ica.method,
                    verbose=ifelse(verbosity>1,TRUE,FALSE),...)

		if(ica.only)return(invisible(a))
	}

	# by trial correction
	# note that this is the only option now
	if(correction[1]=="by.trial"){ 
		S0<-a$S
		if(!trial.cn%in%colnames(x)){
			stop("no trial column. Please add one.\n")
		}
		# using Trial column
		S0<-as.data.frame(S0,stringsAsFactors=FALSE)
		S0[,trial.cn]<-x[,trial.cn]
		trials<-unique(S0[,trial.cn])
		cor.info<-vector("character")
        myinfo<-data.frame(NoiseSignal=vector("character"),
			IC=vector("numeric"),Trial=vector("numeric"),
			Corr=vector("numeric"))
		for (l in noise.sig){
			if(verbosity>0)cat(paste("correcting based on noise signal",
				l,"...\n"))
			pb<-txtProgressBar(min=1,max=length(trials),
				char="=",style=3)		
			# initialize corrected source matrix
			S.corrected<-a$S[-(1:nrow(a$S)),] 
			for(trial in trials){
				if(verbosity>1)setTxtProgressBar(pb,which(trials==trial))
				# set-up temporary source matrix with nrow = one epoch
				S.temp<-S0[S0[,trial.cn]==trial,] 
                S.temp<-S.temp[,-which(colnames(S.temp)==trial.cn)]
                S.temp<-as.matrix(S.temp)
				# set-up temporary x with nrow = one epoch
				x.temp<-x[x[,trial.cn]==trial,]
				for (m in 1:ncol(S.temp)){
					# calculate correlation between each 
					# IC and a given noise sig for a given epoch
					# and store it in variable 'corelations'.
					if(mean(S.temp[,m])!=0)correlation<-cor(x.temp[,l],
						S.temp[,m]) 
					citmp<-paste("noise signal = ",l,"; trial = ",trial,
						"; IC = ",m,"; cor = ",correlation,sep="")
					cor.info<-c(cor.info,citmp)
					if(is.na(correlation))correlation=0 
					if(abs(correlation)>=threshold){
                    	myinfo<-rbind(myinfo,data.frame(NoiseSignal=l,IC=m,
							Trial=trial,Corr=correlation))
						S.temp[,m]<-0
					}
				}
				# add temporary source matrix to corrected source variable
				S.corrected=rbind(S.corrected,S.temp)
			}
			close(pb)
			# replace original source matrix in x 'a' 
			# with corrected source matrix
			a$S=S.corrected
		}
        S0<-S0[,-which(colnames(S0)==trial.cn)]
        S0<-as.matrix(S0)
	}

	if(correct){
		if(verbosity>0)cat("re-mixing ICs ...\n")
		pb<-txtProgressBar(min=1,max=length(channels),char="=",
			style=3)		
		# re-mix corrected data = corrected source matrix * mixing matrix
		tmp=as.data.frame(a$S%*%a$A)
		# give appropriate column names to corrected frame
		colnames(tmp)=channels 
	  
		# replace EEG columns in original data frame with corrected ones
	        # and add channel mean (i.e., de-mean-center channels)
		for(channel in channels){ 
			if(verbosity>1)setTxtProgressBar(pb,which(channels==channel))
			x[,channel]=tmp[,channel]+col.means[channel]
		}
		close(pb)
	}

    rownames(myinfo)<-1:nrow(myinfo)
    myinfo<-myinfo[,c("NoiseSignal","IC",trial.cn,"Corr")]

	options(warn=0)

	# create icac object
	res<-list(data=x,channel=channels,
		noise.sig=noise.sig,threshold=threshold,n.comp=n.comp,
		X=a$X,K=a$K,W=a$W,A=a$A,S=a$S,S0=S0,col.means=col.means,
		correlations=cor.info,correction.info=myinfo)

	class(res)<-"icac"

    #    create proc.time table   #
    # The values presented (user, system, and elapsed) will be 
    # defined by your operating system, but generally, the user 
    # time relates to the execution of the code, the system time 
    # relates to your CPU, and the elapsed time is the difference 
    # in times since you started the stopwatch (and will be equal 
    # to the sum of user and system times if the chunk of code was 
    # run altogether).
    if(proctime){
          ptm2<-proc.time()
          ptm3<-ptm2-ptm1
          ptm4<-as.data.frame(rbind(ptm1[1:3],ptm2[1:3],ptm3[1:3]))
          ptm4<-ptm4/60
          colnames(ptm4)<-c("user","system","elapsed")
          rownames(ptm4)<-c("start time","end time","run time")
          res$proctime<-ptm4
          if(verbosity>3)print(ptm4)
    }
        
	# return output
	return(invisible(res))
}

