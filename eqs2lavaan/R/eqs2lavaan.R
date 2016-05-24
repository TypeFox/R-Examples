eqs2lavaan <-
function(eqs,data=NULL)
{
	data=NULL
	options(warn=-1)
	if(length(data)==1)
	{
		data	<- c(data,data)
	}
	cd		<- FALSE
	if(str_detect(eqs,ignore.case(".out")))
	{
		info	<- out2lavaan(eqs)
		covi	<- info[[1]]
		desc	<- info[[2]]
		cd		<- info
	}
	eqs		<- readLines(eqs, n=-1L)
	x		<- grep("PAGE:",eqs)
	if(length(x)>0)
	{
		eqs	<- eqs[-x]
	}
	x		<- grep("TITLE:",eqs)
	if(length(x)>0)
	{
		eqs	<- eqs[-x]
	}
	
	eq	<- NULL
	kp	<- NULL
	rp	<- NULL

	# LABELS #
	l1	<- grep("/LAB",eqs)[1]+1
	l2	<- grep("/",eqs[l1:length(eqs)])[1]+l1-2
	tmp	<- str_trim(eqs[l1:l2])
	tmp <- unlist(str_split(tmp,";"))
	tmp <- tmp[which(nchar(tmp)>0)]
	tmp	<- str_trim(str_sub(tmp,str_locate(tmp," ")[,1]))
	V	<- unlist(str_split(tmp,"[=]"))
	W	<- V[seq(2,length(V),2)]
	V	<- V[seq(1,length(V),2)]

	# EQUATIONS #
	l1	<- grep("/EQU",eqs)[1]+1
	l2	<- grep("/",eqs[l1:length(eqs)])[1]+l1-2
	tmp	<- str_trim(eqs[l1:l2])
	tmp <- unlist(str_split(tmp,";"))
	tmp <- tmp[which(nchar(tmp)>0)]
	tmp	<- str_trim(str_sub(tmp,str_locate(tmp," ")[,1]))
	for(i in 1:length(tmp))
	{
		x 	<- str_trim(str_split(tmp[i],"[=]")[[1]])
		kp	<- c(kp,word(x[2],-1))
		x	<- str_trim(c(x[1],str_split(x[2],"[+]")[[1]]))
		x	<- x[-length(x)]
		lhs <- str_trim(str_replace(x[-1],"[*]",""))
		pre	<- as.numeric(str_sub(lhs,1,str_locate(lhs,"F")[,1]))
		pre[is.na(pre)] <- ""
		if(length(grep("V",x))==0)
		{
			pre[is.na(pre)] <- ""
			lhs <- str_join(pre,"*",lhs)
			xx	<- which(str_sub(lhs,1,1)=="*")
			if(length(xx)>0)
			{
				lhs[xx] <- str_replace(lhs[xx],"[*]","")
			}
			rhs	<- str_trim(str_replace(x[1],"*",""))
			rp	<- c(rp,list(c(lhs,rhs)))
			lhs <- str_join(lhs,collapse=" + ")
			eq	<- c(eq,str_join(lhs,rhs,sep="  ~ "))
		}
		if(length(grep("V",x))!=0)
		{
			x[1] <- str_replace(x[1],x[1],W[which(V==x[1])])
			lhs	 <- str_trim(str_sub(x[2],str_locate(x[2],"F")[1]))
			pre  <- as.numeric(str_sub(x[2],1,str_locate(x[2],"F")[1]-1))
			rhs	 <- str_trim(x[1])
			rp	 <- c(rp,list(rhs))
			if(!is.na(pre))
			{
				rhs	<- str_trim(str_join(pre,x[1],sep="*"))
			}
			eq	<- c(eq,str_join(lhs,rhs,sep=" =~ "))
		}
	}
	
	kp	<- str_c(kp," ")
	#rp	<- str_c(rp," ")

	# VARIANCES #
	l1	<- grep("/VAR",eqs)[1]+1
	l2	<- grep("/",eqs[l1:length(eqs)])[1]+l1-2
	tmp	<- str_trim(eqs[l1:l2])
	tmp <- unlist(str_split(tmp,";"))
	tmp <- tmp[which(nchar(tmp)>0)]
	for(i in 1:length(tmp))
	{
		x <- str_trim(str_sub(tmp[i],str_locate(tmp[i]," ")[1]))
		if(length(grep(" to ",x))>0)
		{
			x <- str_split(x,"to")[[1]]
			x <- str_trim(c(x[1],str_split(x[2],"[=]")[[1]]))
			if(x[3]=="*")
			{
				x1	<- which(kp==str_c(x[1]," "))
				x2	<- which(kp==str_c(x[2]," "))
				if(length(x1)>0 & length(x2)>0)
				{
					for(j in x1:x2)
					{
						eq <- c(eq,str_c(rp[[j]],rp[[j]],sep=" ~~ "))
					}
				}
			}
			if(is.numeric(x[3]))
			{
				x1	<- which(kp==str_c(x[1]," "))
				x2	<- which(kp==str_c(x[2]," "))
				if(length(x1)>0 & length(x2)>0)
				{
					for(j in x1:x2)
					{
						eq <- c(eq,str_c(rp[[j]],str_c(x[3],"*",rp[[j]]),sep=" ~~ "))
					}
				}
			}
		}
		if(length(grep(",",x))>0)
		{
			x <- str_split(x,",")[[1]]
			x <- str_trim(c(x[1],str_split(x[2],"[=]")[[1]]))
			x1	<- which(kp==str_c(x[1]," "))
			x2	<- which(kp==str_c(x[2]," "))
			if(length(x1)>0 & length(x2)>0)
			{
				if(x[3]=="*")
				{
					eq <- c(eq,str_join(rp[[x1]],rp[[x2]],sep=" ~~ "))			
				}
				if(is.numeric(x[3]))
				{
					eq <- c(eq,str_join(r[[x1]],
					str_c(x[3],"*",rp[[x2]]),sep=" ~~ "))			
				}
				if(length(grep("F",x))==0)
				{
					if(x[3]=="*")
					{
						eq <- c(eq,str_join(rp[[x1]],rp[[x2]],sep=" ~~ "))			
					}
					if(is.numeric(x[3]))
					{
						eq <- c(eq,str_join(rp[[x1]],
						str_c(x[3],"*",rp[[x2]]),sep=" ~~ "))			
					}
				}
				if(length(grep("F",x))>0)
				{
					if(x[3]=="*")
					{
						eq <- c(eq,str_join(x[1],x[2],sep=" ~~ "))			
					}
					if(is.numeric(x[3]))
					{
						eq <- c(eq,str_join(x[1],
						str_c(x[3],"*",x[2]),sep=" ~~ "))			
					}
				}
			}
		}
		if(length(x)<2)
		{
			x   <- str_trim(str_split(x,"[=]")[[1]])
			x1	<- which(kp==str_c(x[1]," "))
			if(length(x1)>0)
			{
				if(x[2]=="*")
				{
					eq <- c(eq,str_join(rp[[x1]],rp[[x1]],sep=" ~~ "))			
				}
				if(is.numeric(x[2]))
				{
					eq <- c(eq,str_join(rp[[x1]],
					str_c(x[2],"*",rp[[x1]]),sep=" ~~ "))			
				}
			}

			if(length(grep("F",x))>0)
			{
				if(x[2]=="*")
				{
					eq <- c(eq,str_join(x[1],x[1],sep=" ~~ "))			
				}
				if(is.numeric(x[2]))
				{
					eq <- c(eq,str_join(x[1],
					str_c(x[2],"*",x[1]),sep=" ~~ "))			
				}
			}
		}
	}

	# COVARIANCES #
	l1	<- grep("/COV",eqs)[1]+1
	l2	<- grep("/",eqs[l1:length(eqs)])[1]+l1-2
	tmp	<- str_trim(eqs[l1:l2])
	tmp <- unlist(str_split(tmp,";"))
	tmp <- tmp[which(nchar(tmp)>0)]
	for(i in 1:length(tmp))
	{
		x <- str_trim(str_sub(tmp[i],str_locate(tmp[i]," ")[1]))
		if(length(grep(" to ",x))>0)
		{
			x <- str_split(x,"to")[[1]]
			x <- str_trim(c(x[1],str_split(x[2],"[=]")[[1]]))
			if(x[3]=="*")
			{
				x1	<- as.numeric(str_sub(x[1],2))
				x2	<- as.numeric(str_sub(x[2],2))
				for(j in x1:x2)
				{
					for(k in j:x2)
					{
						eq <- c(eq,str_join(str_join("F",j),
						        str_join("F",k),sep=" ~~ "))
					}
				}
			}
			if(is.numeric(x[3]))
			{
				x1	<- as.numeric(str_sub(x[1],2))
				x2	<- as.numeric(str_sub(x[2],2))
				for(j in x1:(x2-1))
				{
					for(k in (j+1):x2)
					{
						eq <- c(eq,str_join(str_join("F",j),
						        str_join(x[3],"*F",k),sep=" ~~ "))
					}
				}
			}
		}
		if(length(grep(",",x))>0)
		{
			x  <- str_split(x,",")[[1]]
			x  <- str_trim(c(x[1],str_split(x[2],"[=]")[[1]]))
			x1 <- which(kp==str_c(x[1]," "))
			x2 <- which(kp==str_c(x[2]," "))
			if(length(grep("F",x))==0)
			{
				if(x[3]=="*")
				{
					eq <- c(eq,str_join(rp[[x1]],rp[[x2]],sep=" ~~ "))			
				}
				if(is.numeric(x[3]))
				{
					eq <- c(eq,str_join(rp[[x1]],
					str_c(x[3],"*",rp[[x2]]),sep=" ~~ "))			
				}
			}
			if(length(grep("F",x))>0)
			{
				if(x[3]=="*")
				{
					eq <- c(eq,str_join(x[1],x[2],sep=" ~~ "))			
				}
				if(is.numeric(x[3]))
				{
					eq <- c(eq,str_join(x[1],
					str_c(x[3],"*",x[2]),sep=" ~~ "))			
				}
			}
		}
	}
	for(i in 1:length(kp))
	{
		#eq	<- str_replace_all(str_c(eq," "),kp[i],rp[i])
	}
	eq	<- str_trim(unique(eq))

	if(length(cd)>1 & length(data)<=1)
	{
		loc	<- agrep("CASES=",eqs)[1]
		x	<- str_locate(eqs[loc],ignore.case("CASES"))[,2]
		loc	<- str_sub(eqs[loc],x+1)
		loc	<- str_replace(loc,"=","")
		loc	<- as.numeric(str_trim(str_replace(loc,";","")))
		info<- c(info,list(loc))
		M	<- as.vector(info[[2]][,1],"numeric")
		n	<- as.numeric(info[[3]])
		lav	<- lavaan(toupper(eq),sample.cov=covi,sample.mean=M,sample.nobs=n,mimic="EQS")
	}
	if(length(data)>1)
	{
		lav	<- lavaan(eq,data=data,mimic="EQS")
	}
	i	<- intersect(grep("CHI-SQUARE",eqs),grep("DEGREES",eqs))
	if(length(i)>0)
	{
		chsq	<- as.numeric(str_trim(str_split(eqs[i[2]]," ")[[1]]))
		chsq	<- chsq[which(!is.na(chsq))]
		if(abs(chsq[1]-lav@Fit@test[[1]]$stat)>0.3 || chsq[2]!=lav@Fit@test[[1]]$df)
		{
			lav		<- noquote("The results for this translation seem to not be measuring something correctly.  Please report the error and send a copy of the .out file to craigmk@my.uri.edu so this error can be corrected.")
		}
	}
	if(length(i)==0)
	{
		lav		<- noquote("EQS file does not have a Chi-Square and likely has a matrix that is not positive definite.")
	}
	return(lav)
}