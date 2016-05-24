DimSelector <-
function(out,n,m,p,model){

#size adjustments
eff_n=n
if (n>m*p){
	eff_n=m*p
}
eff_m=m
if (m>n*p){
	eff_m=n*p
}
eff_p=p
if (p>m*n){
	eff_p=m*n
}
#P,Q,R : maximum number of dimensions 
P=max(out[,1])
Q=max(out[,2])
R=max(out[,3])
# Sorting by total number of components (5th column of out)
out=out[ord(out[,5])$a,]
# Step 2 Retain only the best-fitting solutions
tmp1=nrow(out)
fpi=min(out[,5])
for (i in 1:tmp1){
    if (out[i,5] != max(fpi)){
		fpi=rbind(fpi,out[i,5])	
    }
}
tmp3=length(fpi)
out_best=matrix(0,tmp3,5)
for (i in 1:tmp3){
    for (j in 1:tmp1){
        if (fpi[i]==out[j,5]){
            if (out_best[i,4]<out[j,4]){
               out_best[i,]=out[j,]
            }
        }
    }
}
# Step 4 Clean for f
out_best_best =out_best[1,]
for (i in 2:tmp3){
    if (max(out_best_best[4])< out_best[i,4]){
        out_best_best=rbind(out_best_best, out_best[i,])
    }
}
# Steps 5 and 6 Testing triplets of adjacent solutions
tmp4=nrow(out_best_best)
out_best_best_best=out_best_best[1,]
i=2
while (i < (tmp4-1)){
    f1=out_best_best[(i-1),4]
    f2=out_best_best[i,4]
    f3=out_best_best[(i+1),4]
	#Tucker3,Tucker2,Tucker1
	#Candecomp/Parafac
	if (model==1){ 
		fp1=(eff_n+eff_m+eff_p-2)*out_best_best[(i-1),1]
		fp2=(eff_n+eff_m+eff_p-2)*out_best_best[i,1]
		fp3=(eff_n+eff_m+eff_p-2)*out_best_best[(i+1),1]
	} else{ 
		fp1=eff_n*out_best_best[(i-1),1]+eff_m*out_best_best[(i-1),2]+eff_p*out_best_best[(i-1),3]+out_best_best[(i-1),1]*out_best_best[(i-1),2]*out_best_best[(i-1),3]-out_best_best[(i-1),1]^2-out_best_best[(i-1),2]^2-out_best_best[(i-1),3]^2
		fp2=eff_n*out_best_best[i,1]+eff_m*out_best_best[i,2]+eff_p*out_best_best[i,3]+out_best_best[i,1]*out_best_best[i,2]*out_best_best[i,3]-out_best_best[i,1]^2-out_best_best[i,2]^2-out_best_best[i,3]^2
		fp3=eff_n*out_best_best[(i+1),1]+eff_m*out_best_best[(i+1),2]+eff_p*out_best_best[(i+1),3]+out_best_best[(i+1),1]*out_best_best[(i+1),2]*out_best_best[(i+1),3]-out_best_best[(i+1),1]^2-out_best_best[(i+1),2]^2-out_best_best[(i+1),3]^2
	}	
	if (LineCon(f1,f2,f3,fp1,fp2,fp3)==0){
		out_best_best=rbind(out_best_best[1:(i-1),],out_best_best[(i+1):tmp4,])
		tmp4=nrow(out_best_best)
		i=1
    }
	i=i+1
}
# Step 7 Compute St
tmp4=nrow(out_best_best)
st=cbind(vector(mode="numeric",length=tmp4))
for (i in 2:(tmp4-1)){
    fi=out_best_best[i,4]
    fi_p=out_best_best[(i-1),4]
    fi_n=out_best_best[(i+1),4]
	fpi_p=eff_n*out_best_best[(i-1),1]+eff_m*out_best_best[(i-1),2]+eff_p*out_best_best[(i-1),3]+out_best_best[(i-1),1]*out_best_best[(i-1),2]*out_best_best[(i-1),3]-out_best_best[(i-1),1]^2-out_best_best[(i-1),2]^2-out_best_best[(i-1),3]^2
	fpi=eff_n*out_best_best[i,1]+eff_m*out_best_best[i,2]+eff_p*out_best_best[i,3]+out_best_best[i,1]*out_best_best[i,2]*out_best_best[i,3]-out_best_best[i,1]^2-out_best_best[i,2]^2-out_best_best[i,3]^2
	fpi_n=eff_n*out_best_best[(i+1),1]+eff_m*out_best_best[(i+1),2]+eff_p*out_best_best[(i+1),3]+out_best_best[(i+1),1]*out_best_best[(i+1),2]*out_best_best[(i+1),3]-out_best_best[(i+1),1]^2-out_best_best[(i+1),2]^2-out_best_best[(i+1),3]^2
	
	if (st[i-1]<Inf){
		st[i]=((fi-fi_p)/(fpi - fpi_p))/((fi_n - fi)/(fpi_n - fpi))
	} else{
		st[i]=Inf
	}
}
out_st=cbind(out_best_best,st)

# Step 8 Select the solution with highest st-value
imax=SUM(out_st)$valueMaxc[6]	
st_max = out_st[imax,]

# Output

# Table 
if (model==1){
cat("Table: Goodness-of-fit values (%) and Scree test values st",fill=TRUE)
cat("of the solutions on the higher boundary of the convex hull",fill=TRUE)
} else{
cat("Table: Goodness-of-fit values (%), Total number of components tnc and Scree test",fill=TRUE)
cat("values st of the solutions on the higher boundary of the convex hull",fill=TRUE)
}

#CANDECOMP/PARAFAC
if (model==1){
	tab=matrix(0,nrow=tmp4,ncol=2)
	for (i in 1:tmp4){
		tab[i,1]=noquote(paste("     (",out_st[i,1],") "))
		if (i==1){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"      -"))
		}
		if (i==tmp4){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"      -"))
		} else{
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",round(out_st[i,6],digits=2)))
		}
		if (i==imax){
			suggestion=noquote(paste("The hull heuristic indicates the selection of the CP model of complexity (",out_st[i,1],")"))
		}
	}
	rownames(tab)=rep("CP",length=tmp4)
	colnames(tab)=c("   Complexity","   Fit (%)      st")
	print(tab,quote=FALSE)
	cat(suggestion,fill=TRUE)
}
#TUCKER3
if (model==2){
	tab=matrix(0,nrow=tmp4,ncol=2)
	for (i in 1:tmp4){
		tab[i,1]=noquote(paste("   (",out_st[i,1],out_st[i,2],out_st[i,3],") "))
		if (i==1){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		}
		if (i==tmp4){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		} else{
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"    ",round(out_st[i,6],digits=2)))
		}
		if (i==imax){
			suggestion=noquote(paste("The hull heuristic indicates the selection of the T3 model of complexity (",out_st[i,1],out_st[i,2],out_st[i,3],")"))
		}
	}
	rownames(tab)=rep("T3",length=tmp4)
	colnames(tab)=c("   Complexity","   Fit (%)      tnc    st")
	print(tab,quote=FALSE)
	cat(suggestion,fill=TRUE)
}
#TUCKER2
if (model==3){
	tab=matrix(0,nrow=tmp4,ncol=2)
	for (i in 1:tmp4){
		tab[i,1]=noquote(paste("   (",out_st[i,1],out_st[i,2],out_st[i,3],") "))
		if (i==1){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		}
		if (i==tmp4){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		} else{
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"    ",round(out_st[i,6],digits=2)))
		}
		if (i==imax){
			suggestion=noquote(paste("The hull heuristic indicates the selection of the T2 model of complexity (",out_st[i,1],out_st[i,2],out_st[i,3],")"))
		}
	}
	rownames(tab)=rep("T2",length=tmp4)
	colnames(tab)=c("   Complexity","   Fit (%)      tnc    st")
	print(tab,quote=FALSE)
	cat(suggestion,fill=TRUE)
}
#TUCKER1
if (model==4){
	tab=matrix(0,nrow=tmp4,ncol=2)
	for (i in 1:tmp4){
		tab[i,1]=noquote(paste("   (",out_st[i,1],out_st[i,2],out_st[i,3],") "))
		if (i==1){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		}
		if (i==tmp4){
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"      -"))
		} else{
			tab[i,2]=noquote(paste("  ",round(out_st[i,4],digits=4),"    ",out_st[i,5],"    ",round(out_st[i,6],digits=2)))
		}
		if (i==imax){
			suggestion=noquote(paste("The hull heuristic indicates the selection of the T1 model of complexity (",out_st[i,1],out_st[i,2],out_st[i,3],")"))
		}
	}
	rownames(tab)=rep("T1",length=tmp4)
	colnames(tab)=c("   Complexity","   Fit (%)      tnc    st")
	print(tab,quote=FALSE)
	cat(suggestion,fill=TRUE)
}
}
