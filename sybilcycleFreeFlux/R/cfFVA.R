cfFVA <- function(model,rxnList,solver = SYBIL_SETTINGS("SOLVER") ){
	#when a loop exists set rxn(s) of no essential flux to 0: except rxn to be maximized 
	# if more than one rxn can be a problem!!!
	## send obj value to lrFBA
#function to get reaction equation
	getRxnEqn <- function(rxn=1){
  	m=S(model)[,rxn]
  	# input ==> output
  	mcf=ifelse(abs(m)>1,paste("(",abs(m),")",sep=""),"")
  	eqn=paste(gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m<0)],sep=""))),
  		ifelse(react_rev(model)[rxn],"<==>","-->"),
  	    gsub(" "," + ",Reduce(paste,paste(mcf[m<0],met_id(model)[which(m>0)],sep=""))) ,sep=" ")
  	
  return(eqn)
}

	sf=optimizeProb(model,solver=solver);
	bmrxn=which(obj_coef(model)==1);
	objVal=lp_obj(sf)#sf$obj;
	
#####-------------------------########----------------------#######
maxFlx=NULL;
lowbnd(model)[bmrxn]=objVal;
uppbnd(model)[bmrxn]=objVal;

#excReact = findExchReact(model)[1];# 1 is position
#excReactPos=react_pos(excReact$exchange);
excReact = findExchReact(model);
excReactPos=which(react_id(model) %in% react_id(excReact));
	
#mod=sysBiolAlg(model,solver=solver);
#prob_r=problem(mod)
model_r=model

nCols=react_num(model);
# for tracing^^
rxneqns=sapply(c(1:nCols),function(x) getRxnEqn(x))
res=NULL;
for (rxn in rxnList){
	for(dirxn in c(1,-1)){
	#mod=sysBiolAlg(model,solver=solver);
    #prob=problem(mod)
    model=model_r
	sn=which(react_id(model)==rxn);
	#ocf=obj_coef(model)
	#backupProb(mod) to restore original copy, other solver may have time issues
	lb=lowbnd(model);	ub=uppbnd(model);
	#changeColsBnds(prob,j=c(1:nCols),lb=lb,ub=ub)
	for(i in 1:5){# max no of loops
		print(c("rxn:",rxn," iteration:",i));	
		#1- Max flux through rxn
		obj_coef(model)[bmrxn]=0;		obj_coef(model)[sn]=dirxn;#min/max
		#solL=optimizeProb(model,solver=solver);# loopy max
		
		#ocf[bmrxn]=0; ocf[sn]=1;
		#solL=optimizeProb(mod,react=c(1:nCols),obj_coef=ocf,resetChanges =TRUE)
		#changeObjCoefs(prob,c(sn,bmrxn),c(dirxn,0))
		#solL=optimizeProb(mod);#,solver=solver:solver defined in opt_obj
		solL=optimizeProb(model,solver=solver)
		F1=getFluxDist(solL)#$fluxes#[solL$fldind];
		#print(c("i",i))
		#print(F1[react_id(model)=="R_EX_o2_e_"])
					
		if(i==1)F1org=F1;
		#model=model_r;
		# 2-F1 without loops
		obj_coef(model)[sn]=0;		obj_coef(model)[bmrxn]=1;
		#changeObjCoefs(prob,c(sn,bmrxn),c(0,1))
		#print(getObjDir(prob))
		if(abs(F1[sn])>1e-7){
			if(!(sn %in% excReactPos)){#send xchng rxn list to avoid recalc.
					solNL=cfFBA(model,wtflux=F1,excReactPos=excReactPos,solver=solver,objVal=objVal,retOptSol=FALSE);	
			}else{#fix exchange except one rxn
			        excReactPos1=excReactPos[-which(sn %in% excReactPos)]
					solNL=cfFBA(model,wtflux=F1,excReactPos=excReactPos1,fixExchRxn=TRUE,solver=solver,objVal=objVal,retOptSol=FALSE);	
			}
			F2=solNL$fluxes;#F2 represents F1 without loops
		}else{F2=F1}
			
		#3-remove all loops except the one containing [rxn](Identify loop)
		if( abs(F1[sn]-F2[sn]) >0.01){# if a loop exists
			print(c("diff:" , F1[sn],F2[sn])) ;
			#lowbnd(model)[sn]=F1[sn];	uppbnd(model)[sn]=F1[sn]; # Force rxn to go with max flux to identify loop
			if(!(sn %in% excReactPos)){#fix rxn sn
			         sol1L=cfFBA(model,wtflux=F1,fixExchRxn=TRUE,excReactPos=c(sn,excReactPos),solver=solver,objVal=objVal,retOptSol=FALSE);
			}else{
					excReactPos1=excReactPos[-which(sn %in% excReactPos)]
					sol1L=cfFBA(model,F1,excReactPos=excReactPos1,solver=solver,objVal=objVal,retOptSol=FALSE);
			}
			print(c("sol 3:",sol1L$stat,sol1L$obj));
			F3=sol1L$fluxes#[sol1L$fldind];Only one Loop
			#write.csv(file=sprintf("fl_iter%d_%s.csv",i,rxn),cbind(rxn=react_id(model),rxneqns,lb,ub,F1,F2,F3))
			
			lp=(abs(F3-F2) >0.01);## detected loop (same flx rot?)
			#flx_rot=abs(F1[sn]-F2[sn])
			#lp=(abs(F3-F2) >=flx_rot)
			print("Rxns in loop: ");
			print(cbind(which(lp),react_id(model)[lp],F3=F3[lp],F2=F2[lp],sFVA=F1[lp]));
			#4- break loop and maximize again
			# set rxns going to 0 to 0 and maximize again
			#model=model_r; undo last change in model but keep changes thru the crnt rxn iteration
			#lowbnd(model)[sn]=lowbnd(model_r)[sn]; uppbnd(model)[sn]=uppbnd(model_r)[sn];
			brk=(abs(F3-F2) >0.01 & (abs(F2) < 1e-7)); # rxn in loop goes to zero
			print(min(abs(F3[lp])))
			if(brk[sn] && sum(brk)==1) {# crnt rxn the only rxn having Zero flux without loop
				maxFlx=rbind(maxFlx,cbind(rxn,iter=i,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp=paste(react_id(model)[lp],collapse=",")
				,brk=paste(react_id(model)[brk],collapse=","),
				loopyflx=paste(F1[lp],collapse=","),NLflx=paste(F2[lp],collapse=","),LPF3flx=paste(F3[lp],collapse=","),
				rvrs=paste(react_rev(model)[lp],collapse=","),flxrot=abs(F1[sn]-F2[sn]),llen=sum(lp),eqn=paste(rxneqns[lp],collapse=",") ));
				break;
			}
			brk[sn]=FALSE;
			if(sum(brk)==0){
				print("No rxn will go to zero in this Loop!");
				break;
			}
			#brk=lp;
			#, exclude rxn itself !
			#print(react_id(model)[brk]);
			
			maxFlx=rbind(maxFlx,cbind(rxn,iter=i,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp=paste(react_id(model)[lp],collapse=",")
			,brk=paste(react_id(model)[brk],collapse=","),
			loopyflx=paste(F1[lp],collapse=","),NLflx=paste(F2[lp],collapse=","),LPF3flx=paste(F3[lp],collapse=","),
			rvrs=paste(react_rev(model)[lp],collapse=","),flxrot=abs(F1[sn]-F2[sn]),llen=sum(lp),eqn=paste(rxneqns[lp],collapse=",") ));
			
			lb[brk]=ub[brk]=F2[brk];
			#changeColsBnds(prob,j=c(1:nCols),lb=lb,ub=ub)
			lowbnd(model)=lb;uppbnd(model)=ub;#Change bounds
	#		lowbnd(model)[brk]=F2[brk]; uppbnd(model)[brk]=F2[brk];
		}else {		 # max value is ok (there is no loop
		 	if (i>1) {print(c(rxn,i,F1[sn],F2[sn]));
			#write.csv(file=sprintf("fl_iter%d_%s.csv",i,rxn),cbind(rxn=react_id(model),rxneqns,lb,ub,F1,F2,F3))
			}
			maxFlx=rbind(maxFlx,cbind(rxn,iter=i,dirxn,initflx=F1org[sn],LoopyMax=F1[sn],NLMax=F2[sn],lp="",brklp="",loopyflx=NA,NLflx=NA,
			LPF3flx="",rvrs=react_rev(model)[sn],flxrot=0,llen=0,eqn=""));
			break;
		}
	}#iteration for loop
	if(dirxn==1){maxVal=F2[sn];sFVAmax=F1org[sn];maxiter=i}
	else {minVal=F2[sn]; sFVAmin=F1org[sn];miniter=i}
 }#dirxn
 res=rbind(res,cbind(rxn=rxn,maxiter,miniter,minVal,maxVal,sFVAmin,sFVAmax))
 }#rxn
 #by():min,max for each rxn
 #
return(list(res,maxFlx));
}