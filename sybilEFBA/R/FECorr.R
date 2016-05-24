################################################
# Function: FECorr: Flux Expression Correlation
 # uses FVA under different conditions to find fluxes that linearly correlates to corresponding gene expression. 

 FECorr <- function (model, nCond, initCond, geneExpressionData=NULL,RuleExpressionData=NULL,
				  pct_objective=100,#Optional:contains list of objective values to be fixed at each condition
				  selected_rxns=NULL,#optional used to select a set of reactions not all, Boolean with the same length react_id(model)
				  only_identified_rules=FALSE,# ignore rxns containing genes with unidentified expression
				  minExprFoldChange=0,#fold change between min expression level and max expression level (applied to rules) min(expr)*minEFC < max(Expr)
				  lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                 		 solver = SYBIL_SETTINGS("SOLVER"),
                   		method = SYBIL_SETTINGS("METHOD"),
                                #solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM")
				  solverParm=data.frame(CPX_PARAM_EPRHS=1e-6),
                  verboseMode = 2){
#PARAMETERS:
#===========
# model                 Sybil model structure (class modelorg)
# nCond					Number of conditions
# geneExpressionData        a data frame: GeneID,cond_id, expr_val
		#column rows are genes and column j+1 is representing gene expression under condition j
## RuleExpressionData   : rxn_id,cond_id, ExpressionVal		
# initCond: rxn_id,cond_id,lb,ub,objcoef : lower and upper bounds of rxns under different conditions that represents 
#        the available nutrients under these conditions
# pct_objective 
#minExprFoldChange: can be used to consider only gene with a significant change in expression level

# RETURN
# geneID,slope,intercept,base_level: OGOR: one gene one rxn
# iFlux: rxn_id,cond_id,lb,ub,objCoef,xpc_flux,fva_min,fva_max,RuleExprVal,iflux

# Challenges
# 1-rxn_enz: rxn_id, enzyme, cond_id, flux,expression
# 2-isoenzymes,multirxn catalyzers,protein complexes are problems
# 3-reversible reactions: in FVA, use absolute flux
##-------------------------------------------------------------------------------##

#Main steps
    #1- Run FVA for all conditions, exclude rxns fixed in all conditions
    #2- Identify ruleExpression for set of rxns remaining from 1,
    #3- Fit Expr to FVA range.
	#4- Run findMDCFlux to find closest genome-scale flux
    #5- Recalculate correlation: Posterior, iFlux, ruleExpr

##--------------------------------------------------------------------------##
init_model<- function(p_cnd){
	model=model_r;
	cnd_cons=initCond[initCond[,"cond_id"]==p_cnd,]
	if(verboseMode>3) print(c(p_cnd," bnd: ",length(cnd_cons[,1])))
	bnd_ind=match(cnd_cons[,"rxn_id"],react_id(model) ) 
	lowbnd(model)[bnd_ind]=as.numeric(cnd_cons[,"lb"])
	uppbnd(model)[bnd_ind]=as.numeric(cnd_cons[,"ub"])
	obj_coef(model)[bnd_ind]=as.numeric(cnd_cons[,"obj"])
	
	return(model);
}
##--------------------------------------------------------------------------##
# check prerequisites 
    if (!is(model, "modelorg")) {
      stop("An object of class modelorg is needed!")
    }

	if(length(RuleExpressionData)!=0 ){
		cnds=unique(as.vector(RuleExpressionData[,"cond_id"]));
	}else  
		if(length(geneExpressionData)!=0){
			cnds=unique(as.vector(geneExpressionData[,"cond_id"]));
		}else{
			stop("geneExpressionData or RuleExpressionData must be not null!");
		}
##-------------------------------------------------------------------------------##
# keep original copy of model
	model_r=model;
	nc=react_num(model);
	
		if(length( selected_rxns)==0){
                  gpr_rxn_ind=(gpr(model)!="" ) # only rxns having GPR rule
         }else { gpr_rxn_ind=selected_rxns;}

#1- Run FVA for all conditions, exclude rxns fixed in all conditions
	all_fva=NULL;
	rxnfva=gpr_rxn_ind;
	nrxnfva=sum(rxnfva)
  #cnds_init=unique(as.vector(initCond[,"cond_id"]));

    for(cnd in cnds){
	    model=init_model(cnd)
		
		fv <- fluxVar(model,react=which(rxnfva),solver=solver,percentage=pct_objective);#solver=solver, error LP_SOLVER		
		fv_min=lp_obj(fv)[1:sum(rxnfva)];fv_max=lp_obj(fv)[(sum(rxnfva)+1):(2*sum(rxnfva))]
		
		all_fva=rbind(all_fva,cbind(cond_id=rep(cnd,nrxnfva),serial=which(rxnfva),rxn_id=react_id(model)[rxnfva],fv_min,fv_max));
	}
	all_fva=as.data.frame(all_fva)
	all_fva[,"fv_min"]=as.numeric(as.vector(all_fva[,"fv_min"]))
	all_fva[,"fv_max"]=as.numeric(as.vector(all_fva[,"fv_max"]))
	#exclude fixed (blocked) rxns
	blocked_rxns=rep(FALSE,nc)
	
	for(r in which(gpr_rxn_ind)){
		umin=unique(round(all_fva[all_fva[,"serial"]==r,"fv_min"],6))
		umax=unique(round(all_fva[all_fva[,"serial"]==r,"fv_max"],6))
		if(length(umin)==1 && length(umax)==1 && (umin[1]==umax[1])) blocked_rxns[r]=TRUE;#gpr_rxn_ind[r]=FALSE;
	}
	gpr_rxn_ind[blocked_rxns]=FALSE;
	if(verboseMode>2)	print(sprintf("No of non fixed rxns remaining after FVA: %d, number of blocked rxns: %d",sum(gpr_rxn_ind),sum(blocked_rxns)));
	#print(all_fva)
	
   if(verboseMode>3){
	write.csv(file=format(Sys.time(), "all_fva_%Y%m%d_%H%M.csv"),all_fva);
	}
##-------------------------------------------------------------------------------##
#2- Identify ruleExpression for set of rxns remaining from 1,
# one gene to one rxn
# isoenzymes
# protein complex
# others: multiple complex rule
gs=RuleExpressionData
if (length(gs)==0){
	#process geneExpressionData to get RuleExpression
	print("Processing geneExpressionData to get RuleExpression")
	gs=NULL;
	ag=allGenes(model)
	for(cnd in cnds){		
		v_selected_rxns=gpr_rxn_ind;
		geneExpr=geneExpressionData[geneExpressionData$cond_id==cnd,];
		if(only_identified_rules){
			for(r in which(v_selected_rxns)){
				rxnGenes=ag[(rxnGeneMat(model)[r,])]
				if(sum(!(rxnGenes %in% geneExpr$GeneID))>0){
					v_selected_rxns[r]=FALSE;
					if(verboseMode>3) {
						print(sprintf("Rule %s of rxn: %d has unidentified genes and will be excluded.",gpr(model)[r],r))
						print(rxnGenes[!(rxnGenes %in% geneExpr$GeneID)])
					}
				}
			}
		}
		if(sum(v_selected_rxns)>0){
			v_tmp=gene2Rule(model,geneExpr,v_selected_rxns)
			v_tmp=v_tmp[!is.na(v_tmp[,"expr_val"]),,drop=FALSE]
			if(length(v_tmp[,1])>0)
				gs=rbind(gs,cbind(rxn_id=v_tmp[,"rxn_id"],cond_id=rep(cnd,length(v_tmp[,"rxn_id"])),
													expr_val=v_tmp[,"expr_val"],cnt=v_tmp[,"cnt"],gpr=v_tmp[,"gpr"]))
		}
	}
}#end get gpr expression
## differential expression: minExprFoldChange < max(Expr)/min(Expr)
 print(is(gs))
 gs=as.data.frame(gs)
 gs[,"expr_val"]=as.numeric(as.vector(gs[,"expr_val"]))
 de_gpr=rep(TRUE,nc)
 if(minExprFoldChange>0){
	
	mnmx=aggregate(expr_val ~ rxn_id,data=gs,
			FUN=function(x)	{	max(x)-minExprFoldChange*2*min(x)	}	);
	#maxxpr=aggregate(expr_val~rxn_id,data=gs,FUN=max)
	de_gpr[mnmx[mnmx$expr_val<0,"rxn_id"]]=FALSE
 }
 gpr_rxn_ind=gpr_rxn_ind & de_gpr;
 if(verboseMode>3) print(sprintf("Number of rxns after considering differential expressions: %d",sum(gpr_rxn_ind)))
##-------------------------------------------------------------------------------##
 
 #3-
 ### --------------------------------Calculate expected (fitted) flux function ------ ###
FVA_CF <- function(p,a){#FVA cost function
# p: contains line parameters m,c (line eqn: mx+c)
#a: contains 3 columns (x,range): value(expression),min (FVA min),max (FVA max)
# cost function =0 when value is in range and square of distance o.w.
# weights can be used to penalize cases when line is outside from the min-side. 
# Reversibility:if(x[1]
#

########### adjust for reversibility ###############
	 a1=apply(a[,c(2,3)],1,function(x){
		 if(x[1] <0 && x[2]<0){fmn=abs(x[1]);fmx=abs(x[2]);}# both -ve
		 else if(sign(x[1]*x[2])!=1){fmn=0; fmx=max(abs(x[1]),abs(x[2]))} # diff signs (inc 0 & others)
		 else {fmn=x[1];fmx=x[2]}# both +ve
		# print(sprintf("mn=%f,mx=%f,fmn=%f,fmx=%f",x[1],x[2],fmn,fmx))
		 return (cbind(fmn,fmx))
		 })
	a[,c(2,3)]=t(a1)
############################ ###########################

  pnlty=1; # if one no penalty for 
	cf=apply(a,1,function(x) {
                     # function value	                
	                y=ifelse((p[1]==0),p[2],
	                        ifelse(x[1]>(-p[2]/p[1]), p[1]*x[1]+p[2],0))
	                ifelse(y<x[2],pnlty*(x[2]-y)^2,ifelse(y>x[3],(y-x[3])^2,0))
	                })
	return(sum(cf))
}

 f <- function(x,p) { sapply(x,function(x) max(0,p[1]*x+p[2]))}
 
#-----------------------
all_xpc=NULL;
all_lm=NULL;
#--------------------

 #calc_expected_flx= function(gs,all_xpc){
# for each rxn find the line that best fit flux range and corresponding expression levels
#if not differentially expressed then skip.
#gs: represents expression level for each gpr in all conditions

	 #gs=gs_ALL[gs_ALL[,"MDL"]==mdl,] One model
	 erxns=as.vector(unique(gs[,"rxn_id"]))# fit fluxes for a set of rxns (one gene)
	 print(erxns)
	 gpr_rxn_ind[!(react_id(model) %in% (erxns))]=FALSE;
	 print(sprintf("Number of rxns to be fitted:%d ",sum(gpr_rxn_ind)));
	 if(sum(gpr_rxn_ind)==0){
		print("No rxns satisfying the conditions of flux variabilty and gene differentiability among the set of selected rxns.");
		return(NULL);
	 }
	 for(rxn in react_id(model)[gpr_rxn_ind]){
			xpr=gs[gs[,"rxn_id"]==rxn,c("rxn_id","cond_id","expr_val")]
			print(xpr)
			flxr=all_fva[all_fva[,"rxn_id"]==rxn,]
			x=merge(xpr, flxr, by.x = c("cond_id","rxn_id"), by.y = c("cond_id","rxn_id"))# , all = TRUE, FULL outer join
			#print(is(x));
						
			#NORmalize FVA!!!!: using objective value not important
			#FVA_CF: takes care of -ve values of fv_min, but estimate should also!
			
			#optimize function, start solution
			#rdata=x[,c("expr_val","fv_min","fv_max")]
			#print(rdata)
			
			rdata=cbind(expr_val=as.numeric(as.vector(x[,"expr_val"])),fv_min=as.numeric(as.vector(x[,"fv_min"])),
				fv_max=as.numeric(as.vector(x[,"fv_max"])) )
			print("rdata")
			print(rdata)
			#30/10/2013: choose a start solution using lm
			lnn=lm(rdata[,2]~rdata[,1])#use fv_min points, consider the sign!!
			 
			#lnx=lm(rdata[,3]~rdata[,1])
			#mdpnts=(rdata[,2]+rdata[,3])/2
			#lnd=lm(mdpnts~rdata[,1])
			#if(!is.na(lnn$coefficients[2])){ cfn=FVA_CF(p=c(lnn$coefficients[2],lnn$coefficients[1]),rdata);
			#}else{cfn=1e6}
			##cfx=FVA_CF(p=c(lnx$coefficients[2],lnx$coefficients[1]),rdata)
			#if(!is.na(lnd$coefficients[2])){cfd=FVA_CF(p=c(lnd$coefficients[2],lnd$coefficients[1]),rdata)
			#			}else{cfd=1e6}
						
			#if(cfd<cfn){ iln=c(lnx$coefficients[2],lnx$coefficients[1]);
			#}else{ iln=c(lnm$coefficients[2],lnm$coefficients[1])}
			#if(cfd<cfn){iln=c(lnd$coefficients[2],lnd$coefficients[1])
			#}else{iln=c(lnn$coefficients[2],lnn$coefficients[1])}
			iln=c(lnn$coefficients[2],lnn$coefficients[1])
			if(is.na(iln[2])|| is.na(iln[1])) iln=c(1,-1);
			print(iln)
			#iln=c(1,-1)
			ln=nlm(FVA_CF,iln,a=rdata,ndigit=15,gradtol=1e-10,steptol=1e-10)
			#ln=nlm(FVA_CF,c(1,-1),a=rdata,ndigit=15,gradtol=1e-10,steptol=1e-10)
			p=ln$estimate
			#print(p)

			x=cbind(x,xpc=f(rdata[,"expr_val"],p));
			#print(x)
			all_xpc=rbind(all_xpc,x);
			all_lm=rbind(all_lm,cbind(rxn,gpr=gpr(model)[react_id(model)==rxn],minExpr=min(rdata[,"expr_val"]),maxExpr=max(rdata[,"expr_val"]),
				minimalFlx=min(rdata[,"fv_min"]),maximalFlx=max(rdata[,"fv_max"]),slope=p[1],const=p[2]))
			
			# if(rxn>10) break;
	    }
 
if(verboseMode>2){
	print(format(Sys.time(), "write linearly fitted values to file: all_xpc_%Y%m%d_%H%M"));
	write.csv(file=format(Sys.time(), "all_xpc_%Y%m%d_%H%M.csv"),all_xpc);
	write.csv(file=format(Sys.time(), "all_lm_%Y%m%d_%H%M.csv"),all_lm);
}
 
#-------------------------------------------------------------------------# 
#4-
#N.B: reversible rxns should be able to go in backward with the same value!!
# may need a constraint dp+dn=|wt|, dp+x>=0 & dn-x>=0
all_iflx=NULL;
cnds=unique(as.vector(all_xpc[,"cond_id"]));#consider only conditions with fluxes calculated

for(cnd in cnds){
	#get all xpc values as wtflx
	#run MDC
	print("Condition No:",cnd);
    # apply initial conditions to model
     model=init_model(cnd);
    # get xpc fluxes for a set of rxns 
	 xpc=  all_xpc[all_xpc[,"cond_id"]==cnd,,drop=FALSE]
     r_ind=match(xpc[,"rxn_id"],react_id(model))
	 fvmin=fvmax=rep(NA,nc)
	 fvmin[r_ind]=as.vector(xpc[,"fv_min"]);fvmax[r_ind]=as.vector(xpc[,"fv_max"])
      
	 wtflx=numeric(nc);
	 wtflx=rep(NA,nc)
	 #wtflx[r_ind]=ifelse(abs(as.numeric(as.vector(xpc[,"xpc"])))>200,NA,xpc[,"xpc"])#>200 WILL BE NULL IN FILE
	 wtflx[r_ind]=xpc[,"xpc"]
	 for(r in r_ind){#map expected from linear fit to wildtype flux, consider sign of fva
		if( fvmax[r]<=0 | (wtflx[r] > fvmax[r] & wtflx[r] <= -fvmin[r]) ){
			wtflx[r] = -wtflx[r]
		}
	 }
	 xpr=rep(NA,nc)
	 xpr[r_ind]=as.vector(xpc[,"expr_val"]);
	  # run !!
	 print(c("Count of wtflx",sum(!is.na(wtflx))));
	 #31/10/2013: ignore biomass objective
	 iflx=findMDCFlux(model,wtflx,pct_objective=pct_objective,solver=solver)#,objVal=pct_objective
	 #print(iflx$mdcflx[gpr_rxn_ind])
	 	 all_iflx=rbind(all_iflx,cbind(cnd=rep(cnd,nc),rxn_id=react_id(model),lb=lowbnd(model),ub=uppbnd(model),gpr=gpr(model),
		expr_val=xpr,selected=gpr_rxn_ind,de_gpr,blocked_rxns,fvmin,fvmax,xpc=wtflx,iflx=iflx$mdcflx))#,fv_min=,fv_max
	}   
   # save Fluxes   

 if(verboseMode>3) write.csv(file=sprintf("all_iFlx_%s.csv",format(Sys.time(), "%Y%m%d_%H%M")),all_iflx);
  
 
#--------------Output-------------------------#

        optsol <- list(all_xpc,
	                    all_iflx,
	                    all_fva,
						all_lm,
						gs
#						stat = lp_stat,
	                  )
	    
    return(optsol)
}

 #rxnStat=cbind(all_iflx[gpr_rxn_ind,]
#	xflx=findMDCFlux(model,wtflux=wtflx,objVal=NA,solver=solver#"cplexAPI" 
#		, solverParm=data.frame(CPX_PARAM_EPRHS=1e-9))#,verboseMode = 3
