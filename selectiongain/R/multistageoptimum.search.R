
# package selectiongain

# modified at 28-06-2013, for 1MAS+2PS
# modified at 13.08.2015 by Jose Marulanda for Preselection on Nurseries for an uncorrelated trait. 
# modified at 13.08.2015 by Jose Marulanda for a better use of the testcross seed production information
# All modified lines or added lines end with #JM

# modified at 2015-11-09, 
# 1, for the warning of "budget too small",  changed from stop to warning, make more freedom
# 2, try to use matrix multiplication to relace vactor multimplication such that the RRO3.2.2 and MKL can take advantages.

# 3, added t2 free = FALSE,  if FALSE, the cost of using T3 and T2 testers will be accounted seperately

# 

multistageoptimum.search<-function (maseff=0.4, VGCAandE, 
  VSCA=c(0,0,0,0), CostProd = c(0.5,1,1), CostTest = c(0.5,1,1), 
  Nf = 10, Budget = 10021, N2grid = c(11, 1211, 30), 
  N3grid = c(11, 211, 5), L2grid=c(1,3,1), L3grid=c(6,8,1),
  T2grid=c(1,1,1), T3grid=c(1,1,1),R2=1,R3=1,  alg = Miwa(),detail=FALSE,fig=FALSE,
  alpha.nursery=1,cost.nursery=c(0,0) #JM
  ,t2free= FALSE, # Mi 2015-11-10
  parallel.search=FALSE 
  )

{  
  if (parallel.search)
  {    
    no_cores <- detectCores() - 1
    #no_cores <- 1
    cl <- makeCluster(no_cores);
  #  clusterEvalQ(cl,library(selectiongain))
  # just using for testing grid2
    clusterEvalQ(cl,{library(selectiongain);source("multistageoptimum.grid.R")}) 
  }
  
#  source('H:/2015-Qingdao/2015-08-15-selectiongain/2_selectiongain_2.0.40_Modif_JM/selectiongainv49/R/multistageoptimum.grid.R')
  
# pre-define parameters
  Vgca=VGCAandE
  Vsca=VSCA
	
  L2limit=L2grid
  L3limit=L3grid
  T2limit=T2grid
  T3limit=T3grid
  
  alpha.nur=alpha.nursery # JM
  Cost.nur=cost.nursery # JM
  
  dim=(T3limit[2]-T3limit[1]+1)/T3limit[3]*(T2limit[2]-T2limit[1]+1)/T2limit[3]*(L3limit[2]-L3limit[1]+1)/L3limit[3]*(L2limit[2]-L2limit[1]+1)/L2limit[3]

  gainmatrix=array(0,c(1,18)) # JM  # modifyed by X. Mi v49
  colnames(gainmatrix)<-c("Nf","Nini","alpha.nur","N1","N2","N3","L2","L3","T2","T3","R2","R3","Bini","B1","B2","B3","Budget","Gain") # JM  # modifyed by X. Mi v49


# main function

    if (alpha.nur == 1 )                                           # JM
       {                                                           # JM
         CostProdMod1<-CostProd[1]+Cost.nur[1]                     # JM
	     warning("No nursery is used as alpha.nursery is set to 1. Then cost of production in Nursery added to CostProd[1]")  # JM
       }       
    else
	   {
	   CostProdMod1<-CostProd[1]
	   }

	   
  if (Budget< c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1]))
  {
    warning("budget too small, try value => c(N2grid[1]*L2grid[1]*T2grid[1]+N3grid[1]*L3grid[1]*T3grid[1])")
  }
 
   # modified at 2015-11-09, from stop to warning, make more freedom
  
  if (Budget> c(N2grid[2]*L2grid[2]*T2grid[2]+N3grid[2]*L3grid[2]*T3grid[2]))
  {
    warning("budget too great, try value =>  c(N2grid[2]*L2grid[2]*T2grid[2]+N3grid[2]*L3grid[2]*T3grid[2])")
  }
  
  # modified at 2015-11-09, from stop to warning, make more freedom

 if (length(CostTest)!= 3)
  {
    stop( "dimension of CostTest has to be 3")
  }  
  if (length(CostProd)!= 3)
  {
    stop( "dimension of CostProd has to be 3")
  }  


for (T3 in seq.int(T3limit[1],T3limit[2],T3limit[3]))
{
  for (T2 in seq.int(T2limit[1],T2limit[2],T2limit[3]))
	{
	  for (L3 in L3limit[1]:L3limit[2])
	  # seq(T3limit[1],T3limit[2],T3limit[3])
	  
		{
		   for (L2 in L2limit[1]:L2limit[2])
		   #seq(T3limit[1],T3limit[2],T3limit[3])
		   
			{
			      allocation = c(Nf,0,alpha.nur,0,0,0,L2,L3,T2,T3,R2,R3,0,0,0,0,0,0) # JM # modifyed by X. Mi v49
            gainmatrix=rbind(gainmatrix,allocation)
			} 
		}
	}
}
  
# here begins the loop circle  


theloop<-function(j,gainmatrix, N2grid= N2grid, N3grid= N3grid,maseff,t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=c(0,0) )

{

i=j+1

#alpha.nur=alpha.nursery
Cost.nur=cost.nursery 

alpha.nur=gainmatrix[i,"alpha.nur"]

#BudgetDH=gainmatrix[i,"U-DH"]
#BudgetMAS=10000-BudgetDH
L3=gainmatrix[i,"L3"]
L2=gainmatrix[i,"L2"]

T3=gainmatrix[i,"T3"]
T2=gainmatrix[i,"T2"]
R3=gainmatrix[i,"R3"]
R2=gainmatrix[i,"R2"]
N.fs=gainmatrix[i,"Nf"]

L1=1
T1=1
R1=1

if (!is.na(maseff))
{
corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3),Rep=c(R1,R2,R3),T=c(T1,T2,T3),index=FALSE,maseff=maseff)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)

# CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM
if(t2free) 
{
  CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM
}else
{
  CostProdloop=c(CostProd[1],CostProd[2],CostProd[3]) # Mi
} 

result=multistageoptimum.grid(N.upper = c(100000,N2grid[2],N3grid[2]), N.lower =  c(1,N2grid[1],N3grid[1]), 
                               Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3]), 
							   Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop, Nf = Nf, 
							   detail = FALSE, alg = Miwa(),fig=FALSE # JM
							   ,alpha.nursery = 1,cost.nursery = c(0,0)
							   )
gainmatrix[i,"Budget"]=Budget
gainmatrix[i,"Nini"]= result[1]                               # JM
gainmatrix[i,"Bini"]= result[1]*(Cost.nur[1]+Cost.nur[2])     # JM
gainmatrix[i,"N1"]= result[2]                                 # JM
gainmatrix[i,"B1"]= result[2]*(CostTest[1]+CostProdMod1)      # JM
gainmatrix[i,"N2"]= result[3]                                 # JM
gainmatrix[i,"N3"]= result[4]                                 # JM

gainmatrix[i,"B2"]= result[3]*( (L2*T2*R2*CostTest[2])+(CostProd[2]*T2)) # JM

if(t2free) 
 {
  gainmatrix[i,"B3"]= result[4]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*(T3-T2))) # JM
 }else
 {
  gainmatrix[i,"B3"]= result[4]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*T3))# Mi
 } 

gainmatrix[i,"Gain"]= result[5]                                   # JM
}else
{
corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L2,L3),Rep=c(R2,R3),T=c(T2,T3),index=FALSE,maseff=maseff)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)
#CostProdloop=c(CostProd[1]+(CostProd[2]*T2),CostProd[3]*(T3-T2))  # JM

if(t2free) 
{
  CostProdloop=c(CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM
}else
{
  CostProdloop=c(CostProd[2],CostProd[3]) # Mi
} 


result=multistageoptimum.grid( N.upper = c(N2grid[2],N3grid[2]), N.lower =  c(N2grid[1],N3grid[1]), 
                               Vg=Vgca[1],corr = corr.matrix, width =  c(N2grid[3],N3grid[3]), 
							   Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop, 
							   Nf = Nf, detail = FALSE, alg = Miwa(),fig=FALSE  # JM
							   ,alpha.nursery = alpha.nur,cost.nursery = Cost.nur
							   )
gainmatrix[i,"Budget"]=Budget
gainmatrix[i,"Nini"]= result[1]                              # JM
gainmatrix[i,"Bini"]= result[1]*(Cost.nur[1]+Cost.nur[2])    # JM
#gainmatrix[i,"N1"]= result[1]
#gainmatrix[i,"B1"]= result[1]*(CostTest[1]+CostProd[1])
gainmatrix[i,"N2"]= result[2]                                # JM
gainmatrix[i,"N3"]= result[3]                                # JM

gainmatrix[i,"B2"]= result[2]*( (L2*T2*R2*CostTest[2])+(CostProdMod1+(CostProd[2]*T2)))# JM


if(t2free) 
{
  gainmatrix[i,"B3"]= result[3]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*(T3-T2))) # JM
}else
{
  gainmatrix[i,"B3"]= result[3]*( (L3*T3*R3*CostTest[3])+(CostProd[3]*T3))# Mi
} 



gainmatrix[i,"Gain"]= result[4]                                     # JM

}

gainmatrix[i,]

}
   
if (!parallel.search)
{
  
  for (j in 1:dim )  
  {
    gainmatrix[j+1,]= theloop(j=j,gainmatrix,N2grid= N2grid, N3grid= N3grid,maseff,t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=cost.nursery)

  }
}else if(parallel.search)
{
  #clusterExport(cl, "alpha.nursery")
   resulta<- parSapply(cl=cl, 1:dim, FUN=theloop,gainmatrix,N2grid= N2grid, N3grid= N3grid,maseff=maseff,t2free=t2free,CostTest=CostTest,CostProd=CostProd,Budget=Budget,Nf=Nf,alg=alg,cost.nursery=cost.nursery)
   gainmatrix[1:dim+1,]<-t(resulta)
}  


  
  
Output=round( gainmatrix,digits=1)
Output[,"Gain"]=round( gainmatrix[,"Gain"],digits=3)
output=Output[2:c(dim+1),c("Nf", "Nini","N1", "N2", "N3","L2","L3","T2","T3","R2","R3","Bini","B1","B2","B3","Budget","Gain")]  # JM
 
gainmax=max(output[,"Gain"])
   
gainlocation = which(output[,"Gain"]==gainmax,arr.ind =TRUE)


if (fig==TRUE)
{ 
  i=gainlocation[1]
#BudgetDH=gainmatrix[i,"U-DH"]
#BudgetMAS=10000-BudgetDH
L3=gainmatrix[i,"L3"]
L2=gainmatrix[i,"L2"]

T3=gainmatrix[i,"T3"]
T2=gainmatrix[i,"T2"]
R3=gainmatrix[i,"R3"]
R2=gainmatrix[i,"R2"]
N.fs=gainmatrix[i,"Nf"]

L1=1
T1=1
R1=1

corr.longin.mas.index = multistagecor(VGCAandE=Vgca,VSCA=Vsca,L=c(L1,L2,L3),Rep=c(R1,R2,R3),T=c(T1,T2,T3),index=FALSE,maseff=maseff)

corr.matrix=corr.longin.mas.index

CostTestloop=c(CostTest[1],CostTest[2]*L2*T2*R2,CostTest[3]*L3*T3*R3)
# CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM

if(t2free) 
{
  CostProdloop=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2))  # JM
}else
{
  CostProdloop=c(CostProd[1],CostProd[2],CostProd[3]) # Mi
} 

result=multistageoptimum.grid( N.upper = c(100000,N2grid[2],N3grid[2]), N.lower =  c(1,N2grid[1],N3grid[1]),
                               Vg=Vgca[1],corr = corr.matrix, width =  c(1,N2grid[3],N3grid[3]), 
							   Budget = Budget, CostProd =CostProdloop, CostTest = CostTestloop, Nf = Nf, 
							   detail = detail, alg = Miwa(),fig=TRUE # JM
							   ,alpha.nursery=alpha.nursery,cost.nursery=cost.nursery)
  }

  if (detail==TRUE)
  {
     output
  }else  if (detail!=TRUE )
  {
     output[gainlocation[1],]
	}

}
