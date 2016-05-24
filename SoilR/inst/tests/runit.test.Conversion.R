#
# vim:set ff=unix expandtab ts=2 sw=2:
test.Conversion=function(){
   require(RUnit)
   t_start=0
   t_end=2
   tn=100
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn
   t=seq(t_start,t_end,timestep)
   th=5730
   A=new("BoundLinDecompOp",t_start,Inf,function(t){matrix(
     nrow=1,
     ncol=1,
     c(
        -log(2)/th
     )
   )})
   c0s=c(1)
   f0_Delta14C=ConstFc(1,format="Delta14C")
   f0_AFM=AbsoluteFractionModern(f0_Delta14C)
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=1,
     ncol=1,
     c(
        1
     )
   ))})
   f_D14C=function(t){0.5*t}
   Fc_D14C=BoundFc(f_D14C,t_start,t_end,format="Delta14C")
   Fc_AFM=AbsoluteFractionModern(Fc_D14C)
   k=log(0.5)/th
   
   mod_D14C=GeneralModel_14(
    t 		  		=t,
    A	    			=A,
    ivList			=c0s,
    initialValF	=f0_Delta14C,
    inputFluxes	=inputrates,
    inputFc			=Fc_D14C,
    di					=k,
    solverfunc	=deSolve.lsoda.wrapper
   )
   mod_AFM=GeneralModel_14(
   t 		  		  =t,
   A	    			=A,
   ivList		  	=c0s,
   initialValF	=f0_AFM,
   inputFluxes	=inputrates,
   inputFc			=Fc_AFM,
   di				  	=k,
   solverfunc	  =deSolve.lsoda.wrapper
   )
   C14_D14C=getC14(mod_D14C) 
   C14_AFM=getC14(mod_D14C) 
   F14_AFM=getF14(mod_AFM) 
   F14_D14C=getF14(mod_D14C) 
   Yode=getC(mod_D14C) 
   Rode=getReleaseFlux(mod_D14C) 
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.test.Conversion.pdf",paper="a4")
   m=matrix(c(1,2),2,1,byrow=TRUE)
   layout(m)
   ###C14###
   plot(t,C14_AFM[,1],type="l",lty=lt1,col=1,ylab="14C-Concentrations",xlab="Time",ylim=c(min(C14_AFM),max(C14_AFM)))
   lines(t,C14_D14C[,1],type="l",lty=lt2,col=1)
   ###F14###
   plot(t,F14_AFM[,1],type="l",lty=lt1,col=1,ylab="14C-C ratio ",xlab="Time",ylim=c(min(F14_AFM,F14_D14C),max(F14_AFM,F14_D14C)))
   lines(t,F14_D14C[,1],type="l",lty=lt2,col=1)
   legend(
   "topright",
     c(
     "anylytic sol for pool 1",
     "numeric sol for pool 1"
     ),
     lty=c(lt1,lt2),
     col=c(1,1)
   )
   dev.off()
# end plots 
# begin checks 
   tol=.02*max(C14_AFM)/tn
   checkEquals(
    C14_AFM,
    C14_D14C,
    "compare solution for 14C-Content computed with atmospheric C14 fractopm input data against the same Problem with input in AbsoluteFractionModern",
    tolerance = tol,
   )
   checkEquals(
    F14_AFM,
    F14_D14C,
    "compare solution for F14 fraction computed from atmospheric C14 fraction input data against the same Problem with input in AbsoluteFractionModern",
    tolerance = tol,
   )

 }
