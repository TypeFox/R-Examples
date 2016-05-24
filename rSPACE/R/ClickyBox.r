## tclck2 package - testing for wolverine
enter.parameters<-function(Parameters=NULL){

tt <- tktoplevel()            # Create window
tkwm.title(tt,"Parameters")   # Title

done <- tclVar(0)
if(is.null(Parameters)){
  Parameters=list(
    N               = 50,                 
    lmda            = 0.933,              
    n_yrs           = 10,                 
    MFratio         = c(0.6, 0.4),  
    buffer          = c( 16, 25),         
    moveDist          = c(8.5, 12.5),       
    moveDistQ         = c(0.9, 0.7),        
    maxDistQ           = c(0.95,0.95),          
    grid_size       = 100,                      
    habitat.cutoff = 1,                  
    sample.cutoff   = 0.5,
    n_visits        = 6
    ) }                


         
         
         
tkgrid(tklabel(tt,text=" "))
         
tkgrid(tklabel(tt,text="-------Population simulation-------"))
N.Val <- tclVar(Parameters$N)
enterN <-tk2entry(tt,textvariable=N.Val, width=5)
tkgrid(tklabel(tt,text="Initial population size? "),enterN)
         
lmda.val<-tclVar(Parameters$lmda)
enterL <-tk2entry(tt,textvariable=lmda.val, width=5)
tkgrid(tklabel(tt,text="Population growth rate? "),enterL)

yrs.val<-tclVar(Parameters$n_yrs)
enter.yrs <-tk2entry(tt,textvariable=yrs.val, width=5)
tkgrid(tklabel(tt,text="Number of years? "),enter.yrs)

grps.val<-tclVar(length(Parameters$MFratio))
enterG <-tk2entry(tt,textvariable=grps.val, width=2)
tkgrid(tklabel(tt,text="Number of individual types? "),enterG)

MF.val<-tclVar(Parameters$MFratio)
enterMF<-tk2entry(tt,textvariable=MF.val, width=15)
tkgrid(tklabel(tt,text="Proportion of population by type? "),enterMF)


tkgrid(tklabel(tt,text=" "))
tkgrid(tklabel(tt,text="-----Movement parameters-----"))

buff.val<-tclVar(Parameters$buffer)
enterBuff<-tk2entry(tt,textvariable=buff.val, width=15)
tkgrid(tklabel(tt,text="Buffer distance between activity centers (km)? "),enterBuff)

how.far.val<-tclVar(Parameters$moveDist)
enterHF<-tk2entry(tt,textvariable=how.far.val, width=15)
tkgrid(tklabel(tt,text="Movement radius (km)? "),enterHF)

how.much.val<-tclVar(Parameters$moveDistQ)
enterHM<-tk2entry(tt,textvariable=how.much.val, width=15)
tkgrid(tklabel(tt,text="Proportion of movements in radius? "),enterHM)

maxDistQ.val<-tclVar(Parameters$maxDistQ)
entermaxDistQ<-tk2entry(tt,textvariable=maxDistQ.val, width=15)
tkgrid(tklabel(tt,text="Max proportion of movements to allow (1 = include all) "),entermaxDistQ)


tkgrid(tklabel(tt,text=" "))
tkgrid(tklabel(tt,text="---------------Sampling Design---------------"))

grid.val<-tclVar(Parameters$grid_size)
enter.grid <-tk2entry(tt,textvariable=grid.val, width=5)
tkgrid(tklabel(tt,text="Cell size (in km\u00b2)? "),enter.grid)

HRcutoff.val<-tclVar(Parameters$habitat.cutoff)
enterHRcut<-tk2entry(tt,textvariable=HRcutoff.val, width=5)
tkgrid(tklabel(tt,text="Minimum habitat value for activity centers? "),enterHRcut)

gridcut.val<-tclVar(Parameters$sample.cutoff)
enter.grid2 <-tk2entry(tt,textvariable=gridcut.val, width=5)
tkgrid(tklabel(tt,text="Proportion of cell in habitat? "),enter.grid2)

visit.val <- tclVar(Parameters$n_visits)
enter.visit <-tk2entry(tt,textvariable=visit.val, width=5)
tkgrid(tklabel(tt,text="Maximum visits per year? "),enter.visit)





tkgrid(tklabel(tt,text=" "))

# Subfunction #
check.values<-function(Parameters,n.grps){
  ok.check<-1
    
  if(length(Parameters$buffer) < n.grps) {ok.check <- 0}
  if(length(Parameters$maxDistQ) < n.grps) {ok.check <- 0}
  if(length(Parameters$moveDistQ) < n.grps) {ok.check <- 0}   
  if(length(Parameters$moveDist) < n.grps) {ok.check <- 0}  
  
  
  if(ok.check == 0) {cat("Something's wrong!\n");flush.console()}
  return(ok.check)
  } 
##

# Create two buttons to set the value of done
OK.but <- tkbutton(tt,text="  OK  ",    command=function() tclvalue(done)<-1)
Cancel.but <- tkbutton(tt,text="Cancel",command=function() tclvalue(done)<-2)
tkgrid(OK.but)
tkgrid(Cancel.but)

tkbind(tt,"<Destroy>",function() tclvalue(done)<-2)
tkwait.variable(done)

doneVal <- as.integer(tclvalue(done))


if(doneVal==1){
	Parameters=list(N=as.numeric(tclvalue(N.Val)),
         lmda=as.numeric(tclvalue(lmda.val)),
         n_yrs=as.numeric(tclvalue(yrs.val)),
         MFratio = as.numeric(unlist(strsplit(tclvalue(MF.val),split=' '))),
         buffer = as.numeric(unlist(strsplit(tclvalue(buff.val),split=' '))),
         moveDist = as.numeric(unlist(strsplit(tclvalue(how.far.val),split=' '))),
         moveDistQ = as.numeric(unlist(strsplit(tclvalue(how.much.val),split=' '))),
         maxDistQ = as.numeric(unlist(strsplit(tclvalue(maxDistQ.val),split=' '))),
         grid_size=as.numeric(tclvalue(grid.val)),
         habitat.cutoff=as.numeric(tclvalue(HRcutoff.val)),
         sample.cutoff=as.numeric(tclvalue(gridcut.val)),
         n_visits=as.numeric(tclvalue(visit.val)))
         
  n.grps<-as.numeric(tclvalue(grps.val))
  check.values(Parameters,n.grps)
}   

tkdestroy(tt)

return(Parameters)
#return('Hi!')
}
