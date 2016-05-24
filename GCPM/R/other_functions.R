gcpm.round <- function(a){ifelse(a%%0.5==0&a%%1!=0&trunc(a)%%2==0,round(a,0)+1,round(a,0))}

init <- function(model.type="CRP",link.function="CRP",N,seed,loss.unit,alpha.max=0.9999,
                 loss.thr=Inf,sec.var,random.numbers=matrix(),LHR, max.entries=1e3){
 
  packageStartupMessage("    Generalized Credit Portfolio Model \n    Copyright (C) 2015 Kevin Jakob & Dr. Matthias Fischer

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor,
    Boston, MA  02110-1301, USA.\n")
  
  if(alpha.max<0 || alpha.max>1)
    stop("alpha.max is not between 0 and 1!")
  if(missing(loss.unit)){
    if(model.type=="CRP")
      stop("Please provide a suitable loss unit > 0.")
    else if(model.type=="simulative")
      loss.unit=1
  }
  if(loss.unit<=0)
    stop("Loss.unit has to be positive!\n")
  
  if(!any(model.type==c("CRP","simulative")))
    stop("Wrong specification of model.type! Choose between CRP or simulativ.")
  if(!is.matrix(random.numbers))
    stop("random.numbers has to be a matrix")
  
  if(model.type=="simulative" && empty.matrix(random.numbers)){
    stop("If model.type=simulative, random numbers have to be provided.")    
  }
  if(model.type=="simulative" && missing(LHR)){
    warning("No LHR provided for simulative model, assuming equally likelihood for all szenarios.")
    LHR=rep(1,nrow(random.numbers))
  }
  if(model.type=="CRP" && missing(LHR))
    LHR=0
  if(missing(seed))
    seed=NA
  else
    seed=as.integer(round(seed))
  if(model.type=="simulative" && missing(N))
    N=nrow(random.numbers)
  else if(missing(N))
    N=0
  if(model.type=="simulative" && N<=0)
    stop("Number of simulations N has to be positive!")
  if(N>nrow(random.numbers) && model.type=="simulative"){
    warning("N is greater than the number of provided scenarios in random.numbers. Scenarios will be recycled.")
    scenarios=rep(1:nrow(random.numbers),floor(N/nrow(random.numbers)))
    if(length(scenarios)<N)
      scenarios=c(scenarios,1:(N-length(scenarios)))
  }
  else if(model.type=="simulative"){
    scenarios=1:N
  }
  else
    scenarios=-1
  if(model.type=="simulative" && loss.thr==0)
    warning("Calculating risk contributions while saving all loss szearios leads to a large memory demand. Increasing loss.thr may be appropriate.")
  if(loss.thr==Inf && model.type=="simulative")
    warning("loss.thr is not finite. Risk contributions (to EC, VaR and ES) will be not available.")
  if(!any(link.function==c("CRP","CM")))
    stop("Wrong specification of link.function. Choose between CRP and CM.")
  if(link.function=="CM" && model.type=="CRP")
    warning("If model.type=CRP only CRP link function is available.")
  if(model.type=="CRP" && missing(sec.var))
    stop("No sector variances provided.")
  else if(missing(sec.var))
    sec.var=0
  if(model.type=="CRP" && is.null(names(sec.var)))
    stop("No sector names given as names of sec.var")
  
  return(new("GCPM",N=floor(N),link.function=link.function,loss.thr=loss.thr,seed=as.numeric(seed),
             model.type=model.type,loss.unit=loss.unit,
             alpha.max=alpha.max,sec.var=sec.var,random.numbers=random.numbers,
             LHR=LHR,scenarios=scenarios,max.entries=max(0,floor(max.entries))))

}

fo <- function(x){                                                                                   # function formatting the output of big numbers
  
  s=""
  s1=""
  for(i in 1:length(x)){
    if(is.na(x[i]))
      s1="NA"
    else if(!is.numeric(x[i]))
      s1="NA"
    else if(!is.finite(x[i]))
      s1="Inf"
    
    else if(abs(x[i])>=1e12)
      s1=paste(round(x[i]/1e12,2),"T")
    else if(abs(x[i])>=1e9)
      s1=paste(round(x[i]/1e9,2),"B")
    else if(abs(x[i])>=1e6)
      s1=paste(round(x[i]/1e6,2),"M")
    else if(abs(x[i])>=1e3)
      s1=paste(round(x[i]/1e3,2),"K")
    else
      s1=paste(round(x[i],2))
    if(i==1)
      s=s1
    else
      s=c(s,s1)
  }
  return(s)
}

empty.matrix<-function(M){
  if(!is.matrix(M))
    stop("M is not a matrix!")
  return(nrow(M)==1 && ncol(M)==1 && is.na(M[1,1]))
}

num.geq<-function(x,y,tol=1e-10){
  return((x>y|abs(x-y)<tol))
}
