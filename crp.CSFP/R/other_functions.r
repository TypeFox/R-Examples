crp.round <- function(a){ifelse(a%%0.5==0&a%%1!=0&trunc(a)%%2==0,round(a,0)+1,round(a,0))}

crp.init <- function(PATH.IN, PATH.OUT, PORT.NAME = "portfolio.csv",
                              RISK.NAME = "rating_pd.csv", PDVAR.NAME = "pd_sector_var.csv", SEC.VAR.EST = 5,
                              LOSS.UNIT = 1e+06, NITER.MAX = 0.9999, NITER.MAX.GLOBAL = 1e+05,
                              ALPHA = c(0.999), PLOT.PDF = TRUE, CALC.RISK.CONT = FALSE, PLOT.SCALE = 1e+06,
                              PLOT.RANGE.X = c(0, 0), PLOT.RANGE.Y = c(0, 0), save.memory = FALSE,
                              file.format = "csv", portfolio=data.frame(), risk.matrix=data.frame(),
                              sec.var=data.frame()){
  cat("crp.init() is depricated. Please use init() instead.")
}

init <- function(path.in="",path.out="",port.name="portfolio.csv",
                 rating.scale.name="rating_pd.csv",sec.var.name="pd_sector_var.csv",
                 sec.var.est=5,loss.unit=1e6,Niter.max=0,alpha.max=0.9999,
                 Niter.max.global=1e5,alpha=c(0.999),PLOT.PDF=TRUE,
                 export.to.file=FALSE,calc.rc=FALSE,PLOT.scale=1e6,
                 PLOT.range.x=c(0,0),PLOT.range.y=c(0,0),save.memory=FALSE,
                 file.format="csv",portfolio=data.frame(),rating.scale=data.frame(),
                 sec.var=data.frame()){
 
  packageStartupMessage("    CreditRisk+ portfolio model \n    Copyright (C) 2011  Dr. Matthias Fischer, Kevin Jakob & Stefan Kolb

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
    Boston, MA  02110-1301, USA.
                        
    Please note that this package will not be updated anymore.
    Instead we recommend to use the GCPM package which includes
    the functionality of this package as well as more flexible
    and more powerfull extensions.\n")
  missing.path.in=FALSE
  if(missing(path.in) && (nrow(portfolio)==0 || nrow(rating.scale)==0 || (nrow(sec.var)==0 && sec.var.est==5))){                                                                               
    cat("ERROR: You have to define an input path or pass all input data directly to 'init()'.\n")
    return()
  }
  else if(missing(path.in)){
    path.in=""
    missing.path.in=TRUE
  }
  if(!missing.path.in){
    if(substr(path.in,nchar(path.in),nchar(path.in))!="/" && substr(path.in,nchar(path.in),nchar(path.in))!="\\")
      path.in=paste(path.in,"/",sep="")
  }
  if(!missing(path.out)){
    export.to.file=TRUE
    if(substr(path.out,nchar(path.out),nchar(path.out))!="/" && substr(path.out,nchar(path.out),nchar(path.out))!="\\")
      path.out=paste(path.out,"/",sep="")
  }
  if(missing(path.out) && !missing.path.in)
    path.out=path.in
  else if(missing(path.out) && export.to.file){
    cat("ERROR: Please specify path.out.\n")
    return()
  }
  else if(missing(path.out) && missing(path.in))
    path.in=""
  if(!missing(file.format)){
    if(!(file.format=="csv" || file.format=="csv2")){
      cat("Wrong specification of file.format. Please choose between csv (sep= , dec= . ) and csv2 (sep= ; dec= , ).\n")
      file.format="csv"
    }
  }
  if(alpha.max<0 || alpha.max>1)
    stop("alpha.max is not between 0 and 1!")
  if(round(Niter.max)!=Niter.max || Niter.max<0)
    stop("Niter.max has to be a non-negative integer!\nTo calculate the cdf up to a specific level use alpha.max instead.")
  if(loss.unit<=0)
    stop("Loss.unit has to be positive!\n")
  
  return(new("crp.CSFP",path.in=path.in,path.out=path.out,port.name=port.name,rating.scale.name=rating.scale.name, sec.var.name=sec.var.name,sec.var.est=sec.var.est,loss.unit=loss.unit,Niter.max=Niter.max,alpha.max=alpha.max,Niter.max.global=Niter.max.global,alpha=alpha,PLOT.PDF=PLOT.PDF,export.to.file=export.to.file,calc.rc=calc.rc,PLOT.scale=PLOT.scale, PLOT.range.x=PLOT.range.x,PLOT.range.y=PLOT.range.y,save.memory=save.memory,file.format=file.format,input=list(portfolio=portfolio,rating.scale=rating.scale,sec.var=sec.var)))

}

fo <- function(x){                                                                                   # function formatting the output of big numbers
  
  if(is.na(x))
    return("NA")
  if(!is.numeric(x))
    return(x)
  s=""
  s1=""
  for(i in 1:length(x)){
    if(abs(x[i])>=1e12)
      s1=paste(round(x[i]/1e12,2),"Tril.")
    else if(abs(x[i])>=1e9)
      s1=paste(round(x[i]/1e9,2),"Bil.")
    else if(abs(x[i])>=1e6)
      s1=paste(round(x[i]/1e6,2),"Mio.")
    else if(abs(x[i])>=1e3)
      s1=paste(round(x[i]/1e3,2),"Thd.")
    else
      s1=paste(round(x[i],2))
    if(i==1)
      s=s1
    else
      s=c(s,s1)
  }
  return(s)
}

