lin.convergence<-function(origin, aim, convergence_time_steps=NULL,start_year=NULL, end_year=NULL, before="stable", after="stable") {

  if(!is.magpie(origin)){stop("origin is no magpie object")}
  if(!is.magpie(aim)){
    if(is.numeric(aim)){
      aim<-as.magpie(array(aim,dim=dim(origin),dimnames=dimnames(origin)))
    } else {stop("aim is no magpie object")}}
  if (all(dimnames(aim)[[1]]!=dimnames(origin)[[1]])) stop("regions have to be the same")
  if (dim(origin)[3]!=1) {
    if (identical(dimnames(origin)[[3]], dimnames(aim)[[3]]) == FALSE) {
      stop("If there ist more than one name-column, dimnames have to be the same")    
    }
  }
  
  if(dim(aim)[2]==1) {
    if (is.null(convergence_time_steps)&is.null(end_year)) {
        end_year<-getYears(aim)
        end_year_num<-getYears(aim,as.integer=TRUE)
    }
    aim_new<-new.magpie(dimnames(origin)[[1]],dimnames(origin)[[2]],dimnames(origin)[[3]])
    for (name_x in 1:dim(aim)[3]) {
      aim_new[,,name_x]<-aim[,,name_x]
    }
    aim<-aim_new
    rm(aim_new)
  
  } 
  if (any(dimnames(aim)[[2]]!=dimnames(origin)[[2]])) stop("Objects need the same timesteps, or aim has to have only one timestep")

  if (is.null(start_year)) {start_year<-getYears(aim)[1]}
  if (is.null(end_year)) {
    if(is.null(convergence_time_steps)) { 
      end_year <- getYears(aim)[length(getYears(aim))]
    } else {
      end_year <- getYears(aim)[which(getYears(aim)==start_year)+ convergence_time_steps - 1]
    }
  } else {
     if(!is.null(convergence_time_steps)) {stop("cannot use both convergence_time_steps and end_year")}
  }

  if(isYear(end_year,with_y=TRUE)){end_year_num<-as.numeric(substr(end_year,2,5))}else{stop("wrong year format for convergence aim")}
  if(isYear(start_year,with_y=TRUE)){start_year_num<-as.numeric(substr(start_year,2,5))}else {stop("wrong year format for convergence aim")}  
  convergence_distance<-end_year_num-start_year_num
  dimnames(aim)[[3]]<-dimnames(origin)[[3]]
  if (isYear(before)) {
    before_num<-as.numeric(substr(before,2,5))
    convergence_distance_back<-start_year_num-before_num
  }
  
  converged<-origin 
  
  for (name_x in getNames(converged)) {
    for (year_x in getYears(converged)) {
      year_x_num<-as.numeric(substr(year_x,2,5))
      mix_up    <-  (year_x_num - start_year_num)/(convergence_distance)
      mix_down  <-  1-mix_up    
      if ((after=="stable")&(mix_up>1)) {
        mix_up<-1
        mix_down<-0
      }
      if (before=="stable"){
        if(mix_up<0) {
          mix_up<-0
          mix_down<-1
        }
      } else if (isYear(before)) {
        if(mix_up<0) {
          mix_up    <-  (start_year_num - year_x_num )/(convergence_distance_back)
          mix_down  <-  1-mix_up 
        }
      }
      converged[,year_x,name_x]<-aim[,year_x,name_x]*mix_up + origin[,year_x,name_x]*mix_down
    }
  }
  
  return(converged)
}