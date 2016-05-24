#' Linear gap interpolation
#' 
#' This function linearly interpolates gaps in data series, such as daily
#' temperature records.
#' 
#' The function returns a list with two elements: interp is a new vector, in
#' which all gaps in x have been linearly interpolated. missing is a second
#' vector, which contains information on which values were filled in by
#' interpolation.
#' 
#' @param x a numeric vector, or a vector that can be coerced with as.numeric.
#' Missing values are either NA or non-numeric values.
#' @return \item{interp}{ numeric vector, in which all gaps in x have been
#' linearly interpolated} \item{missing}{ boolean vector of the same length as
#' interp and x, which marks all gaps in x as TRUE}
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @keywords utility
#' @examples
#' 
#' weather<-make_all_day_table(KA_weather)
#' Tmin_int<-interpolate_gaps(KA_weather[,"Tmin"])
#' weather[,"Tmin"]<-Tmin_int$interp
#' weather[,"Tmin_interpolated"]<-Tmin_int$missing
#' 
#' Tmax_int<-interpolate_gaps(KA_weather[,"Tmax"])
#' weather[,"Tmax"]<-Tmax_int$interp
#' weather[,"Tmax_interpolated"]<-Tmax_int$missing
#' 
#' #this function is integrated into the fix_weather function, but it can also be run on its own.
#' 
#' 
#' @export interpolate_gaps
interpolate_gaps <-
function(x)   #x is the vector to work with
            {
            suppressWarnings(x<-as.numeric(as.character(x)))
            gaps<-list()
            gap_length<-c()
            fake1<-FALSE
            fakelast<-FALSE
            if(length(which(!is.na(x)))>0)
            {if(is.na(x[1])) {fake1<-TRUE; x[1]<-x[min(which(!is.na(x)))]}
            if(is.na(x[length(x)])) {fakelast<-TRUE; x[length(x)]<-x[max(which(!is.na(x)))]}
            miss<-which(is.na(x))
            if(length(miss)>0)
            {for (na in 1:length(miss))
                {if(na==1) {curr<-miss[na]
                           curr_gap<-1} else
                   {if(miss[na]==miss[na-1]+1) {curr<-c(curr,miss[na])
                                              curr_gap<-curr_gap+1} else
                                              {curr<-miss[na]
                                               curr_gap<-1}}
                   if(na==length(miss))
                      {gaps[[length(gaps)+1]]<-curr
                       gap_length<-c(gap_length,curr_gap)} else
                   if(!miss[na]==miss[na+1]-1)
                      {gaps[[length(gaps)+1]]<-curr
                       gap_length<-c(gap_length,curr_gap)}}
            for (gg in 1:length(gap_length))
                {start_val<-x[min(gaps[[gg]])-1]
                end_val<-x[max(gaps[[gg]])+1]
                valdiff<-end_val-start_val
                for (gl in 1:gap_length[gg])
                   {x[gaps[[gg]][gl]]<-start_val+(valdiff/(gap_length[gg]+1))*gl}
            }}
            if(fake1) miss<-c(1,miss)
            if(fakelast) miss<-c(miss,length(x))
            if(length(miss)==0) miss<-rep(FALSE,length(x)) else {mi<-rep(FALSE,length(x))
                                                                 mi[miss]<-TRUE
                                                                 miss<-mi}
            } else {warning("no data in dataset! Nothing to interpolate")
                    interpr<-rep(NA,length(x))
                    miss<-rep(TRUE,length(x))}
                return(list(interp=x,missing=miss))
        }
