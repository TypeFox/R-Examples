#' Find attrition stats for TAPS waves
#'
#' Find response attrition statistics for TAPS waves
#'
#' @param outcome A character vector of the names of outcome variables of interest
#'
#' @return A list of attrition rates in population margins. The number in each category corresponds to the percentage of people within that category that left/attrited the panel
#' @docType methods
#' @author David G. Carlson \email{carlson.david@@wustl.edu} Michelle Torres: \email{smtorres@@wustl} Taeyong Park \email{t.park@@wustl.edu}
#' @seealso \code{\link{weightTAPSPACK}} \code{\link{weightTAPS}} \code{\link{variablesTAPS}} \code{\link{weightTAPSoutput}} \code{\link{simpleWeight}} \code{\link{subsetTAPS}} \code{\link{multipleImp}} \code{\link{hotdeckImp}} \code{\link{wavesTAPS}}
#' @rdname attritTAPS
#' @aliases attritTAPS,ANY-method
#' @export
setGeneric(name="attritTAPS",
           def=function(outcome)
           {standardGeneric("attritTAPS")}
)

setMethod(f="attritTAPS",
          definition=function(outcome){
    
            baseline<-subsetTAPS(outcome[1])
            f.agegend.names<-names(table(baseline$agegend))
            f.agegend.vec<-as.numeric(table(baseline$agegend))
            names(f.agegend.vec)<-f.agegend.names
            f.ethm.names<-names(table(baseline$ethm))
            f.ethm.vec<-as.numeric(table(baseline$ethm))
            names(f.ethm.vec)<-f.ethm.names
            f.educat.names<-names(table(baseline$educat))
            f.educat.vec<-as.numeric(table(baseline$educat))
            names(f.educat.vec)<-f.educat.names
            f.regmetro.names<-names(table(baseline$regmetro))
            f.regmetro.vec<-as.numeric(table(baseline$regmetro))
            names(f.regmetro.vec)<-f.regmetro.names
            f.incomcat.names<-names(table(baseline$incomcat))
            f.incomcat.vec<-as.numeric(table(baseline$incomcat))
            names(f.incomcat.vec)<-f.incomcat.names
            f.ppnet.names<-names(table(baseline$ppnet))
            f.ppnet.vec<-as.numeric(table(baseline$ppnet))
            names(f.ppnet.vec)<-f.ppnet.names
            
            final<-subsetTAPS(outcome)
            f.agegend.names<-names(table(final$agegend))
            f.agegend.vec.final<-as.numeric(table(final$agegend))
            names(f.agegend.vec.final)<-f.agegend.names
            f.ethm.names<-names(table(final$ethm))
            f.ethm.vec.final<-as.numeric(table(final$ethm))
            names(f.ethm.vec.final)<-f.ethm.names
            f.educat.names<-names(table(final$educat))
            f.educat.vec.final<-as.numeric(table(final$educat))
            names(f.educat.vec.final)<-f.educat.names
            f.regmetro.names<-names(table(final$regmetro))
            f.regmetro.vec.final<-as.numeric(table(final$regmetro))
            names(f.regmetro.vec.final)<-f.regmetro.names
            f.incomcat.names<-names(table(final$incomcat))
            f.incomcat.vec.final<-as.numeric(table(final$incomcat))
            names(f.incomcat.vec.final)<-f.incomcat.names
            f.ppnet.names<-names(table(final$ppnet))
            f.ppnet.vec.final<-as.numeric(table(final$ppnet))
            names(f.ppnet.vec.final)<-f.ppnet.names
            
            t.agegend<-round((f.agegend.vec-f.agegend.vec.final)/f.agegend.vec,3)
            t.ethm<-round((f.ethm.vec-f.ethm.vec.final)/f.ethm.vec,3)
            t.educat<-round((f.educat.vec-f.educat.vec.final)/f.educat.vec,3)
            t.regmetro<-round((f.regmetro.vec-f.regmetro.vec.final)/f.regmetro.vec,3)
            t.incomcat<-round((f.incomcat.vec-f.incomcat.vec.final)/f.incomcat.vec,3)
            t.ppnet<-round((f.ppnet.vec-f.ppnet.vec.final)/f.ppnet.vec,3)
            
            socio.diffs<-list(t.agegend, t.ethm, t.educat, t.regmetro, t.incomcat, t.ppnet)
            names(socio.diffs)<-c("AgeGender", "Ethnicity", "Schooling", "Region", "Income", "InternetUsers")
            
            return(socio.diffs)
          }
)
