findwtdinteraction <- function(reg, across, by=NULL, at=NULL, acrosslevs=NULL, bylevs=NULL,
atlevs=NULL, weight=NULL, dvname=NULL, acclevnames=NULL, bylevnames=NULL,
atlevnames=NULL, stdzacross=FALSE, stdzby=FALSE, stdzat=FALSE, limitlevs=20,
                               type="response"){
    UseMethod("findwtdinteraction")
}
