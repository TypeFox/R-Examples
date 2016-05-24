"add.lmomco.axis" <-
function(side=1, twoside=FALSE,
         side.type=c("NPP", "RI", "SNV"), otherside.type=c("NA", "RI", "SNV", "NPP"),
         alt.lab=NA,
         NPP.control=NULL, RI.control=NULL, SNV.control=NULL, ...) {

   other.side <- switch(as.character(side), "1"=3, "2"=4, "3"=1, "4"=1)

   side.type      <- match.arg(side.type)
   otherside.type <- match.arg(otherside.type)
   if(otherside.type == "NA") otherside.type <- NA
   if(twoside & ! is.na(otherside.type)) twoside <- FALSE

   "my.nonexceeds" <- function(minors=FALSE) {
      if(minors) {
         F <- c(0.55, 0.65, 0.75, 0.825, 0.850, 0.875,
                0.91, 0.92, 0.93, 0.94, 0.96, 0.97,
                0.9925, 0.996, 0.997)
         return(c(sort(1-F), F))
      } else {
         F <- c(0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999)
         return(c(sort(1-F), 0.5, F))
      }
   }

   if(is.null(NPP.control)) {
      the.label <- ifelse(is.na(alt.lab), "NONEXCEEDANCE PROBABILITY", alt.lab)
      NPP.control <- list(label=the.label,
                          probs=my.nonexceeds(minors=TRUE),
                          probs.label=my.nonexceeds(minors=FALSE),
                          digits=3, line=3, as.exceed=FALSE)
   }
   if(is.null(RI.control)) {
      the.label <- ifelse(is.na(alt.lab), "RECURRENCE INTERVAL, YEARS", alt.lab)
      RI.control <- list(label=the.label,
                         Tyear=c(2, 5, 10, 25, 50, 100, 200, 500), line=2)
   }
   if(is.null(SNV.control)) {
      the.label <- ifelse(is.na(alt.lab), "STANDARD NORMAL VARIATE", alt.lab)
      SNV.control <- list(label=the.label,
                          begin=-5, end=5, by=0.5, line=2)
   }

   NPPf <- function(side, other.side) {
      dots <- list(...)
      tcl <- ifelse("tcl" %in% names(dots), dots$tcl, par()$tcl)
      NPP <- NPP.control$probs;  NPP.lab <- NPP.control$probs.lab
      if(NPP.control$as.exceed) {
         the.true.NPP.lab <- 1 - NPP.lab
      } else {
         the.true.NPP.lab <-     NPP.lab
      }
      qNPP <- qnorm(NPP); qNPP.lab <- qnorm(NPP.lab)
      NPP.lab <- format(the.true.NPP.lab, nsmall=NPP.control$digits)
      # By placing tcl last (after ...), its value will trump that potentially in ...
      Axis(qNPP,     at=qNPP,     labels=NA,      side=side, ..., tcl=0.8*tcl)
      Axis(qNPP.lab, at=qNPP.lab, labels=NPP.lab, side=side, ..., tcl=1.3*tcl)
      mtext(NPP.control$label, line=NPP.control$line, side=side)
      if(twoside) {
         Axis(qNPP,     at=qNPP,     labels=NA,      side=other.side, ..., tcl=0.8*tcl)
         Axis(qNPP.lab, at=qNPP.lab, labels=NPP.lab, side=other.side, ..., tcl=1.3*tcl)
      }
   }

   RIf <- function(side, other.side) {
      F <- 1 - 1/RI.control$Tyear; qF <- qnorm(F); labF <- RI.control$Tyear
      Axis(qF, at=qF, labels=labF, side=side, ...)
      if(twoside) {
         Axis(at=qF,  labels=NA, side=other.side, ...)
      }
      mtext(RI.control$label, line=RI.control$line, side=side)
   }

   SNVf <- function(side, other.side) {
      SNV <- NULL
      try( SNV <- seq(SNV.control$begin, SNV.control$end, by=SNV.control$by) )
      if(is.null(SNV)) {
         warning("Poorly constructed SNV.control, trapping, and using alternative")
         SNV <- seq(-5, 5, by=0.5)
      }
      Axis(SNV, at=SNV, side=side, ...)
      mtext(SNV.control$label, line=SNV.control$line, side=side)
      if(twoside) {
         Axis(SNV, at=SNV, side=other.side, ...)
      }
   }

   NULLf <- function() { return("no axis function made") }

   primary.axis <- switch(side.type, NPP=NPPf, RI=RIf, SNV=SNVf, NULLf)
   primary.axis(side, other.side)

   if(! is.na(otherside.type)) {
      secondary.axis <- switch(otherside.type, NPP=NPPf, RI=RIf, SNV=SNVf, NULLf)
      secondary.axis(other.side, side)
   }
}

