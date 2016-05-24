# 1/6- Definition: to be added to file AllClasses.R

#------------------------------------------------------------------------------#
#                  definition of the class optsol_dynamicFBA                   #
#------------------------------------------------------------------------------#

setClass("optsol_dynamicFBA",
      representation(
	concentrationMatrix="numeric", # Matrix of extracellular metabolite concentrations
	excRxnNames="character",       # Names of exchange reactions for the EC metabolites
	timeVec="numeric",             # Vector of time points
	biomassVec="numeric"           # Vector of biomass values
      ),
      contains = "optsol_optimizeProb",
      package = "sybil"
)

#showClass("optsol_dynamicFBA")


#-------------------------------------------------------------------------------#

# 2/6-Constructor: to be added to file AllClasses-constructors.R
# optsol_dynamicFBAClass
optsol_dynamicFBA <- function(solver, method, nprob,
                #lpdir,
                ncols, nrows, 
                #objf,
                fld,concmat,exRxn,tmVec,bmVec) {
    if (missing(solver) || 
        missing(method) ||
        missing(nprob)  ||
        #missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        #missing(objf)   ||
        missing(fld)    ||
        missing(bmVec) ||
        missing(tmVec)
       ) {
        stop("Not enough arguments for creating an object of class optsol_dynamicFBA!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_dynamicFBA", 
        solver       = as.character(solver),
        method       = as.character(method),
        num_of_prob  = as.integer(nprob),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob),
        lp_ok        = integer(nprob),
        lp_stat      = integer(nprob),
        #lp_dir       = as.character(lpdir),
        #obj_function = as.character(objf),
        fluxdist     = fldist,
     #   num_of_steps_executed=nsteps,
	concentrationMatrix=concmat,
	excRxnNames=exRxn,
	timeVec= tmVec,  
        biomassVec= bmVec
       )
}


#-----------------------------------------------------------------------------#

# 3/6-dynamicFBA: dynamicFBA.R
##              3.1 Check
##              3.2 optimizeProb -> get OptObj
##              3.3 Set Bounds
##              3.4 Call SimpleFBA
##              3.5 Store Solution
##              3.6 Adjust OUTPUT to class optsol_dynamicFBA

#-----------------------------------------------------------------------------#

# 4/6-  accessors: to be added to file optsol_dynamicFBA-accessors.R

#-----------------------------------------------------------------------------#

# 5/6-plot: to be added to file plot-methods.R

# optsol_dynamicFBAClass
setMethod("plot", signature("optsol_dynamicFBA","missing"),
          function(x,y,
                 ylim=50,
                   xlab = "",
                   ylab = "Value",
                   type = "p",
                   pch = 20,
                   col = "black",             
    #               collower, colupper, pchupper, pchlower,
    #               dottedline = TRUE,
                  plotRxns=NULL,
                   baseline = 0,
                    ...) {
                if(missing(plotRxns)){
                   plot(x@timeVec,x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
                   }
                else {
			def.par <- par(no.readonly = TRUE);
			layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
			#layout.show(2);
			# first plot biomass
			 plot(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 1
			                    ,main='Cell density',xlab='Time(hrs)',ylab="X(g/l)",type="l",lwd=2);
			 points(x@timeVec,x@biomassVec, col = "red",lwd=2);

			# plot concentrations
			##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
			## define min/max ()plot(x@timeVec,
			ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
                        ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
			for ( i in 1:length(plotRxns) ){
				plotInd=(x@excRxnNames %in% plotRxns[i]);
				#print( x@concentrationMatrix[plotInd]);
				if(i==1){
                                   plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                                   type="l", col =i,main="Concentrations",ylab = "mmol",xlab='Time(hrs)'
                                   ,ylim=c(ymin,ymax));
                                   }
                                else{
       				  lines(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"), col =i);
       				  }
				}
			legend(1.8,ymax, plotRxns, col=1:length(plotRxns), lty=1);
                }
                  
		#if (!missing(plotRxns)){		   }
          }
)


#-----------------------------------------------------------------------------#

# 6/6-Documentation: optsol_dynamicFBA.Rd


#New version-----------------------------
