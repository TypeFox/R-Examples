##################################################################################
# Set Method: buildFormula
##################################################################################
buildFormula <- function(object, varlist = NULL){}
setMethod('buildFormula', signature(object = 'rasclass'),

function(object, varlist = NULL){
	
	if(is.null(object@samplename)){
		stop('please specify samplename in the rasclass object')
	}
	if(is.null(object@data)){
		stop('please load data into the rasclass object')
	}

	logitFormula = paste(object@samplename, ' ~ ', sep = '')

	# Add variables to the formula
	if(is.null(varlist)){
		for(thisName in names(object@data)){
			if(thisName != object@samplename){
				logitFormula = paste(logitFormula, thisName, ' + ', sep = '')
			}
		}
	}
	else{
		for(thisName in varlist){
			logitFormula = paste(logitFormula, thisName, ' + ', sep = '')
		}		
	}
	
	object@formula <- substr(logitFormula, 0 , nchar(logitFormula)-3)
	
	cat('Built Formula: ')
	show(object@formula)

	object
}
)