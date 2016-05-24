#---------------------------------------------
#											
#				R FUZZY TOOLKIT				 
#											
#---------------------------------------------

# Craig Knott
# Luke Hovell
# Nathan Karimian

# Based on the work started by Jonathan Garibaldi and the IMA group at the 
# University of Nottingham http://www.ima.ac.uk

#---------------------------------------------
#
# Membership functions
#
#---------------------------------------------

gaussbMF <- function(mfName, x, mfParams) {
# Inputs	: mfName(string) representing the name of the membership function, 
#			x (numeric vector) which should be the same as the variable it will be added to (if at all), and 
#			mfParams (numeric vector) representing the left_sigma, left_mean, right_sigma, right_mean, and height.
# Outputs	: Numeric vector holding the values after being applied to a guassion curve with the above inputs.

	# Validate user input.
	mfValidate(mfName, mfParams)
		
	left_sigma	<- mfParams[1]
	left_mean	<- mfParams[2]
	right_sigma	<- mfParams[3]
	right_mean	<- mfParams[4]
	height		<- mfParams[5]
		
	vec_boolA = (x <= left_mean)
	vec_boolB = (x >= right_mean)
	
	# Returns a list of the mfName, mfType, mfX, mfParams and mfVals.	
	list(
		mfName= 	mfName,
		mfType=	"gaussbmf",
		mfX=		x,
		mfParams=	mfParams,
		mfVals =	(exp(-(x-left_mean)^2/(2*left_sigma^2))*vec_boolA + (1-vec_boolA))*
				(exp(-(x-right_mean)^2/(2*right_sigma^2))*vec_boolB + (1-vec_boolB))			
	)
}

gaussMF <- function(mfName, x, mfParams)  { 
# Inputs	: mfName (String) representing the membership function name, 
#			x (numeric vector) such as 1:10 representing bounds/range and 
#			mfParams (numeric vector) of values sigma, mean, and height
# Outputs	: Vector of values after being applied to a gaussian function

	# Validate user input.
	mfValidate(mfName, mfParams)
		
	sigma 	<- mfParams[1]
	mean 	<- mfParams[2]
	height	<- mfParams[3]
	
	# Returns a list of the mfName, mfType, mfX, mfParams and mfVals.
	list(
		mfName=	mfName,
		mfType=	"gaussmf",
		mfX= 		x,
		mfParams=	mfParams,
		mfVals=	(height * exp(-(((x - mean)^2)/(2*(sigma^2)))))
	)
}

trapMF <- function(mfName, x, mfParams) {
# Inputs	: mfName (String) representing membership function name, 
#			x (numeric vector) such as 1:10 representing bounds/range and 
#			mfParams (numeric vector) of values left_foot, left_shoulder, right_shoulder, right_foot and height
# Outputs	: Numeric vector of values after being applied to a trapezoidal function
	
	# Validate user input.
	mfValidate(mfName, mfParams)
		
	left_foot		<- mfParams[1]
	left_shoulder		<- mfParams[2]
	right_shoulder	<- mfParams[3]
	right_foot		<- mfParams[4]
	height			<- mfParams[5]
		
	y = pmax(pmin((x-left_foot)/(left_shoulder-left_foot), 
		1,
		(right_foot-x)/(right_foot-right_shoulder)),0) * height

	y[is.na(y)]= 1; 
	
	# Returns a list of the mfName, mfType, mfX, mfParams and mfVals.
	list(
		mfName=	mfName,
		mfType= 	"trapmf",
		mfX= 		x,
		mfParams= 	mfParams,
		mfVals= 	y
	)
}

triMF <- function(mfName, x, mfParams) {
# Inputs	: mfName (String) representing membership function name, 
#			x (numeric vector) such as 1:10 representing the bounds/range,
#			mfParams (numeric vector) representing the left, mean, right and height values. 
# Outputs	: Vector of values after being applied to a triangular function
	
	# Validate user input.
	mfValidate(mfName, mfParams)
	
	left		<- mfParams[1]
	mean		<- mfParams[2]
	right		<- mfParams[3]
	height		<- mfParams[4]

	y= pmax(pmin( (x-left)/(mean-left), (right-x)/(right-mean) ), 0) * height
	y[is.na(y)]= 1; 
	
	# Returns a list of the mfName, mfType, mfX, mfParams and mfVals.
	list(
		mfName= 	mfName,
		mfType= 	"trimf",
		mfX= 		x,
		mfParams= 	mfParams,
		mfVals= 	y
	)
}


#---------------------------------------------
#
# Fuzzy Inference System functions
#
#---------------------------------------------


newFIS <- function (	FISName, 
			FISType= 	"mamdani", 
			version= 	"1.0",
			andMethod= 	"min", 
			orMethod= 	"max", 
			impMethod= 	"min", 
			aggMethod= 	"max", 
			defuzzMethod= "centroid") {
# Inputs	: Strings FISName, FISType, andMethod, orMethod, impMethod, aggMethod and defuzzMethod
# Outputs	: Creation of a FIS

	# Validate FIS name input
	nameValidate(FISName)
	
	# Returns a FIS
	newFIS <- list(
		name = FISName, 
		type = FISType, 
		version = version,
		andMethod = andMethod, 
		orMethod = orMethod, 
		impMethod = impMethod, 
		aggMethod = aggMethod, 
		defuzzMethod = defuzzMethod, 
		inputList = NULL, 
		outputList = NULL, 
		ruleList = NULL)
}

addMF <- function (FIS, varType, varIndex, mf) {
# Inputs	: FIS (FIS) representing the FIS to have the membership function added to, 
#			varType (string) which can be either "input" or "output", 
#			varIndex (numeric integer) representing the index of the varType to be added to, and 
#			mf(mf) which is the membership function to be added.
# Outputs	: A FIS with the input membership function added to the appropriate variable.
	
	# The following block affects input variables.
	if(varType == "input") {
		if(is.null(FIS$inputList[[varIndex]]$membershipFunctionList[1]) ) {
		} else {
			for(i in 1:length(FIS$inputList[[varIndex]]$membershipFunctionList)) {
				if(mf$mfName == FIS$inputList[[varIndex]]$membershipFunctionList[[i]]$mfName) {
					stop("Membership functions should have a unique name within their variable\n")
				}
			}
		}
		if(varIndex <= length(FIS$inputList)) {
			# If there is no previous membership function created on the input variable,
			# then store this at index 1.
			if(is.null(FIS$inputList[[varIndex]]$membershipFunctionList[1])) {
				FIS$inputList[[varIndex]]$membershipFunctionList[1] = list(mf)
			# Otherwise, append this membership function to the end of the list 
			# of membership functions for this input variable.
			} else {
				FIS$inputList[[varIndex]]$membershipFunctionList = append(FIS$inputList[[varIndex]]$membershipFunctionList, list(mf))
			} 
		}
	# The following block affects output variables.
	} else if(varType == "output") {
		if(is.null(FIS$outputList[[varIndex]]$membershipFunctionList[1]) ) {
		} else {
			for(i in 1:length(FIS$outputList[[varIndex]]$membershipFunctionList)) {
				if(mf$mfName == FIS$outputList[[varIndex]]$membershipFunctionList[[i]]$mfName) {
					stop("Membership functions should have a unique name within their variable\n")
				}
			}
		}
		if(varIndex <= length(FIS$outputList)) {
			# If there is no previous membership function created on the output variable,
			# then store this at index 1.
			if(is.null(FIS$outputList[[varIndex]]$membershipFunctionList[1])) {
				FIS$outputList[[varIndex]]$membershipFunctionList[1] = list(mf)
			# Otherwise, append this membership function to the end of the list 
			# of membership functions for this input variable.
			} else {
				FIS$outputList[[varIndex]]$membershipFunctionList = append(FIS$outputList[[varIndex]]$membershipFunctionList, list(mf))
			}
		}
	# If the user has not specified the varType then output this message.
	} else {
		stop("You can enter only 'input' or 'output'\n")
	}
	FIS
}

addVar <- function (FIS, varType, varName, varBounds) {
# Inputs	: FIS (FIS) representing the FIS to have the variable added to, 
#			varType (String) which can be either "input" or "output", 
#			varName (String) representing the name of the variable, and 
#			varBounds (numeric vector) representing the boundaries of the variable (e.g. 1:10).
# Outputs	: A FIS with the variable added to it.
	
	# Validate FIS name input
	nameValidate(varName)
	
	# This block of code affects the variable type 'input'.
	if(varType == "input") {
		# Creates a list of inputs to be added.
		inputsToBeAdded <- list(inputName=varName, inputBounds=varBounds, membershipFunctionList= NULL)
		# If there are no input variables created for this FIS, then place this variable at index 1.
		if(is.null(FIS$inputList[1])) {
			FIS$inputList[1] <- list(inputsToBeAdded)
		# Otherwise append this variable to the end of the input variable list.
		} else {
			FIS$inputList <- append(FIS$inputList, list(inputsToBeAdded))
		}
	# This block of code affects the variable type 'output'.
	} else if (varType=="output") { 
		# Creates a list of outputs to be added.
		outputsToBeAdded <- list(outputName=varName, outputBounds=varBounds, membershipFunctionList = NULL)
		# If there are no output variables created for this FIS, then place this variable at index 1.
		if(is.null(FIS$outputList[1])) {
			FIS$outputList[1] <- list(outputsToBeAdded)
		# Otherwise append this variable to the end of the output variable list.
		} else {
			FIS$outputList <- append(FIS$outputList, list(outputsToBeAdded))	
		}
	# If the user has not specified the varType then output this message.
	} else {
		stop("You can only enter 'input' or 'output'\n")
	}
	FIS
}

addRule <- function(FIS, inputRule) {
# Inputs 	: A FIS (FIS), and a numeric vector of length m+n+2 where m represents the amount of inputs,
#			n represents the amount of outputs, and 2 represents the weight and connective.
# Outputs	: A FIS with the rule added
  
	# Takes a vector (inputRule) and binds this to a matrix (ruleList)
	FIS$ruleList <- rbind(FIS$ruleList, rbind(inputRule))
	rownames(FIS$ruleList) <- NULL

	FIS  
}

evalMF <- function(x, mfParams, mfType) {
# Inputs	: x(numeric vector) which represents the variable bounds/range, 
#			mfParams (numeric vector) which contains the given parameters of the membership function to be evaluated,
#			mfType (String) which contains the type of membership function (i.e. gaussMF, gaussbMF etc.), 
#			mfName (String) representing the name for the membership function
# Outputs	: Membership function with evaluated values

	# Depending on the type of the Membership Function, create the appropriate Membership Function with the values entered
	switch(mfType,
		"gaussmf" = 	(gaussMF("mfName", x, mfParams))$mfVals,
		"gaussbmf"= 	(gaussbMF("mfName", x, mfParams))$mfVals,
		"trapmf"= 	(trapMF("mfName", x, mfParams))$mfVals,
		"trimf"= 	(triMF("mfName", x, mfParams))$mfVals)	
}


removeVar <- function(FIS, varType, varIndex) {
# Inputs	: FIS(FIS) representing the FIS structure in question, 
#			varType(String) representing the variable type ('input' or 'output'), 
#			varIndex (Integer) representing the idnex value of the variable to be removed.
# Outputs	: A FIS with the variable removed.
	
	# This block of code affects the input variables
	if(varType == "input") {
		# Sets the variable at the specified index to NULL which removes all the data from the FIS
		FIS$inputList[[varIndex]]= NULL
		# If the input list is already empty, then keep it that way.
		if(length(FIS$inputList) == 0) {
			FIS$inputList= NULL
		}
	# This block of code affects the output variables
	} else if (varType == "output") {
		# Sets the variable at the specified index to NULL which removes all the data from the FIS	
		FIS$outputList[[varIndex]] <- NULL
		# If the output list is already empty, then keep it that way.
		if(length(FIS$outputList) == 0) {
			FIS$outputList= NULL
		}
	# If the user has not specified the varType then output this message.
	} else {
		stop("Must be 'input' or 'output' only \n")
	}

	FIS
}

removeMF <- function(FIS, varType, varIndex, mfIndex) {
# Inputs	: FIS(FIS) representing the FIS structure in question, 
#			varType(String) representing the variable type ('input' or 'output'), 
#			varIndex (Integer) representing the idnex value of the variable, and 
#			mfIndex (Integer) representing the index value of the membership function to be removed.
# Outputs	: A FIS with the membership function removed.
	
	# This block of code affects the input variables
	if(varType == "input") {
		# Sets the membership function at the specified index to NULL which removes all the data from the FIS
		FIS$inputList[[varIndex]]$membershipFunctionList[[mfIndex]]= NULL
		# If the input list is already empty, then keep it that way.
		if(length(FIS$inputList[[varIndex]]$membershipFunctionList) == 0) {
			FIS$inputList[[varIndex]]$membershipFunctionList= NULL
		}
	# This block of code affects the input variables
	} else if(varType == "output") {
	# Sets the membership function at the specified index to NULL which removes all the data from the FIS
		FIS$outputList[[varIndex]]$membershipFunctionList[[mfIndex]]= NULL
		# If the output list is already empty, then keep it that way.
		if(length(FIS$outputList[[varIndex]]$membershipFunctionList) == 0) {
			FIS$outputList[[varIndex]]$membershipFunctionList= NULL
		}
	# If the user has not specified the varType then output this message.
	} else {
		stop("Must be 'input' or 'output' only \n")
	}

	FIS
}

#---------------------------------------------
#
# Defuzzification function
#
#---------------------------------------------


defuzz <- function(x, vals, type) {
# Inputs	: x (vector) representing the processed values of the membership function
#			vals (vector) representing the range of values of the membership function
#			type (String) representing the type of defuzzification to be applied
# Outputs	: Defuzzified crisp value for the membership function given

	# Depending on the defuzzification method type, execute appropriately to get the crisp value.
	if(type == "centroid") {
		sum(vals * x) / sum(vals)
	} else if(type == "bisector") {
		tmp = 0
		for(i in 1:length(x)) {
			tmp = tmp + vals[i]
			if(tmp > sum(vals)/2) {
				break
			}
		}
		x[i]
	} else if(type == "mom") {
		mean(x[which(vals == max(vals))])
	} else if(type == "som") {
		x[min(which(vals == max(vals)))]
	} else if(type == "lom") {
		x[max(which(vals==max(vals)))]
	
	} else {
		# If the user has not specified a valid type then output this message.
		stop("Unsupported defuzzification function type\n")
	}
}


#---------------------------------------------
#
# File IO functions
#
#---------------------------------------------


writeFIS <- function(FIS, fileName) {
#Inputs	: FIS (FIS) of the FIS to be written to file, 
#			fileName(String) is an absolute path to the file to be written 
#			(or if non-existing, created and written to)
#Outputs	: A file at the specified directory
	
	# Validate FIS name input
	nameValidate(fileName)
	
	# If the file name does not end in .fis then append .fis to the end
	if(grepl("^.*?\\.fis", fileName) == FALSE) {
		fileName = paste(fileName, ".fis", sep="")
	}
	
	inputCount =	length(FIS$inputList)
	outputCount=	length(FIS$outputList)
	ruleCount = 	length(FIS$ruleList[,1])
	
	# Setting basic FIS information, and setting the buffer -txt- to NULL.
	txt=	  NULL
	txt[1]=  paste("[System]",sep="\n")
	txt[2]=  paste("Name='",FIS$name,"'",sep="")
	txt[3]=  paste("Type='",FIS$type,"'",sep="")
	txt[4]=  paste("Version=",FIS$version,sep="")
	txt[5]=  paste("NumInputs=",inputCount,sep="")
	txt[6]=  paste("NumOutputs=",outputCount,sep="")
	txt[7]=  paste("NumRules=",ruleCount,sep="")
	txt[8]=  paste("AndMethod='",FIS$andMethod,"'",sep="")
	txt[9]=  paste("OrMethod='",FIS$orMethod,"'",sep="")
	txt[10]= paste("ImpMethod='",FIS$impMethod,"'",sep="")
	txt[11]= paste("AggMethod='",FIS$aggMethod,"'",sep="")
	txt[12]= paste("DefuzzMethod='",FIS$defuzzMethod,"'",sep="")
	
	txtc = length(txt)+1
	txt[txtc] = paste("",sep="")
	txtc = length(txt)+1
	
	# If the input list is not null, paste the input variable's data into the buffer.
	if(!is.null(FIS$inputList)){
		for(i in 1:inputCount) {
			mfCount =	length(FIS$inputList[[i]]$membershipFunctionList)
			
			txt[txtc]=	paste("[Input",i,"]", sep="")
			txtc=txtc+1

			txt[txtc]=	paste("Name='",FIS$inputList[[i]]$inputName,"'", sep="")
			txtc= txtc+1

			txt[txtc]=	paste("Range=[",FIS$inputList[[i]]$inputBounds[[1]]," ",
					FIS$inputList[[i]]$inputBounds[[length(FIS$inputList[[i]]$inputBounds)]],"]",sep="")
			txtc=txtc+1

			txt[txtc]=	paste("NumMFs=",mfCount, sep="")
			txtc= txtc+1
			
			# The following block of code will -assuming that membership functions exist- paste 
			# every membership function of an input variable into the buffer.
			if(!is.null(FIS$inputList[[i]]$membershipFunctionList)){
				for(j in 1:mfCount) {			
					segment0= paste("MF",j,"=",sep="")
					segment1= paste("'",FIS$inputList[[i]]$membershipFunctionList[[j]]$mfName,"'", sep="")
					segment2= paste("'",FIS$inputList[[i]]$membershipFunctionList[[j]]$mfType,"'", sep="")
					segment3= paste(FIS$inputList[[i]]$membershipFunctionList[[j]]$mfParams, sep="", collapse=" ")
					txt[txtc]= paste(segment0,segment1, ":", segment2, ",[", segment3, "]", sep="", collapse=" ")
					txtc=txtc+1
				}
			}
			txt[txtc]=""
			txtc=txtc+1
		}
	}
	# If the output list is not null, paste the output variable's data into the buffer.
	if(!is.null(FIS$outputList)){
		for(i in 1:outputCount) {
			mfCount =	length(FIS$outputList[[i]]$membershipFunctionList)
			txt[txtc]=	paste("[Output",i,"]", sep="")
			txtc=txtc+1

			txt[txtc]=	paste("Name='",FIS$outputList[[i]]$outputName,"'", sep="")
			txtc= txtc+1

			txt[txtc]=	paste("Range=[",FIS$outputList[[i]]$outputBounds[[1]]," ",
					FIS$outputList[[i]]$outputBounds[[length(FIS$outputList[[i]]$outputBounds)]],"]",sep="")
			txtc=txtc+1

			txt[txtc]=	paste("NumMFs=",mfCount, sep="")
			txtc= txtc+1
			
			# The following block of code will -assuming that membership functions exist- paste 
			# every membership function of an output variable into the buffer.
			if(!is.null(FIS$outputList[[i]]$membershipFunctionList)){
				for(j in 1:mfCount) {			
					segment0= paste("MF",j,"=",sep="")
					segment1= paste("'",FIS$outputList[[i]]$membershipFunctionList[[j]]$mfName,"'", sep="")
					segment2= paste("'",FIS$outputList[[i]]$membershipFunctionList[[j]]$mfType,"'", sep="")
					segment3= paste(FIS$outputList[[i]]$membershipFunctionList[[j]]$mfParams, sep="", collapse=" ")
					txt[txtc]= paste(segment0, segment1, ":", segment2, ",[", segment3, "]", sep="", collapse=" ")
					txtc=txtc+1
				}
			}
			txt[txtc]=""
			txtc=txtc+1
		}
	}
	
	txt[txtc]=	paste("[Rules]")
	txtc=txtc+1
	# If rules do exist, then paste every rule into the buffer.
	if(!is.null(FIS$ruleList)) {
		for(i in 1:ruleCount) {
			segment1= paste(FIS$ruleList[i,1:inputCount], collapse=" ")
			segment2= paste(FIS$ruleList[i,(inputCount+1):(inputCount+outputCount)],collapse=" ")
			segment3= paste(" (",FIS$ruleList[i,inputCount+outputCount+1], ") : ",FIS$ruleList[i,inputCount+outputCount+2], sep="")
			txt[txtc]= paste(segment1,", ",segment2,segment3,sep="")
			txtc=txtc+1
		}
	}
	
	# If the file does not already exist, output a message saying it will create a new one.
	if(!file.exists(fileName)) {
		cat("The specified file does not exist, creating a new file\n")
	}
	
	# Write the buffer to file.
	fileConn= file(fileName)
	writeLines(txt, fileConn)
	close(fileConn)
}


showFIS <- function(FIS) {
# Inputs	: FIS (FIS) which will be shown via console (as in, all of its data)
# Outputs	: Text to the console in an ordered format

	inputCount =	length(FIS$inputList)
	outputCount=	length(FIS$outputList)
	ruleCount = 	length(FIS$ruleList[,1])

	# Setting basic FIS information, and setting the buffer -txt- to NULL.
	txt=	  NULL
	txt[1]=  paste("[System]",sep="\n")
	txt[2]=  paste("Name='",FIS$name,"'",sep="")
	txt[3]=  paste("Type='",FIS$type,"'",sep="")
	txt[4]=  paste("Version=",FIS$version,sep="")
	txt[5]=  paste("NumInputs=",inputCount,sep="")
	txt[6]=  paste("NumOutputs=",outputCount,sep="")
	txt[7]=  paste("NumRules=",ruleCount,sep="")
	txt[8]=  paste("AndMethod='",FIS$andMethod,"'",sep="")
	txt[9]=  paste("OrMethod='",FIS$orMethod,"'",sep="")
	txt[10]= paste("ImpMethod='",FIS$impMethod,"'",sep="")
	txt[11]= paste("AggMethod='",FIS$aggMethod,"'",sep="")
	txt[12]= paste("DefuzzMethod='",FIS$defuzzMethod,"'",sep="")
	 
	txtc = length(txt)+1
	txt[txtc] = paste("",sep="")
	txtc = length(txt)+1
	
	# If the input list is not null, paste the input variable's data into the buffer.
	if(!is.null(FIS$inputList)){
		for(i in 1:inputCount) {
			mfCount =	length(FIS$inputList[[i]]$membershipFunctionList)
			
			txt[txtc]=	paste("[Input",i,"]", sep="")
			txtc=txtc+1

			txt[txtc]=	paste("Name='",FIS$inputList[[i]]$inputName,"'", sep="")
			txtc= txtc+1

			txt[txtc]=	paste("Range=[",FIS$inputList[[i]]$inputBounds[[1]],":",
					FIS$inputList[[i]]$inputBounds[[length(FIS$inputList[[i]]$inputBounds)]],"]",sep="")
			txtc=txtc+1

			txt[txtc]=	paste("NumMFs=",mfCount, sep="")
			txtc= txtc+1
			
			# The following block of code will -assuming that membership functions exist- paste 
			# every membership function of an input variable into the buffer.
			if(!is.null(FIS$inputList[[i]]$membershipFunctionList)){
				for(j in 1:mfCount) {			
					segment0= paste("MF",j,"=",sep="")
					segment1= paste("'",FIS$inputList[[i]]$membershipFunctionList[[j]]$mfName,"'", sep="")
					segment2= paste("'",FIS$inputList[[i]]$membershipFunctionList[[j]]$mfType,"'", sep="")
					segment3= paste(FIS$inputList[[i]]$membershipFunctionList[[j]]$mfParams, sep="", collapse=" ")
					txt[txtc]= paste(segment0,segment1, ":", segment2, ",[", segment3, "]", sep="", collapse=" ")
					txtc=txtc+1
				}
			}
			txt[txtc]=""
			txtc=txtc+1
		}
	}
	# If the output list is not null, paste the input variable's data into the buffer.
	if(!is.null(FIS$outputList)){
		for(i in 1:outputCount) {
			mfCount =	length(FIS$outputList[[i]]$membershipFunctionList)
			txt[txtc]=	paste("[Output",i,"]", sep="")
			txtc=txtc+1

			txt[txtc]=	paste("Name='",FIS$outputList[[i]]$outputName,"'", sep="")
			txtc= txtc+1

			txt[txtc]=	paste("Range=[",FIS$outputList[[i]]$outputBounds[[1]],":",
					FIS$outputList[[i]]$outputBounds[[length(FIS$outputList[[i]]$outputBounds)]],"]",sep="")
			txtc=txtc+1

			txt[txtc]=	paste("NumMFs=",mfCount, sep="")
			txtc= txtc+1
			
			# The following block of code will -assuming that membership functions exist- paste 
			# every membership function of an output variable into the buffer.
			if(!is.null(FIS$outputList[[i]]$membershipFunctionList)){
				for(j in 1:mfCount) {			
					segment0= paste("MF",j,"=",sep="")
					segment1= paste("'",FIS$outputList[[i]]$membershipFunctionList[[j]]$mfName,"'", sep="")
					segment2= paste("'",FIS$outputList[[i]]$membershipFunctionList[[j]]$mfType,"'", sep="")
					segment3= paste(FIS$outputList[[i]]$membershipFunctionList[[j]]$mfParams, sep="", collapse=" ")
					txt[txtc]= paste(segment0, segment1, ":", segment2, ",[", segment3, "]", sep="", collapse=" ")
					txtc=txtc+1
				}
			}
			txt[txtc]=""
			txtc=txtc+1
		}
	}
	
	txt[txtc]=paste("[Rules]")
	txtc=txtc+1
	
	# If rules do exist, then paste every rule into the buffer.
	if(!is.null(FIS$ruleList)) {
		for(i in 1:ruleCount) {
			segment1= paste(FIS$ruleList[i,1:inputCount], collapse=" ")
			segment2= paste(FIS$ruleList[i,(inputCount+1):(inputCount+outputCount)],collapse=" ")
			segment3= paste(" (",FIS$ruleList[i,inputCount+outputCount+1], ") : ",FIS$ruleList[i,inputCount+outputCount+2], sep="")
			txt[txtc]= paste(segment1,", ",segment2,segment3,sep="")
			txt[txtc]= gsub("\\-","\\!",txt[txtc])
			txtc=txtc+1
		}
	}
txt	
}


readFIS <- function(fileName) {
#Inputs	: fileName (string) which is an absolute path to the file 
#			(with a .fis extension) to be read into memory
#Outputs	: a FIS structure

	# Validate FIS name input
	nameValidate(fileName)
	
	# If the file name does not end in .fis then append .fis to the end
	if(grepl("^.*?\\.fis", fileName) == FALSE) {
		fileName = paste(fileName, ".fis", sep="")
	}
	
	# Read in all the lines from the given file
	txt= readLines(fileName)
	
	# See if 'lc' matches '[System]', otherwise stop.
	lc= charmatch("[System]", txt)
	if (is.na(lc) || lc == 0){
		stop('No \'[System Variables]\' line in file', fileName)
	}
	
	# Format the basic FIS information structure.
	for(i in 1:12) {
		txt[i]= gsub("^.*?=('| |)","", txt[i])
		txt[i]= gsub("'","", txt[i])
	}
	
	# Parse the appropriate variable count into input count, output count, and rule count.
	inputCount = eval(parse(text=txt[5]))
	outputCount= eval(parse(text=txt[6]))
	ruleCount=	eval(parse(text=txt[7]))
	
	# Creates placeholder FIS for the coming values.
	FIS <- newFIS("temp")	
	
	# Overwrite the values in the placeholder FIS with the values from the file.
	FIS$name=	 	txt[2]
	FIS$type=	 	txt[3]
	FIS$version=	 	eval(parse(text=txt[4]))
	FIS$andMethod= 	txt[8]
	FIS$orMethod=	 	txt[9]
	FIS$impMethod= 	txt[10]
	FIS$aggMethod= 	txt[11]
	FIS$defuzzMethod=	txt[12]
	
	# The following line counts how many input variables exist from the file by pattern matching.
	inputLines= grep("\\[Input.\\]", txt)
	# Checks to see if any inputs exist, if not, ignore code block.
	if(length(inputLines > 0)) {
		# The following loop obtains each variable's name and range.
		for(i in 1:length(inputLines)) {
			txtc = inputLines[[i]]
			mfCount = eval(parse(text=txt[txtc+3]))
			txtc= txtc+1
			varHolder= list()
			# The following loop formats the lines appropriately (to meet the set standard for our .fis files).
			for(j in 1:2) {
				txt[txtc]= gsub("\\[|\\]","", gsub("'", "", gsub("^.*?=('|)","",txt[txtc])))
				txt[txtc]= gsub(" ", ":", txt[txtc])
				varHolder= append(varHolder,txt[txtc])
				txtc= txtc+1
			}
			# Adds the variable name and range to the FIS structure currently stored in memory.
			FIS= addVar(FIS, "input", varHolder[[1]], eval(parse(text=varHolder[[2]])))
			
			# The following block reads all the data from each of the variable's membership function from the file,
			# and stores the data in relevant variables for a later addition to the FIS structure in memory which 
			# occurs on every iteration.
			if(mfCount>0) {
				txtc=txtc+1
				mfHolder= list(mfName="", mfType="", mfParams=c())
				for(j in 1:mfCount) {
					mfVal = txt[txtc]
					mfVal = gsub("'.*$","", gsub("^.*?'","", mfVal))
					mfHolder$mfName[[j]] = mfVal
					
					mfVal = txt[txtc]
					mfVal = gsub("'.*$","",gsub("^.*:'","",mfVal))
					mfHolder$mfType[[j]] = mfVal
					
					mfVal = txt[txtc]
					mfVal = gsub("^.*\\[","",gsub("\\]$","",mfVal))
					mfVal = strsplit(mfVal, " ")
					
					vectorStore= c()
					for(k in 1:length(mfVal[[1]])) {
						vectorStore = append(vectorStore, mfVal[[1]][k])
					}
					vectorStore= as.numeric(vectorStore)
					mfHolder$mfParams[[j]] = vectorStore
					
					if(mfHolder$mfType[[j]] == "gaussmf") {
						returnMF = 	gaussMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else if(mfHolder$mfType[[j]] == "gaussbmf") {
						returnMF = 	gaussbMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else if(mfHolder$mfType[[j]] == "trimf") {
						returnMF = 	triMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else {
						returnMF = 	trapMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					}
					
					# Adds the membership function to the input variable.
					FIS = addMF(FIS, "input", i, returnMF)
					txtc=txtc+1
				}
			}			
		}
	}
	
	# The following line counts how many output variables exist from the file by pattern matching.
	outputLines= grep("\\[Output.\\]", txt)
	# Checks to see if any outputs exist, if not, ignore code block.
	if(length(outputLines > 0)) {
		# The following loop obtains each variable's name and range.
		for(i in 1:length(outputLines)) {
			txtc = outputLines[[i]]
			mfCount = eval(parse(text=txt[txtc+3]))
			txtc= txtc+1
			# The following loop formats the lines appropriately (to meet the set standard for our .fis files).
			varHolder= list()
			# Adds the variable name and range to the FIS structure currently stored in memory.
			for(j in 1:2) {
				txt[txtc]= gsub("\\[|\\]","", gsub("'", "", gsub("^.*?=('|)","",txt[txtc])))
				txt[txtc]= gsub(" ", ":", txt[txtc])
				varHolder= append(varHolder,txt[txtc])
				txtc= txtc+1
			}
			# Adds the output variable to the FIS structure in memory.
			FIS= addVar(FIS, "output", varHolder[[1]], eval(parse(text=varHolder[[2]])))
			
			# The following block reads all the data from each of the variable's membership function from the file,
			# and stores the data in relevant variables for a later addition to the FIS structure in memory which 
			# occurs on every iteration.
			if(mfCount>0) {
				txtc=txtc+1
				mfHolder= list(mfName="", mfType="", mfParams=c())
				for(j in 1:mfCount) {
					mfVal = txt[txtc]
					mfVal = gsub("'.*$","", gsub("^.*?'","", mfVal))
					mfHolder$mfName[[j]] = mfVal
					
					mfVal = txt[txtc]
					mfVal = gsub("'.*$","",gsub("^.*:'","",mfVal))
					mfHolder$mfType[[j]] = mfVal
					
					mfVal = txt[txtc]
					mfVal = gsub("^.*\\[","",gsub("\\]$","",mfVal))
					mfVal = strsplit(mfVal, " ")
					#mfVal = as.numeric(mfVal)
					
					vectorStore= c()
					for(k in 1:length(mfVal[[1]])) {
						vectorStore = append(vectorStore, mfVal[[1]][k])
					}
					vectorStore= as.numeric(vectorStore)
					mfHolder$mfParams[[j]] = vectorStore
					
					if(mfHolder$mfType[[j]] == "gaussmf") {
						returnMF = 	gaussMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else if(mfHolder$mfType[[j]] == "gaussbmf") {
						returnMF = 	gaussbMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else if(mfHolder$mfType[[j]] == "trimf") {
						returnMF = 	triMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					} else {
						returnMF = 	trapMF(mfHolder$mfName[[j]], 
								eval(parse(text=varHolder[[2]])), 
								mfHolder$mfParams[[j]])
					}
					# Adds the membership function to the input variable.
					FIS = addMF(FIS, "output", i, returnMF)
					
					txtc=txtc+1
				}
			}
		}
	}
	
	# Sets the buffer line counter to the rules section.
	txtc= grep("\\[Rules\\]", txt)
	txtc= txtc+1
	
	# If rules exist, format and add to the FIS structure rule matrix. If no rules exist, ignore
	# the following code block.
	if(ruleCount > 0) {
		for(i in 1:ruleCount) {
			iVals = txt[txtc]
			iVals = gsub(",.*$", "", iVals)
			iVals = strsplit(iVals, " ")
			
			iValVector = c()
			for(j in 1:length(iVals[[1]])) {
				iValVector = append(iValVector, as.numeric(iVals[[1]][j]))
			}
			
			oVals = txt[txtc]		
			oVals = gsub(".$","",gsub("\\(.*$","",gsub("^.*,.","",oVals)))
			oVals = strsplit(oVals, " ")
			
			oValVector = c()
			for(j in 1:length(oVals[[1]])) {
				oValVector = append(oValVector, as.numeric(oVals[[1]][j]))
			}
			
			wVal = txt[txtc]
			wVal = gsub(".*\\(", "", gsub("\\).*$","", wVal))
			wVal = as.numeric(gsub(" ","",wVal))
			
			fVal = txt[txtc]
			fVal = gsub("^.*:", "", fVal)
			fVal = as.numeric(gsub(" ","", fVal))
			
			f = c ()
			f = append(f, iValVector)
			f = append(f, oValVector)
			f = append(f, wVal)
			f = append(f, fVal)
			FIS = addRule(FIS, f)
			
			txtc= txtc + 1	
		}
	}
	# Return the FIS with all the values from the file.
	FIS
}


#---------------------------------------------
#
# Plotting functions
#
#---------------------------------------------


plotMF <- function(FIS, varType, varIndex) {
#Inputs	: FIS (FIS) representing a FIS which will be used, 
#			varType(String) which can be either "input" or "output" and 
#			varIndex (integer) representing the variable in the input/output whose 
#			membership functions will be plotted to a graph.
#Outputs	: A graphic containing the derived input data in graph format
	
	# Require the library, 'splines' for graphical chart creation.
	require(splines)
	# Set the plot character to nothing so it will not show on the graph.
	pchVal= ""
	# Set the y axis height.
	ylimVal= c(0,1.025)
	# Create a new plot.
	plot.new()
	# The following block corresponds to whether the variable type is 'input'.
	if(varType == "input") {
		# The following block plots all the membership functions of a specified input variable onto a graph
		plot.window(xlim=c(	FIS$inputList[[varIndex]]$inputBounds[1],FIS$inputList[[varIndex]]$inputBounds[length(FIS$inputList[[varIndex]]$inputBounds)]), ylim=ylimVal)
		numMFs= length(FIS$inputList[[varIndex]]$membershipFunctionList)
		for(i in 1:numMFs) {
			colorRGB= runif(3,20,160)
			# mfList is a convenience variable in that it saves a lot of extra code to access the same data.
			mfList= FIS$inputList[[varIndex]]$membershipFunctionList[[i]]
			if(mfList$mfType == 'gaussmf' || mfList$mfType == 'gaussbmf') {
				curvePredict= predict(interpSpline(mfList$mfX, mfList$mfVals))
				lines(curvePredict, col=colorRGB, type="o", xlim=c(1,length(FIS$inputList[[varIndex]]$inputBounds)), ylim=c(0,1.1), ann=FALSE, pch=pchVal)
				text(match(TRUE,mfList$mfVals==max(mfList$mfVals))-1,1.025,mfList$mfName)
			} else {
				lines(mfList$mfX, mfList$mfVals, type="o", col=colorRGB, xlim=c(0,length(FIS$inputList[[varIndex]]$inputBounds)), ylim=c(0,1.1), ann=FALSE, pch=pchVal)
				text(match(TRUE,mfList$mfVals==max(mfList$mfVals))-1,1.025,mfList$mfName)
			}
		}
		title(paste("Membership functions from input variable '",FIS$inputList[[varIndex]]$inputName,"'", sep=""))
	} else if(varType == "output") {
		# The following block plots all the membership functions of a specified output variable onto a graph
		plot.window(xlim=c(FIS$outputList[[varIndex]]$outputBounds[1],FIS$outputList[[varIndex]]$outputBounds[length(FIS$outputList[[varIndex]]$outputBounds)]), ylim=ylimVal)
		numMFs= length(FIS$outputList[[varIndex]]$membershipFunctionList)
		for(i in 1:numMFs) {
			colorRGB= runif(3,20,160)
			# mfList is a convenience variable in that it saves a lot of extra code to access the same data.
			mfList= FIS$outputList[[varIndex]]$membershipFunctionList[[i]]
			if(mfList$mfType == 'gaussmf' || mfList$mfType == 'gaussbmf') {
				curvePredict= predict(interpSpline(mfList$mfX, mfList$mfVals))
				lines(curvePredict, col=colorRGB, type="o", xlim=c(0,length(FIS$outputList[[varIndex]]$outputBounds)), ylim=c(0,1.1), ann=FALSE, pch=pchVal)
				text(match(TRUE,mfList$mfVals==max(mfList$mfVals))-1,1.025,mfList$mfName)
			} else {
				lines(mfList$mfX, mfList$mfVals, type="o", col=colorRGB, xlim=c(0,length(FIS$outputList[[varIndex]]$outputBounds)), ylim=c(0,1.1), ann=FALSE, pch=pchVal)
				text(match(TRUE,mfList$mfVals==max(mfList$mfVals))-1,1.025,mfList$mfName)
			}
		}
		title(paste("Membership functions from output variable '",FIS$outputList[[varIndex]]$outputName,"'", sep=""))
	} else {
		stop("Must be either 'input' or 'output'\n")
	}
	# Plot axes, axex title.
	axis(1)
	axis(2)
	title(xlab="Range")
	title(ylab="Degree of membership")
	box()
}

#---------------------------------------------
#
# FIS Evaluation Function
#
#---------------------------------------------


evalFIS <- function(inputStack, fis, numPoints=101) {
#Inputs	: Input stack, FIS object, and an integer number of points to plot against
#Outputs	: A evaluated and defuzzified crisp value for a FIS
	
	#Set number of points, default = 101
	cat("Setting number of points to ", numPoints, "\n", sep="")
	
	#Initialize FIS variables
	fisType = fis$type 													
	numInputs = length(fis$inputList) 									
	numOutputs = length(fis$outputList) 								
	inputMFCount = NULL 												
	outputMFCount = NULL 												
	numRules = nrow(fis$ruleList) 										
	ruleInputs = fis$ruleList[,1:length(fis$inputList)] 				
	andOr = fis$ruleList[,numInputs+numOutputs+2] 						
	ruleWeight = fis$ruleList[,numInputs+numOutputs+1]					
	ruleConns = fis$ruleList[, (numInputs+1):(numInputs+numOutputs)] 	

	# Collect number of membership functions in each input variable
	for(i in 1:numInputs) {
		inputMFCount[i] = length(fis$inputList[[i]]$membershipFunctionList)
	}

	# Collect number of membership functions in each output variable
	for(i in 1:numOutputs) {
		outputMFCount[i] = length(fis$outputList[[i]]$membershipFunctionList)
	}
	
	outputRange = matrix(0, numOutputs, 2)
	outputMfMatrix = matrix(0, sum(outputMFCount)+1, numPoints)      		
	outputMfMatrix[1,] = 1
	o_index = 1
	
	for(i in 1:numOutputs) {
		for(j in 1: outputMFCount[i]){
		# Evaluates each output membership function from each output
			p = fis$outputList[[i]]$membershipFunctionList[[j]]
			o = fis$outputList[[i]]$outputBounds
			outputRange[i,] = c(o[1], o[length(o)])
			
			outputMfMatrix[o_index+1, ] =
 				evalMF(seq(outputRange[i,1], outputRange[i,2], length=numPoints), 
 					p$mfParams, 
 					p$mfType)

			o_index = o_index + 1
		}
	}
	
	o_index = abs(ruleConns)+matrix((0:(numOutputs-1))*numRules+1, 
		numRules, numOutputs, byrow=TRUE)
	o_index[ruleConns==0] = 1
		
	outputMfVals = outputMfMatrix[t(o_index),]					
	outputMfVals[t(ruleConns) < 0,] = 1 - outputMfVals[t(ruleConns<0)]
	outputMfVals = matrix(t(outputMfVals), numRules, numPoints * numOutputs, byrow=TRUE)
 	
	finalOutputValues = matrix(0, numRules, numPoints * numOutputs)			
	totalOutputValues = matrix(0, numPoints * numOutputs)					

	if (is.vector(inputStack)) {
		inputStack= rbind(inputStack)
	}
	data_n= nrow(inputStack)

	#Assign output stack
	outputStack = matrix(0, data_n, numOutputs)
		
	for(i in 1:data_n) {

		input = inputStack[i,]
		mfVal = rep(0, sum(inputMFCount))
		v_index = 1

		for(j in 1:(numInputs)) {			
			for(k in 1:inputMFCount[[j]]) {
				mfVal[v_index] = evalMF(input[j], 
							fis$inputList[[j]]$membershipFunctionList[[k]]$mfParams, 
							fis$inputList[[j]]$membershipFunctionList[[k]]$mfType)
				v_index = v_index + 1		
			}
		}
			
		mfVal = c(0, mfVal)

		#Create rule matrix
		m_index= matrix(1,numRules,1) %*% cumsum(c(0, inputMFCount[1:(numInputs-1)]))+abs(ruleInputs)+1
		inputMfValues= matrix(mfVal[m_index], numRules, numInputs)
		
		#Assign values of connectives
		inputMfValues[which(((andOr == 1) * (ruleInputs == 0)) == 1)] = 1
		inputMfValues[which(((andOr == 2) * (ruleInputs == 0)) == 1)] = 0

	    	# Deal with inverted rules
   		v_index= which(ruleInputs < 0)
   		inputMfValues[v_index]= 1 - inputMfValues[v_index]
    	
		# Assign weight (firing strength) of rules
		weight = matrix(0, numRules, 1)
	
		if ( numInputs == 1){
			weight = inputMfValues
		} else {
			andex = which(andOr == 1)
			ordex = which(andOr == 2)
				
			weight[andex] = apply(rbind(inputMfValues[andex,]), 1, fis$andMethod)
			weight[ordex] = apply(rbind(inputMfValues[ordex,]), 1, fis$orMethod)
		} # end if ( numInputs == 1)

		weight = weight * ruleWeight
	
		if(fisType == "mamdani") {
			#Convert input stack to matrix if not already
			# Create matrix of weights against number of points defined above
			temp = matrix(weight, nrow(weight), numPoints * numOutputs)

			if( fis$impMethod == 'prod' ){	
				finalOutputValues = temp * outputMfVals
			} else if ( fis$impMethod == 'min' ){
			finalOutputValues = pmin(temp, outputMfVals)
			} else {
				stop("Unsupported implication method\n")
			}

			totalOutputValues = apply(finalOutputValues, 2, fis$aggMethod)

  			# Defuzzify each of the outputs
  			for( j in 1:numOutputs){
				outputStack[i, j] = defuzz(
					seq(outputRange[j,1], outputRange[j,2], length = numPoints),
					totalOutputValues[((j-1)*numPoints+1):(j*numPoints)], 
					fis$defuzzMethod)
      			}
		} else {
			# Room for expansion in future projects
			stop("Unsupported FIS type\n")
		} # end fisType == "mamdani"
	} # end for(i in 1:data_n)
	
	#Return the output stack
	outputStack
 }


meshgrid <- function(a,b) {
#Inputs		: a and b, both sets of points
#Outputs	: Union of a and b
  
	list(x=outer(b*0, a, FUN="+"), y=outer(b, a*0, FUN="+"))

}


gensurf <- function(fis, ix1=1, ix2=2, ox1=1) {
#Inputs		: A FIS structure
#Outputs	: A 3-D graph with two inputs on the x and y axes, and one ouput on the z
	
	i1 = fis$inputList[[ix1]]
	i2 = fis$inputList[[ix2]]
	o1 = fis$outputList[[ox1]]

	i1b = i1$inputBounds
	i2b = i2$inputBounds

	i1_min = i1b[1]
	i1_max = i1b[length(i1b)]

	i2_min = i2b[1]
	i2_max = i2b[length(i2b)]
  
	x_values = seq(i1_min, i1_max, length = 15)
	y_values = seq(i2_min, i2_max, length = 15)
	
	m_values = meshgrid(x_values, y_values)

	o_values = evalFIS(cbind(c(m_values$x), c(m_values$y)), fis)


	z_values = matrix(o_values[,ox1], 15, 15, byrow=TRUE)

  	h_values = (z_values[-15,-15] + z_values[-1,-15] + z_values[-15,-1] + z_values[-1,-1]) / 4
  	h_values= floor((h_values-min(h_values))/(max(h_values)-min(h_values))*14+.5)+1

 	persp(x_values, y_values, z_values, 
    	xlab=i1$inputName, ylab=i2$inputName, zlab=o1$outputName, 
    	theta=-30, phi=30, col=rainbow(15)[16-h_values], ticktype='detailed')
}


#---------------------------------------------
#
# Input Validation Functions
#
#---------------------------------------------


mfValidate <- function(mfName, mfParams) {
#Inputs	: mfName (String) representing name of membership function, 
#			mfParams (numeric vector) of the input parameters
#Outputs	: None 
	
	err_mfParamsNumnerNotMatched = "Incorrect amount for mfParams\n"
	
	# Checks the call stack to see if a membership function is being
	# called manually or within another function. If the former, then
	# validate user input on the name.
	if(is.null(sys.call(-2))) {
		#errt means error text (for storing error messages)
		errt= NULL
		#errc is the line counter for the error text.
		errc= 1
		
		if(nchar(mfName) == 0) {
			errt[errc]= "Name cannot be empty"
			errc= errc+1
		}
		
		if(grepl("\\$|\\[|\\]|\\|", mfName)) {
			errt[errc]= "Illegal characters entered"
			errc= errc+1
		}
		
		if(typeof(mfName) != "character") {
			errt[errc] = "Must enter a string for the name"
			errc= errc+1
		}
		# The following block checks which function was called, and then 
		# ensures that the mfParams are of the correct length.
		if(as.character(sys.call(-1)[[1]]) == "gaussbMF") {
			if(length(mfParams) != 5) {
				errt[errc]= err_mfParamsNumnerNotMatched
				errc= errc+1
			}
		} else if (as.character(sys.call(-1)[[1]]) == "gaussMF") {
			if(length(mfParams) != 3) {
				errt[errc]= err_mfParamsNumnerNotMatched
				errc= errc+1
			}
		} else if (as.character(sys.call(-1)[[1]]) == "trapMF") {
			if(length(mfParams) != 5) {
				errt[errc]= err_mfParamsNumnerNotMatched
				errc= errc+1
			}
		} else if (as.character(sys.call(-1)[[1]]) == "triMF") {
			if(length(mfParams) != 4) {
				errt[errc]= err_mfParamsNumnerNotMatched
				errc= errc+1
			}
		} 
		
		# If any errors are detected, then stop the execution of the functions
		# on the call stack and alert the user of whatever error occured.
		if(!is.null(errt)) {
			errt[errc]= "\n"
			stop(errt)
		}
	}
}

nameValidate <- function(name) {
	
	# errt means error text (for storing error messages)
	errt= NULL
	# errc is the line counter for the error text.
	errc= 1
	
	# If the given name from the user is of length 0, output appropriate message.
	if(nchar(name) == 0) {
		errt[errc]= "Name cannot be empty"
		errc= errc+1
	}
	# If the user enters illegal characters in the name, output appropriate message.
	if(grepl("\\$|\\[|\\]|\\|", name)) {
		errt[errc]= "Illegal characters entered"
		errc= errc+1
	}
	# If any errors are detected, then stop the execution of the functions
	# on the call stack and alert the user of whatever error occured.
	if(!is.null(errt)) {
		errt[errc]= "\n"
		stop(errt)
	}
}

#---------------------------------------------
#
# Misc. functions
#
#---------------------------------------------


tippertest <- function() {
#Inputs	: N/A
#Outputs	: A FIS structure which should be assigned to a variable
	
	FIS= newFIS('tippertest')
	
	FIS= addVar(FIS, "input", "service", 0:10)
	FIS= addVar(FIS, "input", "food", 0:10)
	FIS= addVar(FIS, "output", "tip", 0:30)

	#The following block is going to be input 1's membership functions
	mf1= gaussMF("poor", 0:10, c(1.5, 0, 1))
	mf2= gaussMF("good", 0:10, c(1.5, 5, 1))
	mf3= gaussMF("excellent", 0:10, c(1.5, 10, 1))
	
	#The following block is going to be input 2's membership functions
	mf4= trapMF("rancid", 0:10, c(0, 0, 1, 3, 1))
	mf5= trapMF("delicious", 0:10, c(7, 9, 10, 10, 1))
	
	#The following block is going to be output 1's membership functions
	mf6= triMF("cheap", 0:30, c(0, 5, 10, 1))
	mf7= triMF("average", 0:30, c(10, 15, 20, 1))
	mf8= triMF("generous", 0:30, c(20, 25, 30, 1))
	
	#The following three blocks simply add the membership functions to the relevant variables
	FIS= addMF(FIS, "input", 1, mf1)
	FIS= addMF(FIS, "input", 1, mf2)
	FIS= addMF(FIS, "input", 1, mf3)
	
	FIS= addMF(FIS, "input", 2, mf4)
	FIS= addMF(FIS, "input", 2, mf5)
	
	FIS= addMF(FIS, "output", 1, mf6)
	FIS= addMF(FIS, "output", 1, mf7)
	FIS= addMF(FIS, "output", 1, mf8)

	#The following block adds the rules to the fis structure
	FIS= addRule(FIS, c(1,1,1,1,2))
	FIS= addRule(FIS, c(2,0,2,1,1))
	FIS= addRule(FIS, c(3,2,3,1,2))
	
	#Return the FIS structure
	FIS
}