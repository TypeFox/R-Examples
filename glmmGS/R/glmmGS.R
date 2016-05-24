# Construct list of controls
glmmGS.Control = function(reltol = 1.e-6, abstol = 1.e-6, maxit = 200)
{
	x = list(reltol = reltol, abstol = abstol, maxit = maxit);
	return(x);
}

# Internal function to construct sparse matrix
glmmGS.CreateSparseMatrix = function(values, indices, counts)
{
	# Check types
	if ((is.integer(counts) == FALSE) || (is.vector(counts) == FALSE))
	{
		stop("\'counts\' must be a vector of integers");
	}
	if ((is.double(values) == FALSE) || (is.vector(values) == FALSE))
	{
		stop("\'values\' must be a vector of double");
	}
	if ((is.integer(indices) == FALSE) || (is.vector(indices) == FALSE))
	{
		stop("\'indices\' must be a vector of integers");
	}
	
	# Check length of counts
	if (length(counts) == 0)
	{
		stop("\'counts\' must have length greater than zero");
	}
	
	# Set number of columns
	ncols = length(counts) - 1L;
	
	# Check values of counts
	if (counts[1] != 0)
	{
		stop("\'counts[1]\' must be equal to zero");
	}
	for (j in 1:ncols)
	{
		if (counts[j + 1] < counts[j])
		{
			stop("\'counts\' must be monotonic non-descending");
		}
	}

	# Retrieve total number of non-zero entries
	nz = counts[ncols + 1];
	
	# Check length of values and indices
	if (length(values) != nz)
	{
		stop("Wrong length of \'values\'");
	}
	if (length(indices) != nz)
	{
		stop("Wrong length of \'indices\'");
	}
	
	# Check indices
	for (j in 1:ncols)
	{
		for (p in counts[j]:(counts[j + 1] - 1))
		{
			# Notice that p is zero-based
			if ((indices[p + 1L] < 0L) || (indices[p + 1L] >= ncols))
			{
				stop("\'indices\' values must be in [0, ncols - 1]")
			}
		}
	}
	
	# Set return list
	retval = list(values = values, indices = indices, counts = counts);
	
	# Attach "sparse" attribute to list
	attr(retval, "sparse.matrix") = TRUE;
	
	# Return
	return(retval);
}

# Construct sparse matrix
glmmGS.SparseMatrix = function(...)
{
	# Get argument list
	ls = list(...);
	
	# Initialize return value
	retval = NULL;
	
	# Parse argument list
	if (length(ls) == 1L)
	{
		R = ls[[1]];
		if ((is.null(R) == TRUE) || (is.matrix(R) == FALSE) || (nrow(R) != ncol(R)))
		{
			stop("Invalid matrix");
		}
		else
		{
			# Build sparse matrix
			# 1) Count non-zero elements
			ncols = ncol(R);
			counts = integer(ncols + 1L);
			count.total = 0L;
			counts[1] = 0L;
			for (j in 1:ncols)
			{
				nz = which(R[, j] != 0);
				count.total = count.total + length(nz);
				counts[j + 1] = count.total;
			}
			
			# 2) Set values and indices
			values = double(count.total);
			indices = integer(count.total);
			index = 0L;
			for (j in 1:ncols)
			{
				nz = which(R[, j] != 0);
				lnz = length(nz);
				values[index + 1:lnz] = R[nz, j];
				indices[index + 1:lnz] = nz - 1L; # must be zero-based
				index = index + lnz;
			}
			
			# Create sparse matrix
			retval = glmmGS.CreateSparseMatrix(values, indices, counts);
		}
	}
	else if (length(ls) == 3L)
	{
		# Retrieve data
		values = ls[[1]];
		indices = ls[[2]];
		counts = ls[[3]];
		
		# Create sparse matrix
		retval = glmmGS.CreateSparseMatrix(values, indices, counts);
	}
	else
	{
		stop("Invalid argument list");
	}
	
	return(retval);
}

# Construct covariance models
glmmGS.CovarianceModel = function(type, ...)
{
	V = NULL;
	ls = list(...);
	if (type == "identity")
	{
		V = list(type = "IdentityCovarianceModel");
	}
	else if (type == "precision")
	{
		V = NULL;
		R = ls[[1]];
		if (is.null(R) == TRUE)
		{
			stop("Invalid precision matrix");
		}
		else if ((is.matrix(R) == TRUE) && (nrow(R) == ncol(R)))
		{
			V = list(type = "PrecisionModel", R = R);
		}
		else if (is.null(attr(R, "sparse.matrix", TRUE)) == FALSE)
		{
			V = list(type = "SparsePrecisionModel", R = R);
		}
		else
		{
			stop("Invalid precision matrix");
		}
	}
	else
	{
		stop("Unsupported covariance model");
	}
	return(V);
}

# Fit
glmmGS.Fit = function(control)
{
	# Fit model
	glmmGSAPI.Fit(control$reltol, control$abstol, control$maxit);
	
	# Retrieve results
	fixed.effects = glmmGSAPI.GetFixedEffectsCoefficients();
	random.effects = glmmGSAPI.GetRandomEffectsCoefficients();
	covariance.components = glmmGSAPI.GetCovarianceComponents();
	iterations = glmmGSAPI.GetIterations();
	return(list(fixed.effects = fixed.effects,
		random.effects = random.effects,
		covariance.components = covariance.components,
		iterations = iterations));
}

# Get variable
glmmGS.GetVariable = function(name, data, env)
{
	x = NULL;
	if (!is.null(data))
		x = data[[name]];
	if (is.null(x))
		x = get(name, envir = env);
	if (is.null(x))
		stop(paste("Variable '", name, "' not found", sep = ""));
	return(x);
}

# Add covariance model to model
glmmGS.AddCovarianceModel = function(block, covariance.models, env)
{
	pos = 1;
	pos = glmmGSParser.Find(block, "~", pos) + 1;
	token = glmmGSParser.GetToken(block, "", pos);
	varname = glmmGSParser.Trim(token$text);
	covariance_model = glmmGS.GetVariable(varname, covariance.models, env);
	
	if (covariance_model$type == "IdentityCovarianceModel")
	{
		glmmGSAPI.AddIdentityCovarianceModel();
	}
	else if (covariance_model$type == "PrecisionModel")
	{
		glmmGSAPI.AddPrecisionModel(covariance_model$R);
	}
	else if (covariance_model$type == "SparsePrecisionModel")
	{
		glmmGSAPI.AddSparsePrecisionModel(covariance_model$R);
	}
	else
	{
		stop("Invalid covariance model");
	}
}

# Add variables to block
glmmGS.AddPredictors = function(block, data, env)
{
	token = glmmGSParser.GetToken(block, "(\\|)|(~)", 1);
	block = token$text;
	pos = 1;
	repeat
	{
		pos = glmmGSParser.SkipWhites(block, pos);
		token = glmmGSParser.GetToken(block, "\\+", pos);
		varname = glmmGSParser.Trim(token$text);
		if (varname == "1")
		{
			glmmGSAPI.AddIntercept();
		}
		else
		{
			glmmGSAPI.AddCovariate(glmmGS.GetVariable(varname, data, env));
		}
		if (token$pos == -1)
			break;
		pos = token$pos + 1;
	}
}

# Get factor data
glmmGS.GetFactor = function(block, data, env)
{
	pos = 1;
	pos = glmmGSParser.Find(block, "\\|", pos) + 1;
	token = glmmGSParser.GetToken(block, "~", pos);
	varname = glmmGSParser.Trim(token$text);
	return(glmmGS.GetVariable(varname, data, env));
}

# Add offset to model
glmmGS.AddOffset = function(offset, data, env)
{
	# Get response
	varname = glmmGSParser.Trim(offset$text);
	offset = glmmGS.GetVariable(varname, data, env);
	
	# Add offset
	glmmGSAPI.AddOffset(offset)
}
	
# Add predictor block to model
glmmGS.AddBlock = function(block, data, covariance.models, env)
{
	if (attr(block, "effects") == "fixed")
	{
		# Fixed effects
		glmmGSAPI.BeginFixedEffects();
		if (attr(block, "type") == "global")
		{
			# Global block
			glmmGSAPI.BeginGlobalBlock();
			glmmGS.AddPredictors(block$text, data, env);
			glmmGSAPI.EndGlobalBlock();			
		}
		else if (attr(block, "type") == "stratified")
		{
			# Stratified block
			glmmGSAPI.BeginStratifiedBlock(glmmGS.GetFactor(block$text, data, env));
			glmmGS.AddPredictors(block$text, data, env);
			glmmGSAPI.EndStratifiedBlock();			
		}
		else
		{
			stop("Invalid fixed-effects block type");
		}
		glmmGSAPI.EndFixedEffects();
	}
	else if (attr(block, "effects") == "random")
	{
		# Random effects
		glmmGSAPI.BeginRandomEffects();
		if (attr(block, "type") == "global")
		{
			# Global block
			glmmGSAPI.BeginGlobalBlock();
			glmmGS.AddPredictors(block$text, data, env);
			glmmGS.AddCovarianceModel(block$text, covariance.models, env);
			glmmGSAPI.EndGlobalBlock();			
		}
		else if (attr(block, "type") == "stratified")
		{
			# Stratified block
			glmmGSAPI.BeginStratifiedBlock(glmmGS.GetFactor(block$text, data, env));
			glmmGS.AddPredictors(block$text, data, env);
			glmmGS.AddCovarianceModel(block$text, covariance.models, env);
			glmmGSAPI.EndStratifiedBlock();			
		}
		else
		{
			stop("Invalid random-effects block type");
		}
		glmmGSAPI.EndRandomEffects();		
	}
	else
	{
		stop("Invalid 'effects' attribute");
	}
}

# Add response to model
glmmGS.AddResponse = function(response, family, data, env)
{
	# Response
	glmmGSAPI.BeginResponse(family);
	if (family == "binomial")
	{
		# Get response and counts
		pos = 1;
		pos = glmmGSParser.SkipWhites(response, pos);
		pos = glmmGSParser.Match(response, "\\(", pos);
		
		pos = glmmGSParser.SkipWhites(response, pos);
		token = glmmGSParser.GetToken(response, "\\|", pos);
		pos = token$pos + 1;
		y = glmmGS.GetVariable(glmmGSParser.Trim(token$text), data, env);
		
		pos = glmmGSParser.SkipWhites(response, pos);
		token = glmmGSParser.GetToken(response, "\\)", pos);
		pos = token$pos + 1;
		counts = glmmGS.GetVariable(glmmGSParser.Trim(token$text), data, env);
		
		pos = glmmGSParser.SkipWhites(response, pos);
		token = glmmGSParser.GetToken(response, ".", pos);
		if (token$pos != -1)
			stop("Invalid response formula");
		
		# Add response and counts
		glmmGSAPI.AddResponse(y);
		glmmGSAPI.AddCounts(counts);
	}
	else if (family == "poisson")
	{
		# Get response
		varname = glmmGSParser.Trim(response);
		y = glmmGS.GetVariable(varname, data, env);
		
		# Add response
		glmmGSAPI.AddResponse(y);
	}
	else
	{
		stop("Unsupported family");
	}
	glmmGSAPI.EndResponse();
}

# Main function
glmmGS = function(formula, family, data = NULL, covariance.models = NULL, control = glmmGS.Control())
{
	# Gets environment of calling function
	env = environment(formula);
	
	# Convert formula and family into text
	formula = as.character(formula);
	family = family()$family;
	
	# Split response and predictor string
	response = formula[2];
	predictor = formula[3];

	# Initialize state-machine
	glmmGSAPI.ForceEnd();
	glmmGSAPI.Begin();
	
	# Add response
	glmmGS.AddResponse(response, family, data, env);
	
	# Initialize position counter
	pos = 1;
	
	# Add offset
	offset = glmmGSParser.GetOffset(predictor, pos);
	if (is.null(offset) == FALSE)
	{
		pos = offset$pos;
		glmmGS.AddOffset(offset, data, env);
		pos = glmmGSParser.ParseSeparator(predictor, pos);
	}
	
	# Add predictor blocks
	while (pos != -1)
	{
		# Add block
		block = glmmGSParser.GetNextBlock(predictor, pos);
		pos = block$pos;
		glmmGS.AddBlock(block, data, covariance.models, env);

		# Parse separator and set position to next block 
		pos = glmmGSParser.ParseSeparator(predictor, pos);
	}
	
	# Fit model
	results = glmmGS.Fit(control);
	
	# Terminate state-machine
	glmmGSAPI.End();
	
	# Return results
	return(results);
}
