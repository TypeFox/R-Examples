# R low-level interface with GlmmAPI functions

# Error buffer
error_buffer.size = as.integer(1024);
error.buffer = sprintf("%1024s", "");

# Get error occurred in last API call
glmmGSAPI.GetLastError = function()
{

	# Get error message
	size = nchar(error.buffer[1]);
	out = .C("GlmmGSRAPI_GetLastError", error.buffer[1], size, PACKAGE = "glmmGS");
	
	# Trim trailing spaces
	msg = sub(" +$", "", out[[1]]);
	
	# Check error
	if (msg != "")
		stop(msg);
}

# Debug function: only to be used to debug the application
glmmGSAPI.SetOutputFile = function(filename)
{
	.C("GlmmGSRAPI_SetOutputFile", as.character(filename), PACKAGE = "glmmGS");
}

# Begin glmmGSAPI
glmmGSAPI.Begin = function()
{
	.C("GlmmGSRAPI_Begin", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End glmmGSAPI
glmmGSAPI.End = function()
{
	.C("GlmmGSRAPI_End", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End glmmGSAPI
glmmGSAPI.ForceEnd = function()
{
	.C("GlmmGSRAPI_ForceEnd", PACKAGE = "glmmGS");
}

# Begin response
glmmGSAPI.BeginResponse = function(family)
{
	if (is.character(family) == FALSE)
		stop("Invalid type");
	if (length(family) != 1)
		stop("Invalid size");
	size = as.integer(nchar(family[1]));
	.C("GlmmGSRAPI_BeginResponse", family[1], size, PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End response
glmmGSAPI.EndResponse = function()
{
	.C("GlmmGSRAPI_EndResponse", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Begin fixed-effects
glmmGSAPI.BeginFixedEffects = function()
{
	.C("GlmmGSRAPI_BeginFixedEffects", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End fixed-effects
glmmGSAPI.EndFixedEffects = function()
{
	.C("GlmmGSRAPI_EndFixedEffects", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Begin random-effects
glmmGSAPI.BeginRandomEffects = function()
{
	.C("GlmmGSRAPI_BeginRandomEffects", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End random-effects
glmmGSAPI.EndRandomEffects = function()
{
	.C("GlmmGSRAPI_EndRandomEffects", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Begin global block
glmmGSAPI.BeginGlobalBlock = function()
{
	.C("GlmmGSRAPI_BeginGlobalBlock", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# End global block
glmmGSAPI.EndGlobalBlock = function()
{
	.C("GlmmGSRAPI_EndGlobalBlock", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Begin stratified block
glmmGSAPI.BeginStratifiedBlock = function(levels)
{
	size = as.integer(length(levels));
	if (is.integer(levels))
	{
		.C("GlmmGSRAPI_BeginStratifiedBlock", levels, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid factor type");
	}
}

# End stratified block
glmmGSAPI.EndStratifiedBlock = function()
{
	.C("GlmmGSRAPI_EndStratifiedBlock", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Add response vector
glmmGSAPI.AddResponse = function(values)
{
	size = as.integer(length(values));
	if (is.integer(values))
	{
		.C("GlmmGSRAPI_AddResponseInt", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid type");
	}
}

# Add vector of counts (for binomial response)
glmmGSAPI.AddCounts = function(values)
{
	size = as.integer(length(values));
	if (is.integer(values))
	{
		.C("GlmmGSRAPI_AddCountsInt", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid type");
	}
}

# Add offset
glmmGSAPI.AddOffset = function(values)
{
	size = as.integer(length(values));
	if (is.integer(values))
	{
		.C("GlmmGSRAPI_AddOffsetInt", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else if (is.double(values))
	{
		.C("GlmmGSRAPI_AddOffsetDbl", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid type");
	}
}

# Add intercept
glmmGSAPI.AddIntercept = function()
{
	.C("GlmmGSRAPI_AddIntercept", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Add a covariate
glmmGSAPI.AddCovariate = function(values)
{
	if (is.vector(values))
	{
		size = as.integer(length(values));
		if (is.integer(values))
		{
			.C("GlmmGSRAPI_AddCovariateInt", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
			glmmGSAPI.GetLastError();
		}
		else if (is.double(values))
		{
			.C("GlmmGSRAPI_AddCovariateDbl", values, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
			glmmGSAPI.GetLastError();
		}
		else
		{
			stop("Invalid type");
		}
	}
	else if (is.matrix(values))
	{
		dimensions = dim(values);
		if (is.integer(values))
		{
			.C("GlmmGSRAPI_AddCovariatesInt", values, dimensions, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
			glmmGSAPI.GetLastError();
		}
		else if (is.double(values))
		{
			.C("GlmmGSRAPI_AddCovariatesDbl", values, dimensions, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
			glmmGSAPI.GetLastError();
		}
		else
		{
			stop("Invalid type");
		}
	}
	else
	{
		stop("Invalid type");
	}
}

# Add identity covariance model
glmmGSAPI.AddIdentityCovarianceModel = function()
{
	.C("GlmmGSRAPI_AddIdentityCovarianceModel", PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Add precision model
glmmGSAPI.AddPrecisionModel = function(R)
{
	if ((is.matrix(R) == TRUE) && (is.double(R) == TRUE))
	{
		dim = dim(R);
		.C("GlmmGSRAPI_AddPrecisionModel", R, dim[1], dim[2], DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid type");
	}
}

# Add sparse-precision model
glmmGSAPI.AddSparsePrecisionModel = function(R)
{
	if (is.null(attr(R, "sparse.matrix", TRUE)) == FALSE)
	{
		ncols = as.integer(length(R$counts)) - 1L;
		.C("GlmmGSRAPI_AddSparsePrecisionModel", R$values, R$indices, R$counts, ncols, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	else
	{
		stop("Invalid type");
	}
}

# Fit current model
glmmGSAPI.Fit = function(relative.tolerance, absolute.tolerance, maxiter)
{
	.C("GlmmGSRAPI_Fit", as.double(relative.tolerance), as.double(absolute.tolerance), as.integer(maxiter), PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
}

# Get estimated fixed-effects coefficients
glmmGSAPI.GetFixedEffectsCoefficients = function()
{
	size = integer(1);
	.C("GlmmGSRAPI_GetFixedEffectsSize", size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
	est = double(size);
	err = double(size);
	if (size > 0)
	{
		.C("GlmmGSRAPI_GetFixedEffectsEstimates", est, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
		.C("GlmmGSRAPI_GetFixedEffectsErrors", err, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	return(list(estimates = est, standard.errors = err));
}

# Get estimated random-effects coefficients
glmmGSAPI.GetRandomEffectsCoefficients = function()
{
	size = integer(1);
	.C("GlmmGSRAPI_GetRandomEffectsSize", size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
	est = double(size);
	err = double(size);
	if (size > 0)
	{
		.C("GlmmGSRAPI_GetRandomEffectsEstimates", est, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
		.C("GlmmGSRAPI_GetRandomEffectsErrors", err, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	return(list(estimates = est, standard.errors = err));
}

# Get estimated covariance components
glmmGSAPI.GetCovarianceComponents = function()
{
	size = integer(1);
	.C("GlmmGSRAPI_GetCovarianceComponentsSize", size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
	est = double(size);
	err = double(size);
	if (size > 0)
	{
		.C("GlmmGSRAPI_GetCovarianceComponentsEstimates", est, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
		.C("GlmmGSRAPI_GetCovarianceComponentsErrors", err, size, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
		glmmGSAPI.GetLastError();
	}
	return(list(estimates = est, standard.errors = err));
}

# Get number of iterations
glmmGSAPI.GetIterations = function()
{
	iterations = integer(1);
	.C("GlmmGSRAPI_GetIterations", iterations, DUP = FALSE, NAOK = FALSE, PACKAGE = "glmmGS");
	glmmGSAPI.GetLastError();
	return(iterations);
}

