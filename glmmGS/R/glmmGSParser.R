# Trim
glmmGSParser.Trim = function(string)
{
	# Trim trailing white characters
	return(gsub("(^ +)|( +$)", "", string));
}

# Find
glmmGSParser.Find = function(string, pattern, start)
{
	pos = -1;
	match_len = 0;
	positions = gregexpr(pattern, string)[[1]];
	indices = which(positions >= start);
	n = length(indices)
	if (n > 0)
	{
		for (i in 1:n)
		{
			match_len = attr(positions, "match.length")[indices[i]];
			if (match_len > 0)
			{
				pos = positions[indices[i]];
				break;
			}
		}
	}

	return(pos);
}

# Skip
glmmGSParser.SkipWhites = function(string, start)
{
	positions = gregexpr(" ", string)[[1]];
	indices = which(positions >= start);
	pos = start;
	n = length(indices)
	if (n > 0)
	{
		for (i in 1:n)
		{
			if (positions[indices[i]] == pos)
			{
				pos = positions[indices[i]] + 1;
			}
			else
			{
				break;
			}
		}
	}
	return(pos);
}

# Match
glmmGSParser.Match = function(string, pattern, start)
{
	pos = NULL;
	positions = gregexpr(pattern, string)[[1]];
	index = which(positions == start);
	if (length(index) == 0)
		stop("Wrong format");
	pos = start + attr(positions, "match.length")[index];
	return(pos);
}

# Get token
glmmGSParser.GetToken = function(string, sep, start)
{
	text = NULL;
	stop = NULL;
	pos = glmmGSParser.Find(string, sep, start);
	if (pos < 0)
	{
		stop = nchar(string);
	}
	else
	{
		stop = pos - 1;
	}
	text = substr(string, start, stop);
	return(list(text = text, pos = pos));
}

# glmmGSParser.GetOffset
glmmGSParser.GetOffset = function(predictor, pos)
{
	pos = glmmGSParser.SkipWhites(predictor, pos);
	token = glmmGSParser.GetToken(predictor, "\\(", pos);
	if (token$text != "offset")
		return(NULL)
	
	# token$text == offset
	pos = token$pos + 1;
	token = glmmGSParser.GetToken(predictor, "\\)", pos);
	offset = list(text = token$text, pos = token$pos + 1);
	return(offset);
}

# glmmGSParser.GetNextBlock
glmmGSParser.GetNextBlock = function(predictor, pos)
{
	pos = glmmGSParser.SkipWhites(predictor, pos);
	pos = glmmGSParser.Match(predictor, "\\(", pos);
	token = glmmGSParser.GetToken(predictor, "\\)", pos);
	block = list(text = token$text, pos = token$pos + 1);
	
	if (glmmGSParser.Find(block, "~", 1) > 0)
	{
		attr(block, "effects") = "random";
	}
	else
	{
		attr(block, "effects") = "fixed";
	}
	
	if (glmmGSParser.Find(block, "\\|", 1) > 0)
	{
		attr(block, "type") = "stratified";
	}
	else
	{
		attr(block, "type") = "global";
	}

	return(block);
}

glmmGSParser.ParseSeparator = function(predictor, pos)
{
	pos = glmmGSParser.SkipWhites(predictor, pos);
	token = glmmGSParser.GetToken(predictor, "\\(", pos);
	pos = token$pos;
	text = glmmGSParser.Trim(token$text);
	if (nchar(text) == 0)
	{
		pos = -1;
	}
	else if (text != "+")
	{
		stop("Wrong format");
	}
	return(pos);
}
