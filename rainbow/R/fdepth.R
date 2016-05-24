fdepth <- function(data, type = c("FM", "mode", "RP", "RPD", "radius"), trim = 0.25, alpha, weight)
{
  type = match.arg(type)
  if (type == "FM"){
      output <- depth.FM(data, trim = trim)
  }
  if (type == "mode"){
      output <- depth.mode(data, trim = trim)
  }
  if (type == "RP"){
      output <- depth.RP(data, trim = trim)
  }
  if (type == "RPD"){
      output <- depth.RPD(data, trim = trim)
  }
  if (type == "radius"){
  	  output <- depth.radius(data, alpha, beta = trim, weight)
  }
  return(structure(list(data = data, output = output), class = "fdepth"))
}

