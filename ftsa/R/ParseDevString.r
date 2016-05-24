ParseDevString = function (method)
{
	if (method[1] == "mad")  return (0)
	if (method[1] == "sd")  return (1)
	if (method[1] == "Qn" | method[1] == "qn" )  return (2)
	return (1)
}
