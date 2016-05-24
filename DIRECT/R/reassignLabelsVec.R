reassignLabelsVec <-
function (x)
{
	labels = unique (x)
	x.new = x

	nlabels = length (labels)
	for (i in 1:nlabels)
		x.new[which (x==labels[i])] = i

	return (x.new)
}

