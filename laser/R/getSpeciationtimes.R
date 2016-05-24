`getSpeciationtimes` <-
function(BranchingTimes)
{
	Tmax <- BranchingTimes[1];
  	SpeciationTimes <- Tmax - BranchingTimes;
  	x <- c(SpeciationTimes[1], SpeciationTimes);
  	x;
}

