randomnumGenerator <-
function(randomnumRange,lengthRandomnum,howManyRands)
   {
# this function generates exactly "howManyRands" different integer random numbers
# given the list of the required numbers.

# howManyRands defines how many different integer numbers you want in the defined range

# randomnumRange can be smt like this:
# randomnumRange=1:487
# or
# randomnumRange=6 12 36 67 87 146 267
# therefore length of the randonnumRange will be
# in the first case 487
# in the second case 7 defining the index of the randomnumRange to be swapped.


  	for (i in 1:(howManyRands))
   	{
		integer=sample(lengthRandomnum-1,1)
		tmp=randomnumRange[integer]
		randomnumRange[integer]=randomnumRange[i]
		randomnumRange[i]=tmp
   	}
	return(randomnumRange[1:howManyRands])
   
   }

