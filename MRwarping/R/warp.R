warp <-
function(A,Lambda,R1,R2,x)
{
Wx <- x
for (i in 1:length(A))
    	{
	warp <- comp(A[i],Lambda[i],R1[i],R2[i],Wx)
	Wx <- warp$Wx
	}
return(Wx)
}

