 colors2data<-function(dendat,pcf,lst,paletti=NULL,clusterlevel=NULL,nodes=NULL,
type="regular")
{
if (type=="regular")
return( colorsofdata(dendat,pcf,lst,paletti=paletti,clusterlevel=clusterlevel,nodes=nodes) )
else
return( colorsofdata.adagrid(dendat,pcf,lst,paletti=paletti,clusterlevel=clusterlevel,nodes=nodes) )

}

