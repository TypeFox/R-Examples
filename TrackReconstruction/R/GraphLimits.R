#GraphLimits1=function(infile)
#{
#	maxlat=max(infile$Latitude)
#	minlat=min(infile$Latitude)
#	maxlong=max(infile$Longitude)
#	minlong=min(infile$Longitude)
#	
#	ratio=(maxlat-minlat)/(maxlong-minlong)
#	if(ratio<0.5)
#	{
#		miny=minlat+(maxlat-minlat)/2-(maxlong-minlong)/4-0.025*(maxlong-minlong);
#		maxy=minlat+(maxlat-minlat)/2+(maxlong-minlong)/4+0.025*(maxlong-minlong);
#		minx=minlong-0.05*(maxlong-minlong);
#		maxx=maxlong+0.05*(maxlong-minlong);
#	}else{
#		miny=minlat-0.025*(maxlat-minlat);
#		maxy=maxlat+0.025*(maxlat-minlat);
#		minx=minlong+(maxlong-minlong)/2-(maxlat-minlat)-0.05*(maxlat-minlat);
#		maxx=maxlong-(maxlong-minlong)/2+(maxlat-minlat)+0.05*(maxlat-minlat);
#	}
#	return=list(miny=miny,maxy=maxy,minx=minx,maxx=maxx)
#}



GraphLimits=function(infile)
{
	maxlat=max(infile$Latitude)
	minlat=min(infile$Latitude)
	maxlong=max(infile$Longitude)
	minlong=min(infile$Longitude)
	
	ratio=cos((maxlat+minlat)/2/360*2*pi)
	if(((maxlong-minlong)*ratio) > (maxlat-minlat))
	{
		miny=minlat-((maxlong-minlong)*ratio)/2-0.025*(maxlong-minlong)*ratio;
		maxy=maxlat+((maxlong-minlong)*ratio)/2+0.025*(maxlong-minlong)*ratio;
		minx=minlong-0.025*(maxlong-minlong);
		maxx=maxlong+0.025*(maxlong-minlong);
	}else{
		miny=minlat-0.025*(maxlat-minlat);
		maxy=maxlat+0.025*(maxlat-minlat);
		minx=minlong-((maxlat-minlat)-((maxlong-minlong)*ratio))/2-0.025*(maxlat-minlat)/ratio;
		maxx=maxlong+((maxlat-minlat)-((maxlong-minlong)*ratio))/2-0.025*(maxlat-minlat)/ratio;
	}
	return=list(miny=miny,maxy=maxy,minx=minx,maxx=maxx)
}
