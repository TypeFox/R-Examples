`rst` <-
function (x,phi=getOption("latitude")) 
{
	phi = phi/180*pi;
	
	x = as.eqc(x);

	q = -tan(phi)*tan(x[,2]/180*pi);

	qq = q;
	qq[q>1] = NA;
	qq[q< -1] = NA;

	rise = 24 - acos(qq)*3.819718634205488 + x[,1];
	rise = rise %% 24;
	rise[q>1] = -Inf;
	rise[q< -1] = Inf;

	set = acos(qq)*3.819718634205488 + x[,1];
	set = set %% 24;
	set[q>1] = -Inf;
	set[q< -1] = Inf;

	transit = x[,1];
	transit[q>1] = -Inf;

	aa = sin(x[,2]/180*pi)/cos(phi);
	aa[aa>1] = NA;
	aa[aa< -1] = NA;

	arise = acos(aa)*180/pi;
	aset = 360 - arise;

	arise[q>1] = -Inf;
	arise[q< -1] = Inf;
	aset[q>1] = -Inf;
	aset[q< -1] = Inf;

	res = data.frame(rise=as.vector(rise),transit=as.vector(transit),
		set=as.vector(set),arise=as.vector(arise),aset=as.vector(aset));
	rownames(res) = rownames(x);

	class(res) = c("rst","data.frame");
	class(res$rise) = c("lst","time");
	class(res$set) = c("lst","time");
	class(res$transit) = c("lst","time");
	class(res$arise) = "dms";
	class(res$aset) = "dms";


	return(res);

}

