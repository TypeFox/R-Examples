LSCV.risk.single <- function(h,cases,controls,method,res,edge,inside,ERRV=1000000){
	if(h<=0) return(NA)
	temp.case.pts <- density(cases,sigma=h,edge=edge,dimyx=c(res,res),at="points",leaveoneout=T)/cases$n
	temp.case <- density(cases,sigma=h,edge=edge,dimyx=c(res,res))
	temp.con.pts <- density(controls,sigma=h,edge=edge,dimyx=c(res,res),at="points",leaveoneout=T)/controls$n
	temp.con <- density(controls,sigma=h,edge=edge,dimyx=c(res,res))
	
	ca <- (diff(cases$window$xrange)/res)*(diff(cases$window$xrange)/res)
	
	if(any(is.na(temp.case$v[t(inside)]))||any(is.na(temp.con$v[t(inside)]))) return(NA)
	if(any(temp.case$v[t(inside)]<=0)||any(temp.con$v[t(inside)]<=0)) return(ERRV)
	
	W.im <- as.mask(cases$window)
	casedens.at.cons.indices <- nearest.raster.point(controls$x,controls$y,W.im)
	condens.at.cases.indices <- nearest.raster.point(cases$x,cases$y,W.im)
	
	caseatcon <- diag(temp.case$v[casedens.at.cons.indices$row,casedens.at.cons.indices$col]/cases$n)
	conatcase <- diag(temp.con$v[condens.at.cases.indices$row,condens.at.cases.indices$col]/controls$n)
	
	if(method=="kelsall-diggle"){
		CV <- 2*mean(log(caseatcon/temp.con.pts)/temp.con.pts)-2*mean(log(temp.case.pts/conatcase)/temp.case.pts)-sum((log(temp.case$v/temp.con$v)^2)*ca,na.rm=T)
	} else if(method=="hazelton"){
		CV <- mean((caseatcon/temp.con.pts)^2)-2*mean(temp.case.pts/conatcase)
	}
	return(CV)
}

