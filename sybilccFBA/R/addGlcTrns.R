addGlcTrns<-function(prob,mod2){
Si=0.2# units?
Hxt1	=41*Si/(Si+107)#Kcat*Si/(Si+Km)
Hxt2	= 16.1*Si/(Si+2.9)
Hxt3	= 18.5*Si/(Si+29)
Hxt4	=12*Si/(Si+6.2)
Hxt5	= 14*Si/(Si+10)
Hxt6	=11.4*Si/(Si+1.5)
Hxt7	=11.7*Si/(Si+1.3)
Gal2	=17.5*Si/(Si+1.5)
#Gal2 in existence of Galactose only
colid=getNumCols(lp = problem(prob))+1
trnsCol=NULL
# for(i in (1:7)){
    # trnsCol=rbind(trnsCol,cbind(trns=i,Col=colid))
	# addCols(lp = problem(prob),1)
	# changeColsBnds(lp = problem(prob),colid,lb=0,ub=1000)
	# colid=colid+1;
	# }
	#Hxt4 YHR092C
	#Hxt1 YHR094C
	#Hxt2 YMR011W
		#( YHR092C  or  YLR081W  or  YOL156W  or  YDR536W  or  YHR094C  or  YEL069C  or  YDL245C  or  YJR158W  or  YFL011W  or  YNR072W  or  YMR011W  or  YDR345C  or  YHR096C  or  YDR343C  or  YDR342C  or  YJL214W  or  YJL219W )			
		#( Hxt4 ) or ( Gal2 ) or ( Hxt11 ) or ( Stl1 ) or ( Hxt1 ) or ( Hxt13 ) or ( Hxt15 ) or ( Hxt16 ) or ( Hxt10 ) or ( Hxt17 ) or ( Hxt2 ) or ( Hxt3 ) or ( Hxt5 ) or ( Hxt6 ) or ( Hxt7 ) or ( Hxt8 ) or ( Hxt9 )"
		#Add constraint
rowind=getNumRows(lp = problem(prob))+1
#glcRxn=which(react_id(mod2)=='R_EX_glc_e__b')
glcRxn=which(react_id(mod2)=="R_GLCt1")
addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(trnsCol[1,"Col"],trnsCol[2,"Col"],trnsCol[3,"Col"],trnsCol[4,"Col"],trnsCol[5,"Col"],trnsCol[6,"Col"],
			  trnsCol[7,"Col"],glcRxn)),
              nzval = list(c(-Hxt1,-Hxt2,-Hxt3,-Hxt4,-Hxt5,-Hxt6,-Hxt7,1))
			  ,rnames = "glcTrns"
			 )
#Add to crowding constraint (same budget), Molecular weights required?
return(prob)
}
