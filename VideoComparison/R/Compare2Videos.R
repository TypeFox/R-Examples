Compare2Videos<-function(nv1,nv2,stp=10,fsc=0){
	mm1<-ExtractMotion(nv1)  
	mm2<-ExtractMotion(nv2)
#
	rs<-NULL
	if ( length(mm1) > 0 & length(mm2) > 0 ) {
		res2<-VideoComparison(mm1,mm2,stp)
		rr<-0
		if ( res2$pos1 > 0 & res2$pos2 > 0 ) {
			if ( fsc==0) fsc<-res2$sc
			img1<-ExtractImgPos(nv1)
			img2<-ExtractImgPos(nv2)
#
			p1 <- res2$pos1
			p2 <- res2$pos2
			if ( max(as.numeric(img1)) > max(as.numeric(img2))) {
				p1 <- res2$pos2
				p2 <- res2$pos1
			}
			idx2<-((as.numeric(img2)-p2) > 0 ) &  
  				((as.numeric(img2)-p2) < res2$lngth)
			idx1<-((as.numeric(img1)-p1) > 0 ) &  
  				((as.numeric(img1)-p1) < res2$lngth)
			rimg1<-img1[idx1]  
			rimg2<-img2[idx2] 
			lh1<-apply(as.matrix(rimg1),1,ExtractImgHash,nv1)
			lh2<-apply(as.matrix(rimg2),1,ExtractImgHash,nv2)
#
			rr<-VideoMatch(lh2,lh1,fsc)
		} 
		rs<-list(likelihood=rr,msc=res2$sc,mp1=floor(res2$pos1),
				mp2=floor(res2$pos2),lngth=floor(res2$lngth))
	}
	return(rs)
}
