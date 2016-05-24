find.hazard <- function(t, h.val, h.ranges, type, mode){

	if(type=="increasing"){
		h	<-	find.h.up(t, h.val, h.ranges)
		}
	if(type=="decreasing"){
		h	<-	find.h.down(t, h.val, h.ranges)
		}
	if(type=="unimodal"){
		if(t<mode){
		h	<-	find.h.up(t, h.val, h.ranges)
		}
		if(t>mode){
		h	<-	find.h.down(t, h.val, h.ranges)
		}
		if(t==mode){
		h	<-	Inf
		}}
	if(type=="ushaped"){
		if(t<mode){
		h	<-	find.h.down(t, h.val, h.ranges)
		}
		if(t>mode){
		h	<-	find.h.up(t, h.val, h.ranges)
		}
		if(t==mode){
		h	<-	0
		}}

	return(h)

	}

