
## EventElement.r file
 ##
 ## Author: Jacob T. Ormerod
 ##   (c)2011 JTO Reliability
##


EventElement<-function(element_name,OpLine,EventID,FD,FP1,
                                    FP2,FP3,RD,RP1,RP2,RP3,Seed)  {		
		
		ElementDF<-data.frame("OpLine"=OpLine,"EventID"=EventID,"FD"=FD,"FP1"=FP1,
		  "FP2"=FP2,"FP3"=FP3,"RD"=RD,"RP1"=RP1,"RP2"=RP2,"RP3"=RP3,"Seed"=Seed) 
		
		row.names(ElementDF)<-element_name
		
		return(ElementDF)
	}	
