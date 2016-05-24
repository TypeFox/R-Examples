refine.extraction <-
function(m.data, c.model, s.model)
{
	
	neg.sign.flag <- any(sign(normalize(m.data[which.max(c.model),])-normalize(s.model))==(-1))
	its <- 1
	while(neg.sign.flag)
	{
		substr <- normalize(m.data[which.max(c.model),]) - normalize(s.model)
		substr[substr>0] <- 0
		s.model <- normalize(s.model) + substr
		neg.sign.flag <- any(sign(normalize(m.data[which.max(c.model),])-normalize(s.model))==(-1))
		its <- its + 1
		if(its>10) break
	}
	s.model[s.model<0] <- 0
		
	s.model	
}
