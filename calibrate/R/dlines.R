"dlines" <-
function(SetA,SetB,lin="dotted") {
# 
# Function DLINES connects the rows of SetA to the rows of SetB with lines
# 
# Jan Graffelman
# Universitat Politecnica de Catalunya
# January 2004
#
np<-nrow(SetA)
for(i in 1:np) lines(rbind(SetA[i,1:2],SetB[i,1:2]),lty=lin) 
return(NULL)
}

