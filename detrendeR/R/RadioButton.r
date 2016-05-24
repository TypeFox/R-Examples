RadioButton = function (FRAME, variable= NULL, BUTTON=c("b.r1", "b.r2"), VALUE=c(TRUE,FALSE)){
 BUTTON<-as.vector(BUTTON)
 for (i in 1:length(BUTTON)){
 tkpack(tkradiobutton(FRAME, text = BUTTON[i], value = VALUE[i], variable = variable), anchor = "w")
 }
 }
