variance <-
function(x,y,o.nsims,nsims, qhat){
   
   #get another 200 MC at n yielding 200 teststats
   teststat2<-simulatedata(x,y,o.nsims)
   den2 = ((density(teststat2, from=qhat, to=qhat)$y[1])^2)*nsims
   F.fun =  ecdf(teststat2)
   num2 = F.fun(qhat)*(1-F.fun(qhat))
   var.est = num2/den2
   returnValue = list(teststat2, var.est)
   return(returnValue)
 }
