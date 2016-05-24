`PLL.CoxBoost` <-
function(object, newdata, newtime, newstatus, complexity, ...){
   if(is.list(complexity)){
      PLL <- predict(object, type="logplik", newdata=newdata,
         newtime=newtime, newstatus=newstatus, at.step=complexity$stepno)
      } else {
      PLL <- predict(object, type="logplik", newdata=newdata,
         newtime=newtime, newstatus=newstatus, at.step=complexity)
      }
   PLL
}

