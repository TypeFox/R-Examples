varMarkdown <- function(effect, mod, modsum, statistic='t', B=1, CI=B, pe=2)
{
  results =''

  if (is.numeric(pe)){
    pestat = (abs(coef(modsum)[effect,'t'])^2*1)/(abs(coef(modsum)[effect,'t'])^2*1 + mod$df.residual)
    results = paste0(results, '*p$\\eta^2$* = ', sprintf(paste0('%.',pe,'f'), pestat), ', ')
  }  
  
  if (is.numeric(B)){
    results = paste0(results,'*b* = ', sprintf(paste0('%.',B,'f'),coef(mod)[effect]), ', ')
  }
  
  if (is.numeric(CI)){
    results = paste0(results,'*95% CI(B)* = [', sprintf(paste0('%.',CI,'f'),confint(mod)[effect,1]), ',', sprintf(paste0('%.',CI,'f'),confint(mod)[effect,2]), '], ')
  }
  
  #add t
  if (statistic == 't'){
    results = paste0(results,'*t*(', mod$df.residual, ') = ', sprintf('%.2f',abs(coef(modsum)[effect,'t'])), ', ')
  }  
  
  #add p
  if (round(coef(modsum)[effect,'Pr(>|t|)'],3) < 0.001){
    results = paste0(results, '*p* <  .001')
  }else{
    results = paste0(results, '*p* =  ', sprintf('%.3f',coef(modsum)[effect,'Pr(>|t|)']))
  }
  
    
  
  
  return(results)
}