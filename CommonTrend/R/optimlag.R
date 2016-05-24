optimlag <-
function ( data, k=10,ecdet="none")
{      
  result = matrix(NA,k,1)
  for(j in 2:k)
      { cotest = ca.jo(data, type = "trace", ecdet =ecdet, K=j, spec="transitory")      
        result[ j,]=cotest@teststat[2]
       }
  final=list (   Olag.value=max(result[2:k]),
                 Olag= which(  result==max(result[2:k])  ),
             list.lags=result  
             )
  class(final)<-"optimlag"
  final
}