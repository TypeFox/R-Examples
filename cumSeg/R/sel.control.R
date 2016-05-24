sel.control <-
function (display=FALSE, type=c("bic","mdl","rss"), S=1, Cn="log(log(n))",
      alg=c("stepwise","lasso"), edf.psi=TRUE){
#type=criterio per selezionare il n.psi
#alg=algoritmo da passare a lars().
#Cn=stringa (eventualmente in funzione di n) relativa al criterio bic da usare..
#S: valore soglia (sensato solo per type="rss"
    type <- match.arg(type)
    alg <- match.arg(alg)
    list(display=display, type=type, S=S, Cn=Cn, edf.psi=edf.psi, alg=alg)
    }

