local({
## Vorbereiten
require(klausuR)

## Berechne
klsr.data.obj <- klausur.data(answ=antworten, corr=answers, marks=notenschluessel)
## Drucke Ergebnisse
rk.header("klausuR: Prepare Test Data")
assign("klsr.data.obj", klsr.data.obj, envir=.GlobalEnv)
})
