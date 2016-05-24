# 0.4-13:   Corretto bug in combin_length e ripristinata funzione respdim in C; corretto bug
#           in funzione outllier.replace.
# 0.4-12:   Funzione R respdim inserita anche all'interno di rav().
# 0.4-11:   Funzione respdim riscritta in R per problema di gestione della memoria; bisogna
#           verificare ed eventualmente correggere e ripristinare la vecchia funzione C.
# 0.4-10:   Aggiunta verifica per valore TSS: se zero, non c'e' variabilita' nei dati e nessun
#           calcolo puo' essere svolto.
# 0.4-9:    Aggiornato il calcolo per la standardizzazione dei residui in rav.resid. Funzione
#           'empty.grid' rinominata in 'rav.grid'.
# 0.4-8:    Aggiornata la definzione delle classi al nuovo standard della versione 3 di R. Aggiunte
#           le funzioni rav.single e rav2file.
# 0.4-7:    Corretto bug nella definizione di s.range per metodo L-BFGS-B. Corretto bug nella
#           definizione di N. Implementate le funzioni rav.AIC e rav.BIC. Aggiunto il soggetto 29
#           al dataset 'pasta'. Aggiornato codice per la stampa del messaggio di caricamento.
# 0.4-6:    Corretto bug nel cntrollo if di s.range. Adattata la dimensione della maschera per il
#           fissaggio dei parametri.
# 0.4-5:    Aggiunta una patch per correggere un bug di optim: se optim non converge e da' il
#           messaggio di errore 'ERROR: NO FEASIBLE SOLUTION', allora il valore RSS in output non
#           corrispondera' al vero RSS prodotto dai parametri.
# 0.4-4:    Riordinati gli argomenti di rav. Corretto conflitto tra variabili nella funzione pargen.
# 0.4-3:    Routine di minimizzazione di default modificata in BFGS. Funzioni fit.indexes e
#           rav.indexes rinominate in fit.indices e rav.indices. Funzione rav.sqresid accorpata alla
#           funzione rav.resid. Argomenti type.par e t.out sostituiti con t.par. Eliminato
#           arrotondamento decimale dei valori su datgen e rav.resid. Eliminata procedura di
#           uguagliamento dei pesi in DAM. Correzione sulla funzione che calcola il numero di pesi:
#           un peso non viene considerato nel computo. Eliminata funzione rav.cmd. Rimosso indice
#           chi-quadrato. Riportato in output il metodo di minimizzazione utilizzato.
# 0.4-2     Corretti i nomi e le sigle dei modelli e alcuni bug interni alla funzione rav.
#           Nelle funzioni rav.param, rav.resid e rav.fitted, l'argmento 'which' è stato rinominato
#           in 'whichModel'. Inoltre, l'argomento ora richide una stringa di carattere che specifica
#           la sigla del modello.
#           La funzione pargen ora manda in output valori NA per s0 e w0 se I0=F. Inoltre, attraverso
#           l'argomento type.par, si può specificare quale tipo di parametro si vuole generare, se
#           'w' oppure 't'. I pesi in output sono normalizzati.
#           La funzione datgen ora richiede che venga specificato, tramite l'argomento type.par,
#           se il vettore di parametri in input contiene i pesi in formato 'w' oppure 't'. Inoltre,
#           non è necessario che le venga specificato se lo stato iniziale deve essere presente o
#           assente, perche' lo capisce da sola analizzando il valore di w0: se questo e' NA, allora
#           non considera lo stato iniziale. N.B.: dato che comunque nella generazione degli R sono
#           richiesti dei valori per s0 e w0, questi verranno posti rispettivamente a 0 e a 1e-10.
#           Nella funzione rav, I0 è ora posto a F per default. Inoltre, se I0=F, allora i valori di
#           s0 e w0 saranno restituiti pari a NA.
#           La funzione outlier.remove è stata rinominata in outlier.replace, perché effettivamente
#           non elimina gli outlier ma li sostituisce con un valore. Ora l'utente può passare come
#           valore da sostituire agli NA anche una funzione che, applicata per colonna, calcolera'
#           da se' un valore da sistituire, diverso per ogni colonna della matrice. N.B.: l'argomento
#           'value' prima era denominato 'replace'.
#           Aggiustata la funzione rav.indexes. Come per datgen, non ha piu' bisogno che le venga
#           specificato se lo stato iniziale e' presente o assente, perche' e' in grado di capiro da
#           sola dal valore di w0 (se w0=NA, allora I0=F). Inoltre, grazie all'argomento type.par,
#           e' in grado di capire se i parametri che le vengono passati sono 'w' oppure 't', agendo
#           di conseguenza. Se w0=NA, il valore verra' di nascosto posto pari a 1e-10. Corretto un
#           bug che impediva che il titolo del modello fosse riportato in output.
#           Aggiunta la funzione rav.sqresid, che calcola le medie dei residui al quadrato.
# 0.4-1     Corretti bug della versione 0.4-0. Riunite le funzioni rav e rav.base. Eliminata la
#           funzione rav.sim: dato che i valori di start vengono identificati sequenzialmente, non
#           ha piu' senso stabilirli tramite simulazione Monte-Carlo. In rav, quando I0=F, i valori
#           s0 e w0 non vengono mostrati in output. Per rav, aggiunto l'argomento t.out, che se TRUE
#           mostra in output i t invece dei w.
# 0.4-0     Implementato l'argoritmo in versione 'esponenziale': l'ottimizzazione non avvienee piu'
#           sui parametri w ma sui t. Implementati modelli null, ESM, SAM.
# 0.3-0     Versione ottimizzata della 0.2, modificata parmean, aggiunto argomento sim
# 0.2-2     Corretto bug sul calcolo di TSS
# 0.2-1     Funzioni averaging, residuals, combin e calcolo dei parameteri implementate in C.
#           Aggiunta possibilita' di fissaggio dei pesi, uguagliamento dei 'pesi uguali'
#               entro un certo delta
#           Generale ottimizzazione del codice
# 0.2-0     Estensione a modelli con infiniti fattori e livelli
#           Generale ottimizzazione del codice
#           Eliminata provvisoriamente models_12 dalla libreria
#           Eliminate alcune funzioni importate nel namespace
#           Cambiati i nomi delle funzioni e delle classi
# 0.1-0     Riorganizzato il codice di models_3
# 0.0-12    Aggiunto l'help delle funzioni, Inserimento dataset sperimentazione
#               Verifica dei modelli add + mul per 4 fattori, Verificare ANOVA-RM
#               sui modelli averaging, Creazione dei file di help e del manuale
#               Inserimento file di esempio
# 0.0-11    Correzione per compatibilita con gcc4 e R2.6
# 0.0-10    Correzione di numerosi bug nelle funzioni di stima
# 0.0-09    Inserire funzione che date osservazioni e parametri,
#               calcola indici di fit (2007/07/21)
#               Inserire nativo Chi-Quadro negli indici di bonta
# 0.0-08    Predisposizione per modelli averaging a 4 fattori (2007/07/17)
#               Verifica di tutte le componenti da rendere pubbliche 
# 0.0-07    rAverage reso come pacchetto autonomo da mvUtils (2007/05/15)
# 0.0-06    Inserimento test Chi-Quadro su modelli Average
# 0.0-05    Organizzazione dei parametri, delle funzioni e delle classi
# 0.0-04    Verifica con R CMD CHECK su Windows XP / minGW / R 2.3.0 (2006/05/10)
#               Verifica con R CMD CHECK su Mac Os X - Tiger / R 2.3.0
#               Installazione del pacchetto rAverage
# 0.0-03    Implentazione autonoma dei modelli additivi e moltiplicativi (2006/04/25)
# 0.0-02    Selezione modelli secondo AIC e BIC congiuntamente
# 0.0-02    Funzioni principali scritte in C (2006/02/15)
# 0.0-01    Funziona con modelli 3x3 e 3x3x3x3 (2006/02/01)

.packageName <- "rAverage"
.onAttach <- function(lib, pkg)
{
    where <- match(paste("package:", pkg, sep = ""), search())
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
    packageStartupMessage(paste(title, " (version ", ver, ")", sep = ""))
    #library.dynam("rAverage", pkg, lib)
}

.onLoad <- function(lib, pkg)
    library.dynam("rAverage", pkg, lib)
