fraction.formula <-
function (vecteur, indices) 
{
    vecteur <- gsub(" ", "", vecteur)
    Gauche = vector(length = 0)
    Droite = list()
    if (grepl("~~", vecteur[indices[1]])) {
        signe = "~~"
    }
    else {
        signe = "=~"
    }
    calc.indice <- function(k) {
        equ <- strsplit(vecteur[k], split = signe)[[1]]
        equ <- gsub(" ", "", equ)
        Vgauche <- strsplit(equ[1], split = " ")[[1]]
        Vgauche <- gsub(" ", "", Vgauche)
        Vgauche[Vgauche != ""]
        if (length(Vgauche) != 1) {
            stop("mauvaise de specification du modele")
        }
        Vdroite <- strsplit(equ[2], split = "[+ ]")[[1]]
        Vdroite <- gsub(" ", "", Vdroite)
        Vdroite[Vdroite != ""]
        if (length(Vdroite) < 1) {
            stop("mauvaise de specification du modele")
        }
        Gauche <<- c(Gauche, Vgauche)
        Droite <<- c(Droite, list(Vdroite))
    }
    calc.indice <- Vectorize(calc.indice)
    calc.indice(indices)
    return(list(Gauche, Droite))
}
