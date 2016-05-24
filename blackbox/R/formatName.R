formatName <-
function (charname, format = "expression", Nbstyle = T)
{
    if (class(charname) == "expression")
        stop.redef("(!) from formatName: invalid call with argument of class 'expression'.")
    tmp <- charname
    if (charname == "twoNmu")
        tmp <- AUEformat("2Nmu", "", "2*italic(N)*mu", format = format)
    if (charname == "twoNmuSEQ")
        tmp <- AUEformat("2NmuSEQ", "", "2*italic(N)*mu[SEQ]",
            format = format)
    if (charname == "twoNm")
        tmp <- AUEformat("2Nm", "2Nm", "2*italic(N)*italic(m)",
            format = format)
    if (charname == "g")
        tmp <- AUEformat("g", "g", "italic(g)", format = format)
    if (charname == "D")
        tmp <- AUEformat("Dg/2N", "D[g]/2N", "italic(D)[g]/2*italic(N)",
            format = format)
    if (charname == "T")
        tmp <- AUEformat("Tg/2N", "T[g]/2N", "italic(T)[g]/2*italic(N)",
            format = format)
    if (charname == "pGSM")
        tmp <- AUEformat("pGSM", "p[GSM]", "italic(p)[GSM]",
            format = format)
    if (charname == "twoNancmu")
        tmp <- AUEformat("2Nancmu", "", "2*italic(N)[anc]*mu",
            format = format)
    if (charname == "twoNfoundermu")
        tmp <- AUEformat("2Nfoundermu", "", "2*italic(N)[founder]*mu",
            format = format)
    if (charname == "Nratio")
        tmp <- AUEformat("Nratio", "Nratio", "italic(N)[ratio]",
            format = format)
    if (charname == "NactNfounderratio")
        tmp <- AUEformat("NactNfounder-ratio", "NactNfounder-ratio",
            "italic(N)[act]/italic(N)[founder]-ratio", format = format)
    if (charname == "NfounderNancratio")
      tmp <- AUEformat("NfounderNanc-ratio", "NfounderNanc-ratio",
                       "italic(N)[founder]/italic(N)[anc]-ratio", format = format)
    if (charname == "Nratio")
        tmp <- AUEformat("Nratio", "Nratio", "italic(N)-ratio",
            format = format)
    if (charname == "M1" || charname == "twoNm1")
        tmp <- AUEformat("M1", "", "2*italic(N)[1]*italic(m)[12]",
            format = format)
    if (charname == "M2" || charname == "twoNm2")
        tmp <- AUEformat("M2", "", "2*italic(N)[2]*italic(m)[21]",
            format = format)
    if (charname == "Q1")
        tmp <- AUEformat("Q1", "", "italic(Q)[1]", format = format)
    if (charname == "latt2Ns2")
        tmp <- AUEformat(if (Nbstyle) {
            "Nb"
        }
        else {
            "2Ds2"
        }, if (Nbstyle) {
            "Nb"
        }
        else {
            ""
        }, if (Nbstyle) {
            "Nb"
        }
        else {
            "2*italic(N)*sigma^2"
        }, format = format)
    if (charname == "condS2")
        tmp <- AUEformat("s2cond", "", "sigma~scriptstyle(phantom()[scriptstyle(cond)]^{scriptstyle(2)})",
            format = format)
    return(tmp)
}
