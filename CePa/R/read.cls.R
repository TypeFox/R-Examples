read.cls <-
function (file, treatment, control) {
    label = readLines(file, n = -1)
    label = label[3]
    label = unlist(strsplit(label, " |\t"));
    return(sampleLabel(label = label, treatment = treatment, control = control))
}
