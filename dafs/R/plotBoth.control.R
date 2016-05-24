.genPlots = TRUE
.embedF = TRUE
.embedFoptions = "-sPAPERSIZE=a4 -dPDFSETTINGS=/printer -dCompatibilityLevel=1.3 -dMaxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -dEmbedAllFonts=true -dUseCIEColor=true -sProcessColorModel=DeviceGray"

plotBoth.control = function(genPlots = .genPlots, embedF = .embedF, embedFoptions = .embedFoptions){
    .genPlots = genPlots
    .embedF = embedF
    .embeFoptions = embedFoptions
    return(list(genPlots = genPlots, embedF = embedF, embedFoptions = embedFoptions))
}
