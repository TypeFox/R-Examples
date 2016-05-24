`diversityvariables` <-
function(x, y, digits=8){
    y$richness <- diversityresult(x, index='richness', method='s', digits=digits)[,1]
    y$Shannon <- diversityresult(x, index='Shannon', method='s', digits=digits)[,1]
    y$Simpson <- diversityresult(x, index='Simpson', method='s', digits=digits)[,1]
    y$inverseSimpson <- diversityresult(x, index='inverseSimpson', method='s', digits=digits)[,1]
    y$Logalpha <- diversityresult(x, index='Logalpha', method='s', digits=digits)[,1]
    y$Berger <- diversityresult(x, index='Berger', method='s', digits=digits)[,1]
    y$Jevenness <- diversityresult(x, index='Jevenness', method='s', digits=digits)[,1]
    y$Eevenness <- diversityresult(x, index='Eevenness', method='s', digits=digits)[,1]
    y$richness <- diversityresult(x, index='richness', method='s', digits=digits)[,1]
    return(y)
}

