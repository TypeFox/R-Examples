library('diseasemapping')
data('kentucky')

head(larynx)
10^5*larynxRates[paste(c("M","F"), 50, sep="_")]

kentucky1 = getSMR(kentucky@data, larynxRates)
kentucky1[1:4,c(1,2,grep("expected", names(kentucky1),ignore.case=TRUE))]

kentucky1 = getSMR(kentucky, larynxRates)
kentucky1@data[1:4,c(1,2,grep("expected", names(kentucky1),ignore.case=TRUE))]

if(require('mapmisc', quietly=TRUE)) {
#  kmap = openmap(kentucky)
  col = colourScale(
      kentucky1$expected,
      style='fixed',  
      breaks=c(0:5,max(kentucky1$expected)), 
    dec=0,opacity=c(0.6,1)
  )

  map.new(kentucky)
#  plot(kmap,add=TRUE)
  plot(kentucky1, col=col$plot,add=TRUE)
  legendBreaks('topleft', col)
}

junk = getSMR(kentucky, larynxRates, regionCode='junk')


kentucky2 = getSMR(kentucky@data, larynxRates, 
    larynx, 
    regionCode="County")
kentucky2[1:4,c(1,2,grep("expected|observed", names(kentucky2),ignore.case=TRUE))]


kentucky2 = getSMR(kentucky, 
    larynxRates, 
    casedata=larynx, 
    regionCode="County")
kentucky2@data[1:4,c(1,2,grep("expected|observed", names(kentucky2),ignore.case=TRUE))]

if(require('mapmisc', quietly=TRUE)) {


  col = colourScale(
        kentucky2$observed,
        col='RdYlBu',
        style='quantile', 
        breaks=12, dec=0,opacity=c(0.6,1),
        rev=TRUE
    )
    
    map.new(kentucky)
 #   plot(kmap,add=TRUE)
    plot(kentucky1, col=col$plot,add=TRUE)
    legendBreaks('topleft', col)

    
  col = colourScale(
      kentucky2$expected,
  style='fixed', 
  col=col$col,
  breaks=col$breaks,opacity=c(0.6,1)
  )
  
  map.new(kentucky)
#  plot(kmap,add=TRUE)
  plot(kentucky2, col=col$plot,add=TRUE)
  legendBreaks('topleft', col)
  
}


kentucky3 = getSMR(kentucky@data, 
    model=list(larynxRates, larynxRates*2)
)
kentucky3[1:4,c(1,2,grep("expected|observed", names(kentucky3),ignore.case=TRUE))]

kentucky3 = getSMR(kentucky, 
    model=list('1990'=larynxRates, '1991'=larynxRates*2)
)
kentucky3@data[1:4,c(1,2,grep("expected|observed", names(kentucky3),ignore.case=TRUE))]

modelList = list()
for (D in 3:12) {
  modelList[[
      as.character(D)
      ]] = larynxRates*D/5
}


kentucky4 = getSMR(list(
        '5'=kentucky, '10'=kentucky
    ),
    model=modelList
)
kentucky4[[1]]@data[1:4,c(1,2,grep("expected|observed", 
            names(kentucky4[[1]]),ignore.case=TRUE))
]
kentucky4[[2]]@data[1:4,c(1,2,grep("expected|observed", 
            names(kentucky4[[2]]),ignore.case=TRUE))
]
 
if(require('mapmisc', quietly=TRUE)) {
  
  
  col = colourScale(
      kentucky4[[1]]$expected_6,
      col='RdYlBu',
      style='quantile', 
      breaks=15, dec=0,opacity=c(0.6,1),
      rev=TRUE
  )
  
  map.new(kentucky)
#  plot(kmap,add=TRUE)
  plot(kentucky4[[1]], col=col$plot,add=TRUE)
  legendBreaks('topleft', col)
  
  
  col = colourScale(
      kentucky4[[2]]$expected_11,
      col='RdYlBu',
      style='quantile', 
      breaks=15, dec=0,opacity=c(0.6,1),
      rev=TRUE
  )
  
  map.new(kentucky)
#  plot(kmap,add=TRUE)
  plot(kentucky4[[2]], col=col$plot,add=TRUE)
  legendBreaks('topleft', col)
  
}
