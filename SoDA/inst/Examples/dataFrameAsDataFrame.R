setAs("data.frame", "dataFrame1",
      function(from) 
        new("dataFrame1",
             as.list(from),
             row.names = rownames(from),
             names = colnames(from))
      )
