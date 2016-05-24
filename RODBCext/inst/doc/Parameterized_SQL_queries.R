## ----eval = FALSE--------------------------------------------------------
#  library(RODBC)
#  
#  connHandle <- odbcConnect("cakesDatabase")
#  newData <- read.csv("newData.csv", stringsAsFactors = F)
#  
#  for(row in 1:nrow(newData)){
#    query <- paste0(
#      "UPDATE cakes
#       SET price = ", newData$price[row], "
#       WHERE cake = '", newData$cake[row], "'"
#    )
#    sqlQuery(connHandle, query)
#  }
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBC)
#  
#  connHandle <- odbcConnect('studentsDatabase')
#  newStudents <- read.csv('newStudents.csv', stringsAsFactors = F)
#  
#  for(row in 1:nrow(newStudents)){
#    query <- paste0(
#      "INSERT INTO students (first_name, last_name)
#       VALUES (
#         '", newStudents$first_name[row],"',
#         '", newStudents$last_name[row],"',
#       )"
#    )
#    sqlQuery(P, query)
#  }
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  
#  connHandle <- odbcConnect("cakesDatabase")
#  newData <- read.csv("newData.csv", stringsAsFactors = F)
#  
#  query <- "UPDATE cakes SET price = ? WHERE cake = ?"
#  for(row in 1:nrow(newData)){
#    sqlExecute(connHandle, query, newData[i, ])
#  }
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  
#  connHandle <- odbcConnect("cakesDatabase")
#  newData <- read.csv("newData.csv", stringsAsFactors = F)
#  
#  query <- "UPDATE cakes SET price = ? WHERE cake = ?"
#  sqlExecute(connHandle, query, newData)
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect('EWD') # my sample ODBC database
#  data <- data.frame(1:10000, letters[rep(1:10, 1000)])
#  
#  # Ordinary query - paste0() called in every loop
#  system.time({
#    for(row in 1:nrow(data)){
#      query <- paste0("INSERT INTO my_table VALUES (", data[row, 1], "'", data[row, 2],"')")
#      sqlQuery(connHandle, query)
#    }
#  })
#  #   user  system elapsed
#  #  5.384   2.288  16.397
#  
#  # Ordinary query - paste0() called only once
#  system.time({
#    queries <- paste0(
#      "INSERT INTO my_table VALUES (", data[, 1], "'", data[, 2],"')"
#    )
#    for(query in queries){
#      sqlQuery(connHandle, query)
#    }
#  })
#  #   user  system elapsed
#  #  2.088   2.028   7.255
#  
#  # Parameterized query
#  system.time({
#    sqlExecute(connHandle, "INSERT INTO my_table VALUES (?, ?)", data)
#  })
#  #   user  system elapsed
#  #  0.300   0.232   3.935
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect('EWD') # my sample ODBC database
#  
#  pupils = sqlQuery(
#    connHandle, "SELECT id_obserwacji FROM obserwacje LIMIT 10000",
#    stringsAsFactors = F
#  )[, 1]
#  
#  # Ordinary query - paste0() called in every loop
#  system.time({
#    for(i in pupils){
#      query <- paste0(
#        "SELECT count(*)
#         FROM testy_obserwacje JOIN testy USING (id_testu) JOIN arkusze USING (arkusz)
#         WHERE id_obserwacji = ", pupils[i]
#      )
#      tmp <- sqlQuery(connHandle, query)
#      # some other computations here
#    }
#  })
#  #   user  system elapsed
#  # 10.896   1.508  61.424
#  
#  # Ordinary query - paste0() called only once
#  system.time({
#    queries <- paste0(
#      "SELECT count(*)
#       FROM testy_obserwacje JOIN testy USING (id_testu) JOIN arkusze USING (arkusz)
#       WHERE id_obserwacji = ", pupils
#    )
#    for(query in queries){
#      tmp <- sqlQuery(connHandle, query)
#      # some other computations here
#    }
#  })
#  #   user  system elapsed
#  # 11.016   1.108  51.766
#  
#  # Parameterized query
#  system.time({
#    query = "
#      SELECT count(*)
#      FROM testy_obserwacje JOIN testy USING (id_testu) JOIN arkusze USING (arkusz)
#      WHERE id_obserwacji = ?"
#    sqlPrepare(connHandle, query)
#    for(i in pupils){
#      tmp = sqlExecute(connHandle, NULL, pupils[i], fetch=T)
#      # some other computations here
#    }
#  })
#  #   user  system elapsed
#  # 12.140   0.312  26.468

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect("myDatabase")
#  
#  # good old RODBC call
#  data <- sqlQuery(connHandle, "SELECT * FROM myTable WHERE column = 'myValue'")
#  # RODBCext equivalent
#  data <- sqlExecute(connHandle, "SELECT * FROM myTable WHERE column = ?", 'myValue', fetch = TRUE)
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect("myDatabase")
#  
#  filterData <- data.frame('column1' = 1:5, column2 = c('a', 'b', 'c', 'd', 'e'))
#  data <- sqlExecute(connHandle, "SELECT * FROM myTable WHERE column1 = ? AND column2 = ?", filterData, fetch = TRUE)
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect("myDatabase")
#  
#  sqlExecute(connHandle, "SELECT * FROM myTable WHERE column = ?", 'myValue', fetch = FALSE)
#  data <- sqlGetResults(connHandle, max = 10) # fetch no more than 10 first rows
#  # data processing comes here
#  data <- sqlGetResults(connHandle) # fetch all other rows
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect("myDatabase")
#  
#  sqlExecute(
#    connHandle, "SELECT * FROM myTable WHERE column = ?", 'myValue',
#    fetch = TRUE, stringsAsFactors = FALSE, dec = ",", max = 50, as.is = TRUE
#  )
#  
#  odbcClose(connHandle)

## ----eval = FALSE--------------------------------------------------------
#  library(RODBCext)
#  connHandle <- odbcConnect('myDatabase')
#  
#  sqlPrepare(connHandle, "SELECT * FROM myTable WHERE column = ?") # prepare query
#  
#  # for some reason (e.g. resources limits) data must be processed sequentialy
#  foreach(i in observations){
#    data = sqlExecute(connHandle, NULL, i$column, fetch=T)
#    # data processing for a given observations goes here
#  }
#  odbcClose(connHandle)

