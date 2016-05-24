# 0. Run a program to generate outputs
setwd("C:/aErer"); source("r072sunSJAF.r", echo = FALSE)

# A1. Export two text tables to two seperate files
write.table(x = table.3, file = "OutInsTable3.csv", sep = ",")
write.table(x = table.4, file = "OutInsTable4.csv", sep = ",")

# A2. Export two text tables to one file
# A warning message will pop up "append = TRUE"; fine in this case.
write.table(x = table.3, file = "OutInsTableAll.csv", sep = ",",
  append = FALSE)
write.table(x = table.4, file = "OutInsTableAll.csv", sep = ",",
  append = TRUE)

# A3. Use listn() and write.list() to export multiple tables to one file
# listn() assigns list names automatically.
out.a <- list(table.3 = table.3, table.4 = table.4)
out.b <- listn(table.3, table.4)
identical(out.a, out.b)  # TRUE

# OutInsTableAll2 is similar to OutInsTableAll (small format differences)
out <- listn(table.3, table.4)
write.list(z = out, file = "OutInsTableAll2.csv")

# B. Export multiple tables to Excel
library(xlsx)  # load the library
write.xlsx(x = table.3, file = "OutInsTableAll.xlsx",
  sheetName = "table.3", row.names = FALSE, append = FALSE)
write.xlsx(x = table.4, file = "OutInsTableAll.xlsx",
  sheetName = "table.4", row.names = FALSE, append = TRUE)

# Many tables to Excel with a loop
name <- paste("table", c(3, 4), sep = ".")
for (i in 1:length(name)) {
  write.xlsx(x = get(name[i]), file = "OutInsTableAll.xlsx",
    sheetName = name[i], row.names = FALSE, append = as.logical(i - 1))
}

# C. Export graphs
png(file = "testFigure1.png")
  plot(1:100)
dev.off()