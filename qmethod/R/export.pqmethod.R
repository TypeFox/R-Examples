export.pqmethod <- function(dataset, study.name, study.description, col.range, filename='Q_data_forPQmethod', left.zeros, right.zeros, statements) {
  # Warning: no checks over the data are made
  message('Note that this function does not check the data: whether there are duplicate Q-sorts, whether there are Q-sorts not conforming to the distribution, etc.')
  # Description max. 68 characters
  study.description <- substr(study.description, 1, 68)
  # Basic info
  n.qsorts <- ncol(dataset)
  n.stat <- nrow(dataset)
  filename.full <- paste0(filename, '.dat')
  filename.full.sta <- paste0(filename, '.sta')
  # Create empty object to place the data
  pqm.data <- matrix(as.character(NA), nrow=n.qsorts+2)
  # First line: number of Q sorts, of statements, and name and description of study
  pqm.data[1] <- paste0('  0',
                        sprintf("% 3d",n.qsorts),
                        sprintf("% 3d",n.stat),
                        ' ', study.name, ': ',
                        study.description)
  # Second line: still wondering what it is about
  distro <- as.vector(table(sort(dataset[,1]))) # grid array
  pqm.data[2] <- paste0(c(sprintf("% 3d",col.range),
                          rep("  0", left.zeros),
                          sprintf("% 3d",distro),
                          rep("  0", right.zeros)),
                        collapse='')
  
  
  # Third line onwards: Q-sorts, with id in the first ten characters, and values after (two characters for each statement)
  for (i in 1:ncol(dataset)) {
    q.sort.name <- sprintf("%- 10s",colnames(dataset)[i])
    q.sort.values <- paste0(sprintf("%- 2d",dataset[,i]), collapse='')
    pqm.data[2+i] <- paste0(q.sort.name, q.sort.values)
  }
  # Write file of Q-sorts
  fileConn <- file(filename.full)
  writeLines(pqm.data, fileConn)
  close(fileConn)
  # Write file of statements
  fileConn <- file(filename.full.sta)
  writeLines(substr(statements, 1, 50), fileConn)
  close(fileConn)
  # Final info
  message(paste0('----------------------------------------------\nData exported to PQMethod format (*.dat file).\nName of study: ',study.name,'\nNumber of Q-sorts: ', n.qsorts, '\nNumber of statements: ', n.stat))
}