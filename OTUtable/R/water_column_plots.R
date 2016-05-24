# Functions for analyzing the metadata

make_do_matrix <- function(sampleID, field_data){
  # Pull out only entries for sample subset
  find <- grep(sampleID, field_data$Sample_Name)
  field_data <- field_data[find, ]
  # Get date and depth information
  dates <- unique(extract_date(field_data$Sample_Name))
  depth <- sort(unique(field_data$Depth))
  # Set up matrix
  data_matrix <- matrix(0, nrow = length(depth), ncol = length(dates))
  
  # Put depths in rows and dates in columns, put data point in correct coordinates
  for(i in 1:length(depth)){
    row <- field_data[which(field_data$Depth == depth[i]), ]
    for(j in 1:length(dates)){
      col <- row[which(extract_date(row$Sample_Name) == dates[j]), ]
      data_matrix[i, j] <- col$DO[1]
    }
  }
  # Add labels to matrix and order by date
  rownames(data_matrix) <- depth
  colnames(data_matrix) <- as.character(dates)
  data_matrix <- data_matrix[, order(dates)]
  
  # Fill in missing datapoints by averaging depth above and below
  for(i in 1:(dim(data_matrix)[1] - 1)){
    for(j in 1:dim(data_matrix)[2]){
      if(is.na(data_matrix[i, j]) == T){
        data_matrix[i, j] <- mean(c(data_matrix[i - 1, j], data_matrix[i + 1, j]))
      }
    }
  }
  #If missing point is the deepest point, take the second deepest point
  i = dim(data_matrix)[1]
  for(j in 1:dim(data_matrix)[2]){
    if(is.na(data_matrix[i, j]) == T){
      data_matrix[i, j] <- data_matrix[i - 1, j]
    }
  } 
  return(data_matrix)
}

make_temp_matrix <- function(sampleID, field_data){
  # Pull out only entries for sample subset
  find <- grep(sampleID, field_data$Sample_Name)
  field_data <- field_data[find, ]
  # Get date and depth information
  dates <- unique(extract_date(field_data$Sample_Name))
  depth <- sort(unique(field_data$Depth))
  # Set up matrix
  data_matrix <- matrix(0, nrow = length(depth), ncol = length(dates))
  
  # Put depths in rows and dates in columns, put data point in correct coordinates
  for(i in 1:length(depth)){
    row <- field_data[which(field_data$Depth == depth[i]), ]
    for(j in 1:length(dates)){
      col <- row[which(extract_date(row$Sample_Name) == dates[j]), ]
      data_matrix[i, j] <- col$Temperature[1]
    }
  }
  # Add labels to matrix and order by date
  rownames(data_matrix) <- depth
  colnames(data_matrix) <- as.character(dates)
  data_matrix <- data_matrix[, order(dates)]
  
  # Fill in missing datapoints by averaging depth above and below
  for(i in 1:(dim(data_matrix)[1] - 1)){
    for(j in 1:dim(data_matrix)[2]){
      if(is.na(data_matrix[i, j]) == T){
        data_matrix[i, j] <- mean(c(data_matrix[i - 1, j], data_matrix[i + 1, j]))
      }
    }
  }
  # If missing point is the deepest point, take the second deepest point
  i = dim(data_matrix)[1]
  for(j in 1:dim(data_matrix)[2]){
    if(is.na(data_matrix[i, j]) == T){
      data_matrix[i, j] <- data_matrix[i - 1, j]
    }
  }
  
  return(data_matrix)
}

# Plot contour water columns

rotate <- function(data_matrix) t(apply(data_matrix, 2, rev))


plot_column <- function(data_matrix, title){
    filled.contour(rotate(data_matrix), ylab = "Depth (m)", main = paste(title),
                 color.palette = colorRampPalette(c("blue", "cyan", "green", "yellow", "red")),
                 plot.axes = {
                   axis(2, at = seq(from=0, to = 1, length.out = length(rownames(data_matrix))), labels = rev(rownames(data_matrix)))
                   axis(1, at = seq(from=0, to = 1, length.out = length(colnames(data_matrix))), labels =colnames(data_matrix), cex.axis=0.9, las=2)
                 }
  )
}
