organize_replicates <- function(dataframe, min_range, max_range) {
  # Removing intensity columns
  intensity_columns = c(2)
  pt = 0
  for (i in 3:length(dataframe)) {
    if (pt == 2) {
      intensity_columns = append(intensity_columns, i)
      pt = 0
    } else {
      pt = pt + 1
    }
  }; remove(i, pt)
  dataframe = dataframe[,-intensity_columns]
  
  # Retrieving and correcting column names
  colnames_list = strsplit(colnames(dataframe),'...', fixed = TRUE)
  colnames_cath = c()
  for (i in 1:length(colnames_list)) {
    colnames_cath = append(colnames_cath, colnames_list[[i]][1])
  }; remove(i, colnames_list)
  
  # Storing Labels
  label = c()
  for (i in 1:length(colnames_cath)) {
    if (i %% 2 == 0) {
      label = append(label, colnames_cath[i])
    }
  }; remove(i)
  
  # Defining m/z range
  range_start = min_range
  range_end = max_range
  
  # Creating dataframe to contain results
  df_result = matrix(ncol = ((range_end-range_start)+2), nrow = 0)
  df_result = data.frame(df_result)
  
  # Loop that organizes the data
  numb_list = c()
  pt_label = 0
  for (i in 1:ncol(dataframe)) {
    # Selecting columns
    if (i %% 2 != 0) {
      # Capturing mz values, removing NA and converting to integers
      mz_value = unname(unlist(dataframe[,i]))
      mz_value = mz_value[!is.na(mz_value)]
      mz_value = as.integer(mz_value)
      
      # Capturing relative values and removing NA
      relative_value = unname(unlist(dataframe[,i+1]))
      relative_value = relative_value[!is.na(relative_value)]
      
      # Loop that sorted the m/z and frequencies of the sample according to the range
      lista = rep(0.00,range_end)
      for (i in 1:length(mz_value)) {
        lista[mz_value[i]] = relative_value[i]
      }
      
      # Adding list to result table
      df_result = rbind(df_result, lista[range_start:range_end])
    }
  }
  
  # Editing column names
  colnames(df_result) = c(range_start:range_end)
  
  # Adding column with labels
  df_result = cbind(label, df_result)
  
  # Sorting the rows by the label column
  df_result = df_result[order(df_result$label),]
  
  # Check number of replicates
  labels_unique = unique(df_result[,1])
  label_list = c()
  n_label = c()
  warning_label = c()
  for (e in 1:length(labels_unique)) {
    pt_count = 0
    for (i in 1:nrow(df_result)) {
      if (labels_unique[e] == df_result[i,1]) {
        pt_count = pt_count + 1
      }
    }
    if (pt_count < 3) {
      warning_label = append(labels_unique[e], warning_label)
    }
    label_list = append(labels_unique[e], label_list)
    n_label = append(pt_count, n_label)
  }
  
  # Showing number of replicas per lebal
  df_n_label = data.frame(label_list,n_label)
  colnames(df_n_label) <- c('Label', 'n_replicates')
  print(df_n_label)
  
  # Showing lebal with less than 3 replicates
  if (warning_label != 0) {
    writeLines("\nWarning!!\nThe labels below contain less than 3 replicates:\n")
    print(warning_label)
  }
  return(df_result)
}


means_replicates <- function(dataframe) {
  # identifying unique labels
  label = unique(dataframe[,1])
  
  # creating dataframe to contain results
  df_result = matrix(ncol = ncol(dataframe), nrow = 0)
  df_result = data.frame(df_result)
  
  # calculating the average of replicates
  for (i in 1:length(label)) {
    label_i =  label[i]
    df_temp = subset(dataframe, dataframe[,1] == label_i)
    media_temp = colMeans(df_temp[2:ncol(df_temp)])
    media_temp = round(media_temp, 3)
    df_result = rbind(df_result, media_temp)
  }; remove(i,label_i,df_temp,media_temp)
  
  # renaming result table (with averages) rows and columns
  colnames(df_result) = colnames(dataframe[2:length(dataframe)])
  
  # adding column with labels
  df_result = cbind(label, df_result)
  
  # sorting the rows by the label column
  df_result = df_result[order(df_result$label),]
  
  # removing unnecessary variables
  remove(label)
  
  return(df_result)
}