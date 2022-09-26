
# Write a function to calculate purity
calculate_purity_values <- function(one_row_df, main_df){
  col_names <- colnames(one_row_df)
  # cat("Column names of small df:", col_names, "\n")
  for(i in 1:length(col_names)){
    current_col <- col_names[i]
    # cat("Current column:", current_col, "\n")
    # cat("Summing column:", current_col, "Sum of current main df col sum:", sum(main_df[, current_col]), "\n")
    one_row_df[1, i] <- round((one_row_df[1, i] / sum(main_df[, current_col])) * 100, 2)
  }
  return(one_row_df)
}



process_each_sequence <- function(uniq_asv_df, top_n){
  if(top_n > ncol(uniq_asv_df)){
    print(paste("Warning: Given top_n value is more than the total num of plate-well combinations. Therefore, generating output for maximum possible combinations,", ncol(uniq_asv_df)))
    top_n <- ncol(uniq_asv_df)
  }
  final_df <- as.data.frame(matrix(nrow=0, ncol=top_n+1))
  colnames(final_df) <- c("ASV", c(paste0("top_count_", 1:top_n)))
  for(i in 1:nrow(uniq_asv_df)){
    # print(i)
    options(warn=-1)
    sorted_row <- sort(uniq_asv_df[i, ], decreasing = TRUE)
    one_row_df <- as.data.frame(sorted_row[,1:top_n])
    seq_purity <- calculate_purity_values(one_row_df, uniq_asv_df)
    final_row <- c(rownames(one_row_df), paste(colnames(one_row_df), one_row_df[1,], seq_purity[1,], sep = " | "))
    final_df[nrow(final_df) + 1,] <- final_row
    # print(one_row_df)
    # print(seq_purity)
    options(warn=0)
  }
  return(final_df)
}


get_taxonomy <- function(uniq_asv){
  taxonomy_added <- assignTaxonomy(uniq_asv, "~/Documents/silva_db/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
  return(taxonomy_added)
}

