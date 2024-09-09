
##6-Part ITS## JJM ##


#Example Input
# df <- dat_1
# ILRModel1 <- cph1_vpa 
# adjust_col1 <- 3
# adjust_col2 <- 2
# column_names
# head(compo)
# raw_column_names <- c("mpa", "lpa", "vpa",  "slp"  , "stand" , "sit")

ITS_6_part <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
  # Calculate the mean composition and convert to minutes
  compomean <- as.vector(mean(acomp(compo)))
  compomeanmins <- as.vector(compomean * 60 * 24)
  
  # Initialize data frames for the adjusted compositions
  CompoIT_plus <- as.data.frame(t(compomeanmins))
  CompoIT_minus <- as.data.frame(t(compomeanmins))
  duprow_plus <- CompoIT_plus[1, ]
  duprow_minus <- CompoIT_minus[1, ]
  
  # Create 301 rows for adjustments from -300 to +300
  for (i in 1:300) {
    CompoIT_plus <- rbind(CompoIT_plus, duprow_plus)
    CompoIT_minus <- rbind(CompoIT_minus, duprow_minus)
  }
  
  # Adjust the selected behaviors in both directions
  for (i in 1:300) {
    CompoIT_plus[i + 1, adjust_col1] <- CompoIT_plus[i + 1, adjust_col1] + i  
    CompoIT_plus[i + 1, adjust_col2] <- CompoIT_plus[i + 1, adjust_col2] - i  
    CompoIT_minus[i + 1, adjust_col1] <- CompoIT_minus[i + 1, adjust_col1] - i  
    CompoIT_minus[i + 1, adjust_col2] <- CompoIT_minus[i + 1, adjust_col2] + i  
  }
  
  # Combine datasets and filter out invalid rows
  CompoIT_combined <- rbind(CompoIT_minus, CompoIT_plus)
  CompoIT_combined_raw <- CompoIT_combined[!apply(CompoIT_combined < 0 | is.na(CompoIT_combined), 1, any), ]
  
  # Debugging: Print structure of CompoIT_combined_raw
  print(paste("Number of columns in CompoIT_combined_raw:", ncol(CompoIT_combined_raw)))
  
  # Normalize back to proportions and convert to data frame
  CompoIT_combined <- as.matrix(CompoIT_combined_raw / (60 * 24))
  CompoIT_combined <- as.data.frame(CompoIT_combined)
  
  
  # ILR transformation
  New_Z <- LR_ILR_6part(CompoIT_combined, base_rot = 1)
  CompoIT_combined_raw
  colnames(New_Z) <- column_names
  colnames(CompoIT_combined_raw) <- raw_column_names
  CompoIT_combined_raw
  
 
 # Define a list to store the predictions
 predictions_list <- list()
 
 # Loop through each row in New_Z
 for (i in 1:nrow(New_Z)) {
   # Create a list to store arguments for rms::Predict dynamically
   predict_args <- list(ILRModel1, fun = exp, ref.zero = TRUE)
   
   # Dynamically add arguments based on column names
   for (col in column_names) {
     predict_args[[col]] <- New_Z[[col]][i]
   }
   
   # Call rms::Predict with the dynamically created argument list
   P2 <- do.call(rms::Predict, predict_args)
   
   # Store the prediction for this iteration
   predictions_list[[i]] <- P2
 }
 
 
results <- do.call(rbind, predictions_list)
results <- cbind(results, CompoIT_combined_raw)
  
return(results)
}



# ## Raw function 
# predictions_list <- list()
# 
# # Loop through each row in New_Z
# for (i in 1:nrow(New_Z)) {
#   # Use rms::Predict for the ith row
#   P2 <- rms::Predict(
#     ILRModel1, 
#     VPAvsAll = New_Z$VPAvsAll[i], 
#     SLPvsRemaining = New_Z$SLPvsRemaining[i],
#     StandvsRemaining = New_Z$StandvsRemaining[i],
#     SitvsRemaining = New_Z$SitvsRemaining[i],
#     MPAvsLPA = New_Z$MPAvsLPA[i],
#     fun = exp,
#     ref.zero = TRUE
#   )
#   
#   # Store the prediction for this iteration
#   predictions_list[[i]] <- P2
# }