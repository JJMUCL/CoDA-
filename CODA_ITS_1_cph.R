##Isotemporal Substitution 

##ITS1 - original function - 4 part composition with linear regression

##ITS - 4-part composition - 4 part composition with GLM

##ITS_6_part - 6 part composition fit with Cox Prop Hazards


ITS1 <- function(df, ILRModel1) { 
  ###ILR  - No ADJ1 SAMPLE - UNADJUSTED
  compo <- df[,2:5]
  # CompCog <-df$Composite_z
  # data <- cbind(compo, CompCog)
  # # First rotation model [1] MVPA [2] SB [3] LIPA [4] SLEEP 
  # # compute ilr 
  # Zi<-LR_ILR(compo,1) 
  # #make the model
  # model_data<-cbind(Zi,CompCog)
  # fiti<-lm(CompCog~.,data=model_data)
  # summary(fiti)
  # confint(fiti,'mvpa', level=0.95)
  # confint(fiti,'sedent', level=0.95)
  # confint(fiti,'LIPA', level=0.95)  
  # summary(fiti)
  
  compomean <- as.vector(mean(acomp(compo)))
  compomeanmins <- as.vector(compomean*60*24)
  CompoIT<- as.data.frame(t(compomeanmins))
  CompoIT<- dplyr::rename(CompoIT, mvpa=V1,
                          sedent = V2,
                          LIPA = V3,
                          SLEEP = V4)
  duprow = CompoIT[1,]
  for(i in 1:300)
  {
    CompoIT = rbind(CompoIT,duprow)
  }
  for(i in 1:300) {
    CompoIT$mvpa[i] <- CompoIT$mvpa[i]+i
    CompoIT$sedent[i] <- CompoIT$sedent[i]-i}
  CompoIT <- as.matrix(CompoIT/60/24)
  CompoIT <- as.data.frame(CompoIT)
  
  
  newdata<-as.data.frame(CompoIT) #put them together into one dataframe
  names(newdata)<-c("mvpa","sedent","LIPA", "SLEEP")
  
 
  
  New_Z<-LR_ILR(newdata)
  # names(New_Z)<-c("", "LIPA","sedent","SLEEP")
  
  
  #use model to make prediction
  P1<-predict(ILRModel1, New_Z, interval="confidence") 
  P1<- as.data.frame(P1)
  print(P1)
  return(as.data.frame(P1))
}

# 
# 
# ITS <- function(df, ILRModel1, adjust_col1, adjust_col2) {
#   # Use column indices for the four behaviors
#   compo <- df[, 2:5]
#   
#   # Calculate the mean composition and convert it to minutes
#   compomean <- as.vector(mean(acomp(compo)))
#   compomeanmins <- as.vector(compomean * 60 * 24)
#  
#   # Convert to a data frame and create multiple rows for different adjustments
#   CompoIT <- as.data.frame(t(compomeanmins))
#   duprow <- CompoIT[1,]
#   
#   for (i in 0:600) {
#     CompoIT <- rbind(CompoIT, duprow)
#   }
#   
#   # Adjust the selected behaviors
#   for (i in 0:600) {
#     CompoIT[i, adjust_col1] <- CompoIT[i, adjust_col1] + i  # Add minutes to the first behavior
#     CompoIT[i, adjust_col2] <- CompoIT[i, adjust_col2] - i  # Subtract minutes from the second behavior
#   }
#   
#   # Normalize back to proportions
#   CompoIT <- as.matrix(CompoIT / (60 * 24))
#   CompoIT <- as.data.frame(CompoIT)
#   
#   # ILR transformation
#   New_Z <- LR_ILR(CompoIT, base_rot = 1)  # Assuming base rotation 1 for now
#   colnames(New_Z) <- colnames(ILRModel1$model)[-1]
#   # Predict using the model
#   P2 <- predict(ILRModel1, New_Z, interval = "confidence")
#   P2 <- as.data.frame(P2)
#   
#   
#   colnames(New_Z) <- colnames(fiti$model)[-1]
#   # Predict using the model
#   P2 <- predict(fiti, New_Z, interval = "confidence")
#   P2 <- as.data.frame(P2)
#   
#   
#   return(P2)
# }


ITS <- function(df, ILRModel1, adjust_col1, adjust_col2) {
  # Use column indices for the four behaviors
  # Calculate the mean composition and convert it to minutes
  compomean <- as.vector(mean(acomp(compo)))
  compomeanmins <- as.vector(compomean * 60 * 24)
  
  # Convert to a data frame and create multiple rows for different adjustments (both + and -)
  CompoIT_plus <- as.data.frame(t(compomeanmins))
  CompoIT_minus <- as.data.frame(t(compomeanmins))
  duprow_plus <- CompoIT_plus[1,]
  duprow_minus <- CompoIT_minus[1,]
  
  # Create 301 rows for adjustments from -300 to +300
  for (i in 1:300) {
    CompoIT_plus <- rbind(CompoIT_plus, duprow_plus)
    CompoIT_minus <- rbind(CompoIT_minus, duprow_minus)
  }
  
  # Adjust the selected behaviors in both directions
  for (i in 1:300) {
    CompoIT_plus[i + 1, adjust_col1] <- CompoIT_plus[i + 1, adjust_col1] + i  # Add minutes to the first behavior
    CompoIT_plus[i + 1, adjust_col2] <- CompoIT_plus[i + 1, adjust_col2] - i  # Subtract minutes from the second behavior
    
    CompoIT_minus[i + 1, adjust_col1] <- CompoIT_minus[i + 1, adjust_col1] - i  # Subtract minutes from the first behavior
    CompoIT_minus[i + 1, adjust_col2] <- CompoIT_minus[i + 1, adjust_col2] + i  # Add minutes to the second behavior
  }
  
  # Combine the datasets from both directions
  CompoIT_combined <- rbind(CompoIT_minus, CompoIT_plus)
  CompoIT_combined_raw <- CompoIT_combined[!apply(CompoIT_combined < 0 | is.na(CompoIT_combined), 1, any), ]

  # Normalize back to proportions
  CompoIT_combined <- as.matrix(CompoIT_combined_raw / (60 * 24))
  CompoIT_combined <- as.data.frame(CompoIT_combined)
  
  # ILR transformation
  New_Z <- LR_ILR(CompoIT_combined, base_rot = 1)  # Assuming base rotation 1 for now

  # Convert ILR output to a data frame with the correct names
  colnames(New_Z) <- colnames(ILRModel1$model)[-1]  # Exclude the Outcome column
  colnames(CompoIT_combined_raw) <- colnames(ILRModel1$model)[-1]  # Exclude the Outcome column
  P2 <- predict(ILRModel1, New_Z, type = "link", se.fit = TRUE)
  
  # Calculate the 95% confidence intervals
  ci_lower <- P2$fit - 1.96 * P2$se.fit
  ci_upper <- P2$fit + 1.96 * P2$se.fit
  predictions <- exp(P2$fit) / (1 + exp(P2$fit))  # Convert log-odds to probabilities
  ci_lower_prob <- exp(ci_lower) / (1 + exp(ci_lower))
  ci_upper_prob <- exp(ci_upper) / (1 + exp(ci_upper))
  
  # Combine predictions and confidence intervals with input values
  results <- cbind(CompoIT_combined_raw, predictions, ci_lower_prob, ci_upper_prob)
  colnames(results) <- c(colnames(compo), "Prediction", "CI_Lower", "CI_Upper")
  return(results)
}

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


ITS_5_part_cph <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
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
  New_Z <- LR_ILR_5part(CompoIT_combined, base_rot = 1)
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



ITS_5_part_cph_baseR <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
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
  
  # Normalize back to proportions and convert to data frame
  CompoIT_combined <- as.matrix(CompoIT_combined_raw / (60 * 24))
  CompoIT_combined <- as.data.frame(CompoIT_combined)
  
  # ILR transformation
  New_Z <- as.data.frame(compositions::ilr(CompoIT_combined))
  colnames(New_Z) <- column_names
  colnames(CompoIT_combined_raw) <- raw_column_names
  
  # Define a list to store the predictions
  predictions_list <- list()
  mean_age <- mean(df$age, na.rm = TRUE)
  
  # Add fixed columns (sex and age) to New_Z
  New_Z$sex <- "Male"
  New_Z$age <- mean_age

  # Loop through each row in New_Z
  for (i in 1:nrow(New_Z)) {
    # Prepare the new_data for the current row
    new_data <- New_Z[i, ]
    
    # Call predict() using the new data with type = "lp" for the linear predictor
    P2 <- predict(ILRModel1, newdata = new_data, type = "lp", se.fit = TRUE)
   
    linear_predictor <- P2$linear.predictors[1]
    se_fit <- P2$se.fit[1]
    
    predictions_list[[i]] <- data.frame(
     Prediction = linear_predictor,  # Exponentiated linear predictor (HR)
      SE = se_fit              # Standard error of the linear predictor
    )
  }
  # Combine the predictions into a single data frame
  results <- do.call(rbind, predictions_list)
  results <- cbind(results, CompoIT_combined_raw)
  results <- results %>%
    mutate(
      CI_Lower = Prediction - 1.96 * SE,   # Lower 95% CI on the log scale
      CI_Upper = Prediction + 1.96 * SE,   # Upper 95% CI on the log scale
      yhat = exp(Prediction),    # Exponentiate the prediction
      lower = exp(CI_Lower),        # Exponentiate the lower CI
      upper = exp(CI_Upper)         # Exponentiate the upper CI
    )
  return(results)
}















ITS_5_part_cph <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
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
  New_Z <- LR_ILR_5part(CompoIT_combined, base_rot = 1)
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





ITS_5_part_cph_ilr <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
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
  CompoIT_combined
  
  # ILR transformation
  CompoIT_combined
  New_Z <- as.data.frame(compositions::ilr(CompoIT_combined))
  colnames(New_Z) <- column_names
  colnames(CompoIT_combined_raw) <- raw_column_names
  CompoIT_combined_raw
  New_Z
  
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




ITS_5_part <- function(df, ILRModel1, adjust_col1, adjust_col2) {
  # Use column indices for the four behaviors
  # Calculate the mean composition and convert it to minutes
  compomean <- as.vector(mean(acomp(compo)))
  compomeanmins <- as.vector(compomean * 60 * 24)
  
  # Convert to a data frame and create multiple rows for different adjustments (both + and -)
  CompoIT_plus <- as.data.frame(t(compomeanmins))
  CompoIT_minus <- as.data.frame(t(compomeanmins))
  duprow_plus <- CompoIT_plus[1,]
  duprow_minus <- CompoIT_minus[1,]
  
  # Create 301 rows for adjustments from -300 to +300
  for (i in 1:300) {
    CompoIT_plus <- rbind(CompoIT_plus, duprow_plus)
    CompoIT_minus <- rbind(CompoIT_minus, duprow_minus)
  }
  
  # Adjust the selected behaviors in both directions
  for (i in 1:300) {
    CompoIT_plus[i + 1, adjust_col1] <- CompoIT_plus[i + 1, adjust_col1] + i  # Add minutes to the first behavior
    CompoIT_plus[i + 1, adjust_col2] <- CompoIT_plus[i + 1, adjust_col2] - i  # Subtract minutes from the second behavior
    
    CompoIT_minus[i + 1, adjust_col1] <- CompoIT_minus[i + 1, adjust_col1] - i  # Subtract minutes from the first behavior
    CompoIT_minus[i + 1, adjust_col2] <- CompoIT_minus[i + 1, adjust_col2] + i  # Add minutes to the second behavior
  }
  
  # Combine the datasets from both directions
  CompoIT_combined <- rbind(CompoIT_minus, CompoIT_plus)
  CompoIT_combined_raw <- CompoIT_combined[!apply(CompoIT_combined < 0 | is.na(CompoIT_combined), 1, any), ]
  
  # Normalize back to proportions
  CompoIT_combined <- as.matrix(CompoIT_combined_raw / (60 * 24))
  CompoIT_combined <- as.data.frame(CompoIT_combined)
  
  # ILR transformation
  New_Z <- LR_ILR_5part(CompoIT_combined, base_rot = 1)  # Assuming base rotation 1 for now
  
  # Convert ILR output to a data frame with the correct names
  colnames(New_Z) <- colnames(ILRModel1$model)[-1]  # Exclude the Outcome column
  colnames(CompoIT_combined_raw) <- colnames(ILRModel1$model)[-1]  # Exclude the Outcome column
  P2 <- predict(ILRModel1, New_Z, type = "link", se.fit = TRUE)
  
  # Calculate the 95% confidence intervals
  ci_lower <- P2$fit - 1.96 * P2$se.fit
  ci_upper <- P2$fit + 1.96 * P2$se.fit
  predictions <- exp(P2$fit) / (1 + exp(P2$fit))  # Convert log-odds to probabilities
  ci_lower_prob <- exp(ci_lower) / (1 + exp(ci_lower))
  ci_upper_prob <- exp(ci_upper) / (1 + exp(ci_upper))
  
  # Combine predictions and confidence intervals with input values
  results <- cbind(CompoIT_combined_raw, predictions, ci_lower_prob, ci_upper_prob)
  colnames(results) <- c(colnames(compo), "Prediction", "CI_Lower", "CI_Upper")
  return(results)
}





ITS_4_part_cph <- function(df, ILRModel1, adjust_col1, adjust_col2, column_names, raw_column_names) {
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
  New_Z <- as.data.frame(compositions::ilr(CompoIT_combined))
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




