LR_ILR_6part <- function(compo, base_rot=1) {
  D <- 6  # six-part composition
  
  
  ### Rotate through the basis
  if (base_rot == 1) {
    data <- compo
  } else if (base_rot == 2) {
    data <- cbind(compo[,2:6], compo[,1])
  } else if (base_rot == 3) {
    data <- cbind(compo[,3:6], compo[,1:2])
  } else if (base_rot == 4) {
    data <- cbind(compo[,4:6], compo[,1:3])
  } else if (base_rot == 5) {
    data <- cbind(compo[,5:6], compo[,1:4])
  } else if (base_rot == 6) {
    data <- cbind(compo[,6], compo[,1:5])
  }
  
  # Initialize z to have the correct number of columns for a 6-part composition
  z <- data[, 1:5]  # Initialize with 5 columns for 6-part composition
  
  # Compute the isometric log-ratio transformations
  z[,1] <- sqrt((D - 1) / D) * log(data[,1] / (data[,2] * data[,3] * data[,4] * data[,5] * data[,6])^(1 / (D - 1)))
  z[,2] <- sqrt((D - 2) / (D - 1)) * log(data[,2] / (data[,3] * data[,4] * data[,5] * data[,6])^(1 / (D - 2)))
  z[,3] <- sqrt((D - 3) / (D - 2)) * log(data[,3] / (data[,4] * data[,5] * data[,6])^(1 / (D - 3)))
  z[,4] <- sqrt((D - 4) / (D - 3)) * log(data[,4] / (data[,5] * data[,6])^(1 / (D - 4)))
  z[,5] <- sqrt((D - 5) / (D - 4)) * log(data[,5] / (data[,6])^(1 / (D - 5)))
  
  return(z)
}
