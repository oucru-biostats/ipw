#check if exposure is a binary variable of (0 and 1) or factor char variable 
check_exposure <- function(exposure) {
  # Check if exposure has exactly 2 unique values
  unique_vals <- unique(exposure)
  num_unique_vals <- length(unique_vals)
  
  if (num_unique_vals > 2) {
    stop(paste(
      "Error: The 'exposure' variable must be binary, with only two unique values.\n",
      "Currently, it has", num_unique_vals, "unique values:",
      paste(unique_vals, collapse = ", "), "."
    ))
  } else if (num_unique_vals == 1) {
    warning("Warning: The 'exposure' variable contains only one unique value.\nEnsure this is intentional.")
  }
  
  # Check for numeric type with values 0 and 1 only
  if (is.numeric(exposure)) {
    non_binary_vals <- unique_vals[!unique_vals %in% c(0, 1)]
    if (length(non_binary_vals) > 0) {
      stop(paste(
        "Error: Numeric 'exposure' values must be either 0 or 1.\n",
        "Found invalid values:",
        paste(non_binary_vals, collapse = ", "), "."
      ))
    }
  }
  
  # If non-numeric variable, check factor requirement
  else if (is.factor(exposure)) {
    # Extract factor levels
    factor_levels <- levels(exposure)
    
    if (length(factor_levels) != 2) {
      stop("Error: The factor 'exposure' must have exactly two levels.")
    }
    
    lower_level <- factor_levels[1]
    higher_level <- factor_levels[2]
    
    warning(paste0(
      "Notice: The factor levels of 'exposure' will be recoded.\n",
      "The lower level (", lower_level, ") will be assigned as unexposed, and the higher level (", higher_level, ") will be assigned as exposed."
    ))
    
    # Convert factor to 0 and 1 based on levels
    exposure <- as.numeric(exposure) - 1
  } else {
    stop("Error: non-numeric 'exposure' variable must be converted to a factor before use.")
  }
  
  # Return the modified exposure variable
  return(exposure)
}