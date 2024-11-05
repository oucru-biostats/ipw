check_exposure <- function(exposure, family) {
  unique_vals <- unique(exposure)
  num_unique_vals <- length(unique_vals)
  
  if (family == "binomial") {
    # Binomial check: binary values (0, 1) or two levels if factor
    if (num_unique_vals > 2) {
      stop(paste(
        "Error: For family 'binomial', the 'exposure' variable must be binary with only two unique values.\n",
        "Currently, it has", num_unique_vals, "unique values:", 
        paste(unique_vals, collapse = ", "), "."
      ))
    } else if (num_unique_vals == 1) {
      warning("Warning: The 'exposure' variable contains only one unique value.\nEnsure this is intentional.")
    }
    if (is.numeric(exposure)) {
      non_binary_vals <- unique_vals[!unique_vals %in% c(0, 1)]
      if (length(non_binary_vals) > 0) {
        stop(paste(
          "Error: Numeric 'exposure' values must be either 0 or 1.\n",
          "Found invalid values:", paste(non_binary_vals, collapse = ", "), "."
        ))
      }
    } else if (is.factor(exposure)) {
      if (length(levels(exposure)) != 2) {
        stop("Error: The factor 'exposure' must have exactly two levels.")
      }
      exposure <- as.numeric(exposure) - 1
    } else {
      stop("Error: non-numeric 'exposure' variable must be converted to a factor before use.")
    }
  }
  
  else if (family == "survival") {
    # Survival check: status indicator (0 for alive, 1 for dead, etc.)
    valid_status_values <- c(0, 1, TRUE, FALSE, 2, 3)
    
    if (!is.numeric(exposure) && !is.logical(exposure) && !is.factor(exposure)) {
      stop("Error: For family 'survival', 'exposure' should be a numeric, logical, or factor variable representing the status indicator.")
    }
    
    if (is.numeric(exposure)) {
      # Check for valid numeric status indicators
      non_valid_vals <- unique_vals[!unique_vals %in% valid_status_values]
      if (length(non_valid_vals) > 0) {
        stop(paste(
          "Error: For family 'survival', 'exposure' values must be 0, 1, TRUE, FALSE, or 2/3 for interval censoring.\n",
          "Found invalid values:", paste(non_valid_vals, collapse = ", "), "."
        ))
      }
    } else if (is.logical(exposure)) {
      # Convert logical TRUE/FALSE to numeric
      exposure <- as.numeric(exposure)
      valid_status_values <- as.numeric(valid_status_values)
      non_valid_vals <- unique(exposure[!exposure %in% valid_status_values])
      if (length(non_valid_vals) > 0) {
        stop("Error: For family 'survival', logical exposure must represent 0 (alive) or 1 (dead).")
      }
    } else if (is.factor(exposure)) {
      # Convert factor to numeric and check levels
      factor_levels <- levels(exposure)
      if (!all(factor_levels %in% as.character(valid_status_values))) {
        stop("Error: The factor 'exposure' must have levels corresponding to valid status indicators (0, 1, TRUE, FALSE, 2, 3).")
      }
      exposure <- as.numeric(exposure)
    }
  }
  
  else if (family == "multinomial") {
    # Multinomial check: categorical variable with more than 2 levels
    if (!is.factor(exposure)) {
      stop("Error: For family 'multinomial', 'exposure' should be a factor with more than two levels.")
    }
    if (num_unique_vals <= 2) {
      stop("Error: For family 'multinomial', 'exposure' must have more than two unique values.")
    }
  }
  
  else if (family == "ordinal") {
    # Ordinal check: ordered factor variable with more than 2 levels
    if (!is.ordered(exposure)) {
      stop("Error: For family 'ordinal', 'exposure' should be an ordered factor.")
    }
    if (num_unique_vals <= 2) {
      stop("Error: For family 'ordinal', 'exposure' must have more than two unique values.")
    }
  }
  
  else if (family == "gaussian") {
    # Gaussian check: continuous numeric variable
    if (!is.numeric(exposure)) {
      stop("Error: For family 'gaussian', 'exposure' should be a continuous numeric variable.")
    }
  }
  
  else {
    stop("Error: Unknown family. Choose one of 'binomial', 'survival', 'multinomial', 'ordinal', or 'gaussian'.")
  }
  
  # Return the modified or validated exposure variable
  return(exposure)
}
