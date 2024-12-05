custom_restore <- function() {
  # Load required libraries
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")  # Install renv if not already installed
  }
  library(renv)
  
  # Restore environment
  errors <- list()  # Initialize an empty list to store errors
  
  tryCatch({
    # Attempt to restore the environment
    renv::restore(prompt = FALSE)
  }, error = function(e) {
    # Handle any errors during restore
    errors <<- append(errors, list(e$message))
  })
  
  # Handle individual package errors
  missing_packages <- setdiff(
    renv::dependencies()$Package,
    installed.packages()[, "Package"]
  )
  
  if (length(missing_packages) > 0) {
    cat("Attempting to install missing packages individually...\n")
    for (pkg in missing_packages) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        errors <<- append(errors, paste("Failed to install", pkg, ":", e$message))
      })
    }
  }
  
  # Display errors, if any
  if (length(errors) > 0) {
    cat("\nThe following errors occurred during restoration:\n")
    cat(paste(errors, collapse = "\n"))
  } else {
    cat("\nEnvironment successfully restored!\n")
  }
}
