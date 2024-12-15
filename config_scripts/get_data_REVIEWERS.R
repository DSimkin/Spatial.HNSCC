if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)
source(here("config_scripts/restore_r_env_func.R"))
message("\n==== Installing all necessary packages (restoring from `renv.lock` file) ====\n")
Sys.sleep(5)
custom_restore()


# Specify the directory to store the untared unziped downloaded data
data_dir <- here("data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
file.rename(file.path(here(), list.files(here(), pattern = "^GSE268014")), file.path(here("data"), list.files(here(), pattern = "^GSE268014")))

# Extract files if compressed
if(grepl(".tar$", list.files(data_dir))) {
  untar(file.path(data_dir, grep(".tar$", list.files(data_dir), value = TRUE)), exdir = here("data"), extras = NULL)
}
if(any(grepl(".tar.gz$", list.files(here("data"))))) {
  tar_zip_files <- list.files(here("data"), pattern = "\\.tar.gz$", full.names = TRUE)
  lapply(tar_zip_files, function(x) untar(x, exdir = here("data"), extras = NULL))
} else if(any(grepl(".gz$", list.files(here("data"))))) {
  zip_files <- list.files(here("data"), pattern = "\\.gz$", full.names = TRUE)
  lapply(zip_files, R.utils::gunzip, remove = TRUE)
}

# Arrange files in sub-directorires according to their sample name
files <- list.files(here("data"), pattern = "^GSM8283.*_.*\\..*|^GSM8348.*_.*\\..*")
for (file in files) {
  # Extract the sample name (portion before the first underscore)
  sample <- sub("^[^_]*_([^_]*)_.*", "\\1", file)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(file.path(here("data"), sample))) {
    dir.create(file.path(here("data"), sample))
  }
  if (!dir.exists(file.path(here("data"), sample, "spatial"))) {
    dir.create(file.path(here("data"), sample, "spatial"))
  }
  
  # Clean file name from GSM and sample-name prefixes
  new_name <- sub("^GSM[0-9]+_", "", file)    # remove GEO prefix
  new_name <- sub("^P[0-9]{4}_", "", new_name)
  
  # Move expression files to the directory & move image, scales and coordinate files to a spatial sub-directory for each sample
  if(grepl(".csv$|.png$|.jpg$|.json$", file, ignore.case = TRUE)) {
    file.rename(file.path(here("data"), file), file.path(here("data"), sample, "spatial", new_name))
  } else {
    file.rename(file.path(here("data"), file), file.path(here("data"), sample, new_name))
  }
}

# Sort samples by tumor site (laryngeal, oral, oropharyngeal)
samp_param <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
if (!dir.exists(file.path(here("data"), "Laryngeal_CA"))) {
  dir.create(file.path(here("data"), "Laryngeal_CA"))
}
if (!dir.exists(file.path(here("data"), "Oral_CA"))) {
  dir.create(file.path(here("data"), "Oral_CA"))
}
if (!dir.exists(file.path(here("data"), "Oropharynx_CA"))) {
  dir.create(file.path(here("data"), "Oropharynx_CA"))
}
existing_samps <- list.files(here("data"))[grepl("^P[0-9]{4}", list.files(here("data")))]

for(samp in existing_samps) {
  if(!samp %in% samp_param$Sample) {
    unlink(file.path(here("data"), samp), recursive = TRUE)
  } else {
    file.rename(file.path(here("data"), samp),
                file.path(here("data"), samp_param$Site[samp_param$Sample == samp], samp))
  }
}

message("Data files successfully downloaded and oraganized in the project data directory: ", here("data"))
