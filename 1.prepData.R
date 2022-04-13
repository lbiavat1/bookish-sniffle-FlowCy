rm(list = ls())

library(CATALYST)
library(tidyverse)
library(diffcyt)

### Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory
setwd(PrimaryDirectory)
if(!dir.exists("data"))
  dir.create("data")
setwd("data/")
InputDirectory <- getwd()
if(!dir.exists("fcs"))
  dir.create("fcs")
setwd("fcs/")
fcsDir <- getwd()
fcsDir
setwd(PrimaryDirectory)

### Set 'metadata' directory
setwd(PrimaryDirectory)
if(!dir.exists("metadata"))
  dir.create("metadata")
setwd("metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory
if(!dir.exists("output"))
  dir.create("output", showWarnings = FALSE)
setwd("output")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)


# prep and read in metadata

filenames <- list.files(InputDirectory, pattern = ".csv")
filenames <- as_tibble(filenames) %>% rename(., filenames = value) %>%
  mutate(Well = str_match(filenames, pattern = "Well_[:digit:][:digit:][:digit:]")) 
colnames(filenames)
filenames

# # write_csv(as.data.frame(Filename), file = paste(MetaDirectory, "sample.details.csv", sep = "/"))
# metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
# sample_details <- read_csv(metadata_file)
# #add cells per sample
# sample_details$`Cells per sample` <- map_dbl(data.list, nrow)
# sample_details
# write_csv(sample_details, file = metadata_file)

metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
sample_details <- read_csv(metadata_file)
tmp <- sample_details$Filename[45]
str_match(tmp, pattern = "Well_[:digit:][:digit:][:digit:]")

metadata <- sample_details %>% 
  mutate(Well = str_match(Filename, pattern = "Well_[:digit:][:digit:][:digit:]")) %>% 
  full_join(., filenames, by = c("Well"))

