rm(list = ls())

library(CATALYST)
library(flowCore)
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

# filenames <- list.files(InputDirectory, pattern = ".csv")
# filenames <- as_tibble(filenames) %>% rename(., filenames = value) %>%
#   mutate(Well = str_match(filenames, pattern = "Well_[:digit:][:digit:][:digit:]"))
# colnames(filenames)
# filenames
# 
# # # write_csv(as.data.frame(Filename), file = paste(MetaDirectory, "sample.details.csv", sep = "/"))
# # metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
# # sample_details <- read_csv(metadata_file)
# # #add cells per sample
# # sample_details$`Cells per sample` <- map_dbl(data.list, nrow)
# # sample_details
# # write_csv(sample_details, file = metadata_file)
# 
# metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
# sample_details <- read_csv(metadata_file)
# tmp <- sample_details$Filename[45]
# str_match(tmp, pattern = "Well_[:digit:][:digit:][:digit:]")
# 
# metadata <- sample_details %>%
#   mutate(Well = str_match(Filename, pattern = "Well_[:digit:][:digit:][:digit:]")) %>%
#   full_join(., filenames, by = c("Well")) %>%
#   select(filenames, MRN, Timepoint, Sample, Disease, Group)
# 
# metadata %>% write_csv(file = file.path(MetaDirectory, "metadata.csv"))
# metadata



metadata <- read_csv(file = file.path(MetaDirectory, "metadata.csv"))
metadata

sample_md <- metadata
sample_md

CSVfiles <- sample_md$filenames

# convert csv files to fcs
csvTofcs <- function(file.names, dest){
  # create an empty list to start
  DataList <- list()

  for(file in file.names){
    tmp <- read_csv(file.path(file))
    file <- gsub(".csv", "", file)
    DataList[[file]] <- tmp
  }
  rm(tmp)

  filenames <- names(DataList)
  head(DataList)

  # convert csv to fcs

  for(i in c(1:length(filenames))){
    data_subset <- DataList[i]
    data_subset <- data.table::rbindlist(as.list(data_subset))
    file_name <- names(DataList)[i]

    metadata <- data.frame(name = dimnames(data_subset)[[2]], desc = "")

    # create FCS file metadata
    # metadata$range <- apply(apply(data_subset, 2, range), 2, diff)
    metadata$minRange <- apply(data_subset, 2, min)
    metadata$maxRange <- apply(data_subset, 2, max)


    # data as matrix by exprs
    data_subset.ff <- new("flowFrame", exprs = as.matrix(data_subset),
                          parameters = AnnotatedDataFrame(metadata))

    head(data_subset.ff)
    write.FCS(data_subset.ff, paste0(dest, "/", file_name, ".fcs"), what = "numeric")
  }
}
setwd(InputDirectory)
# csvTofcs(CSVfiles, fcsDir)

# select samples to keep for analysis
sample_md
sample_md <- sample_md %>% 
  # dplyr::filter(Disease != "CLL") %>%
  # dplyr::filter(Disease != "HD") %>%
  # dplyr::filter(Disease != "PCL") %>%
  dplyr::filter(Disease == "AML") %>%
  dplyr::filter(Group != "HD") %>%
  dplyr::filter(Timepoint == "Baseline") %>%
  dplyr::filter(Group != "no_harvest")

fcsFiles <- list.files(path = fcsDir, pattern = ".fcs")
fcsFiles %in% gsub(".csv", ".fcs", sample_md$filenames)
fcsToLoad <- fcsFiles[fcsFiles %in% gsub(".csv", ".fcs", sample_md$filenames)]
# read fcs files as flowSet and add $CYT keyword
fs <- read.flowSet(files = fcsToLoad, path = fcsDir, truncate_max_range = FALSE)


# create tibble sample_md
sample_md
sample_md %>% rename(Filename = filenames)

sample_md <- sample_md %>% rename(Filename = filenames) %>%
  select(Filename, Sample, Disease, Group) %>%
  mutate(Filename = gsub(".csv", "", Filename)) %>%
  dplyr::filter(Filename %in% gsub(".fcs","", fcsToLoad)) %>%
  mutate(file_name = paste0(Filename, ".fcs")) %>%
  mutate(patient_id = Sample) %>%
  mutate(condition = Group) %>%
  mutate(sample_id = paste(Sample, Disease, sep = "_")) %>%
  select(file_name, patient_id, condition, sample_id) %>%
  as.data.frame()

# create tibble panel_md
fcs_colname <- colnames(fs)
fcs_colname
antigen <- fcs_colname
antigen[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[2])
fluorochrome <- fcs_colname
fluorochrome[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[1])
marker_class <- fcs_colname
marker_class[ c(1:10, 12, 33:36)] <- "none"
marker_class[c(11, 13:32)] <- "type"
panel_md <- as_tibble(cbind(fcs_colname, antigen, fluorochrome, marker_class))
as.data.frame(panel_md)

# use CATALYST::prepData() to create SCE object from flowSet
sce <- prepData(fs, panel_md, sample_md, FACS = TRUE)
assay(sce, "exprs") <- assay(sce, "counts")

subsetSCE <- function(x, n_cells){
  cs <- split(seq_len(ncol(x)), x$sample_id)
  cs <- unlist(lapply(cs, function(.) sample(., min(n_cells, length(.)))))
  x <- x[, cs]
  return(x)
}

sub_sce <- subsetSCE(sce, 15000)


pbMDS(sub_sce, dims = c(2,3), fun = "median", features = type_markers(sce),
      size_by = TRUE, by = "sample_id")

sce <- sub_sce

p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 9
p

# n_events <- min(n_cells(sce))
# n_events
# n_cells(sce)
# plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")



saveRDS(sce, file = file.path(OutputDirectory, "sce_Tcell.rds"))


