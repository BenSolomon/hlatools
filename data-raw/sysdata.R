## code to prepare `sysdata` dataset goes here

drb345_knn <- readRDS("data-raw/drb345_knn.RDS")
hla_tree <- readRDS("data-raw/hla_tree.RDS")

usethis::use_data(drb345_knn, hla_tree, overwrite = TRUE, internal = T)
