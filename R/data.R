#' Counts of reads aligned to the human DRB loci for 346 subjects
#'
#' @format A data frame with 346 rows and 5 columns:
#' \describe{
#'   \item{sample}{A unique sample name}
#'   \item{DRB1, DRB3, DRB4, DRB5}{Raw counts of reads aligned to each locus}
#' }
#' @source https://github.com/BenSolomon/hla_benchmark
"drb345_data"

#' Genotype predictions from multiple algorithms
#'
#' @format A data frame with 942 rows and 5 columns:
#' \describe{
#'   \item{sample}{A unique sample name}
#'   \item{locus}{An HLA gene locus including "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5"}
#'   \item{field}{An HLA genotype field designation including "field_1", "field_2", "field_3"}
#'   \item{genotyper}{An HLA genotyping algorithm including "arcasHLA", "hlaminer", "optitype", "phlat"}
#'   \item{allele}{An HLA allele prediction, with fields separated by "_" (e.g. "26_01" or "24" or "01_03_01")}
#' }
#' @source https://github.com/BenSolomon/hla_benchmark
"composite_data"
