#' Read HLA locus alignment statistics from arcasHLA
#'
#' @param dir Path to arcasHLA output directory
#'
#' @return data frame
#' @export
#'
#' @examples
#' # No example
read_arcasHLA_reads <- function(dir){
  arg_col <- checkmate::makeAssertCollection()
  checkmate::assertDirectoryExists(dir)

  files <- list.files(dir, pattern = ".genotype.log")
  if (dir.exists(dir) & length(files)==0) {arg_col$push("No arcasHLA genotype.log files detected")}
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

  tibble::tibble(file = files) %>%
    dplyr::mutate(data = purrr::map(file, function(f){
      sample_path <- sprintf("%s/%s", dir, f)
      df <- suppressMessages(data.table::fread(sample_path, sep = "?", col.names = "lines"))
      if (any(grepl("error", df$lines, ignore.case = T))){ return(NA) }
      else {
        df <- df %>%
          dplyr::mutate(lines = gsub("\t", "", lines)) %>%
          dplyr::filter(grepl("^HLA", lines)) %>%
          tidyr::separate(lines, into = c("locus", "abundance", "reads", "classes"), sep = " +")
        return(df)
      }
      return(sample_path)
    })) %>%
    tidyr::unnest(data) %>%
    dplyr::select(!(dplyr::contains("data"))) %>%
    tidyr::separate(locus, into = c(NA, "locus"), sep = "-") %>%
    dplyr::mutate_at(c("reads", "classes"), as.numeric) %>%
    dplyr::mutate(file = gsub("\\.genotype\\.log", "", file)) %>%
    dplyr::rename(sample = file)
}

# read_arcasHLA_reads("/labs/khatrilab/solomonb/covid/isb/arcasHLA")
# "/labs/khatrilab/solomonb/covid/isb/arcasHLA/INCOV001-CV.genotype.log"

#' Create DRB345 read ratios
#'
#' @param drb A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "DRB1", "DRB3", "DRB4", "DRB5" - raw read counts for each loci
#'
#' @return A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "DRB3", "DRB4", "DRB5" - ratio of read counts for each loci to those of DRB1
#'
#' @examples
#' data(drb345_data)
#' hlatools:::make_drb345_ratios(drb345_data)
#'
make_drb345_ratios <- function(drb){
  # Assert compatible columns in drb
  arg_col <- checkmate::makeAssertCollection()
  drb_columns <- c("sample", "DRB1", "DRB3", "DRB4", "DRB5")
  checkmate::assertNames(names(drb), must.include = drb_columns, .var.name = "drb column names", add = arg_col)
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print); checkmate::reportAssertions(arg_col)}

  # Calculate ratios
  drb %>%
    dplyr::mutate_at(dplyr::vars(DRB3:DRB5), ~ ./DRB1) %>%
    dplyr::select(-DRB1)
}

#' Predict DRB345 copy numbers
#'
#' @param drb A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "DRB1", "DRB3", "DRB4", "DRB5" - raw read counts for each loci
#'
#' @return A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "DRB3", "DRB4", "DRB5" - copy number of each loci represented as value of 0,1,2
#'
#' @import kknn
#'
#' @export
#'
#' @examples
#' data(drb345_data)
#' predict_drb345(drb345_data)
#'
predict_drb345 <- function(drb) {
  # Assertions for drb are contained in make_drb345_ratios

  # Convert DRB345 counts in ratios
  drb <- make_drb345_ratios(drb)

  # drb345_knn model is available from sysdata.rda

  # Format drb and apply kNN
  drb <- tidyr::expand_grid(drb, locus = c("DRB3", "DRB4", "DRB5")) %>%
    dplyr::mutate(locus = factor(locus))

  # Critical that kknn sourced from github pull #24
  # This updates the contrasts function
  # Within kknn::kknn,
  #   originally, c(unordered = "contr.dummy", ordered = "contr.ordinal")
  #   from pull, contrasts=c(unordered = kknn::contr.dummy, ordered = kknn::contr.ordinal)

  knn_prediction <- suppressWarnings(stats::predict(drb345_knn$fit, drb))

  drb %>%
    dplyr::mutate("copy_number" = as.numeric(as.character(knn_prediction))) %>%
    dplyr::select(sample, locus, copy_number) %>%
    tidyr::pivot_wider(names_from = "locus", values_from = "copy_number")

}
