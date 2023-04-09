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
    dplyr::mutate_at(dplyr::vars(DRB3:DRB5), dplyr::funs(./DRB1)) %>%
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
#' @export
#'
#' @examples
#' # data(drb345_data)
#' # predict_drb345(drb345_data)
#'
predict_drb345 <- function(drb) {
  # Assertions for drb are contained in make_drb345_ratios

  # Convert DRB345 counts in ratios
  drb <- make_drb345_ratios(drb)

  # drb345_knn model is available from sysdata.rda

  # Format drb and apply kNN
  drb <- tidyr::expand_grid(drb, locus = c("DRB3", "DRB4", "DRB5"))

  drb %>%
    dplyr::bind_cols(stats::predict(drb345_knn, new_data = drb)) %>%
    dplyr::mutate("copy_number" = as.numeric(as.character(.pred_class))) %>%
    dplyr::select(sample, locus, copy_number) %>%
    tidyr::pivot_wider(names_from = "locus", values_from = "copy_number")
}
