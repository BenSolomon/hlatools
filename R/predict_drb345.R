#' Create DRB345 read ratios
#'
#' @param df A data frame with columns: \cr
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
make_drb345_ratios <- function(df){
  # Assert compatible columns in df
  arg_col <- checkmate::makeAssertCollection()
  df_columns <- c("sample", "DRB1", "DRB3", "DRB4", "DRB5")
  checkmate::assertNames(names(df), must.include = df_columns, .var.name = "df column names", add = arg_col)
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print); checkmate::reportAssertions(arg_col)}

  # Calculate ratios
  df %>%
    dplyr::mutate_at(dplyr::vars(DRB3:DRB5), dplyr::funs(./DRB1)) %>%
    dplyr::select(-DRB1)
}



#' Predict DRB345 copy numbers
#'
#' @param df A data frame with columns: \cr
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
predict_drb345 <- function(df) {
  # Assertions for df are contained in make_drb345_ratios

  # Convert DRB345 counts in ratios
  df <- make_drb345_ratios(df)

  # drb345_knn model is available from sysdata.rda

  # Format df and apply kNN
  df <- tidyr::expand_grid(df, locus = c("DRB3", "DRB4", "DRB5")) %>%

  df %>%
    dplyr::bind_cols(stats::predict(drb345_knn, new_data = df)) %>%
    dplyr::mutate("copy_number" = as.numeric(as.character(.pred_class))) %>%
    dplyr::select(sample, locus, copy_number) %>%
    tidyr::pivot_wider(names_from = "locus", values_from = "copy_number")
}
