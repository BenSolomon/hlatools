#' Convert genotype predictions into tree input format
#' Essentially one-hot encoding of the genotype column
#'
#' @param hla A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "locus" - one of the following HLA loci: "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5" \cr
#' "field" - the resolution of the hla allele in the format of "field_1", "field_2", "field_3" \cr
#' "genotyper" - one of the following genotyping algorithms "arcasHLA", "hlaminer", "optitype", "phlat" \cr
#' "allele" - a character representing the predicted allele with fields separated by "_" (e.g. "01", "01_02", "01_02_01") \cr
#'
#' @return A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "locus" - one of the following HLA loci: "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5" \cr
#' "field" - the resolution of the hla allele in the format of "field_1", "field_2", "field_3" \cr
#' "present_arcasHLA", "present_hlaminer", "present_optitype", "present_phlat" - a T/F value indicating if a prediction was present for the indicated genotyping algorithm
#'
#' @examples
#' data(composite_data)
#' hlatools:::make_tree_input(composite_data)
make_tree_input <- function(hla){
  # Assert compatible columns and values in hla
  arg_col <- checkmate::makeAssertCollection()
  accept_column <- c("sample", "locus", "field", "genotyper", "allele")
  accept_locus <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
  accept_field <- c("field_1", "field_2", "field_3")
  accept_genotyper <- c("arcasHLA", "hlaminer", "optitype", "phlat")
  checkmate::assertNames(names(hla), must.include = accept_column, .var.name = "Data frame columns", add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}
  checkmate::assertNames(unique(hla$locus), subset.of = accept_locus, .var.name = "Locus column", add = arg_col)
  checkmate::assertNames(unique(hla$field), subset.of = accept_field, .var.name = "Field column", add = arg_col)
  checkmate::assertNames(unique(hla$genotyper), subset.of = accept_genotyper, .var.name = "Genotyper column", add = arg_col)
  if (arg_col$isEmpty()==F) {map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

  # Format hla for tree model
  hla %>%
    dplyr::select(sample, locus, field, genotyper, allele) %>%
    tidyr::nest(allele = allele) %>%
    dplyr::mutate(present = purrr::map_lgl(allele, function(x) !all(is.na(x)))) %>% # ID predictions that only generated NAs
    dplyr::select(sample, locus, field, genotyper, present) %>%
    tidyr::pivot_wider(names_from = "genotyper", values_from = "present", names_prefix = "present_", values_fill = FALSE) %>%
    dplyr::mutate_if(is.character, factor) %>%
    dplyr::mutate_if(is.logical, factor) %>%
    dplyr::mutate(sample = paste(sample, 1:dplyr::n(), sep = "_")) %>%
    tibble::column_to_rownames("sample")
}

#' Predict the best genotyping algorithm for each allele prediction
#'
#' @param hla A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "locus" - one of the following HLA loci: "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5" \cr
#' "field" - the resolution of the hla allele in the format of "field_1", "field_2", "field_3" \cr
#' "genotyper" - one of the following genotyping algorithms "arcasHLA", "hlaminer", "optitype", "phlat" \cr
#' "allele" - a character representing the predicted allele with fields separated by "_" (e.g. "01", "01_02", "01_02_01") \cr
#'
#' @return A data frame with columns: \cr
#' "sample" - unique sample name \cr
#' "locus" - one of the following HLA loci: "A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5" \cr
#' "field" - the resolution of the hla allele in the format of "field_1", "field_2", "field_3" \cr
#' "best_genotyper" - the best genotyping out of "arcasHLA", "hlaminer", "optitype", "phlat" as determined by the tree model
#'
#' @export
#'
#' @examples
#' data(composite_data)
#' predict_composite(composite_data)
predict_composite <- function(hla){
  hla <- make_tree_input(hla)

  # hla_tree model is available from sysdata.rda

  hla %>%
    dplyr::select(locus:field) %>%
    dplyr::bind_cols(parsnip::predict.model_fit(hla_tree$model, new_data = hla)) %>%
    dplyr::rename("best_genotyper" = .pred_class) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(sample = gsub("_[0-9]*$", "", sample))
}
