#' Read genotype output from arcasHLA
#'
#' @param dir Path to arcasHLA output directory
#'
#' @return A data frame
#' @export
#'
#' @examples
#' # No example
read_arcasHLA <- function(dir){
    arg_col <- checkmate::makeAssertCollection()
    run_arcas_message <- sprintf("Valid directory, but no genotype.tsv file.
  If directory is correct, consider arcashHLA merge in shell:
  \'ARCAS_DIR=%s;
  arcasHLA merge --i $ARCAS_DIR --o $ARCAS_DIR\'",dir)
    path <- sprintf("%s/genotypes.tsv", dir)

    checkmate::assertDirectoryExists(dir)
    if (dir.exists(dir) & !file.exists(path)) {arg_col$push(run_arcas_message)}
    if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

    suppressMessages(readr::read_tsv(path)) %>%
      tidyr::pivot_longer(-subject, names_to = "allele_id", values_to = "allele") %>%
      tidyr::drop_na() %>%
      dplyr::mutate(allele_id = stringr::str_sub(allele_id, -1),
                    genotyper = "arcasHLA") %>%
      dplyr::rename(sample = subject) %>%
      tidyr::separate(allele, into = c("locus"), sep = "\\*", remove = F, extra = "drop") %>%
      dplyr::select(sample, locus, allele_id, allele, genotyper)
}
# read_arcasHLA("/labs/khatrilab/solomonb/covid/isb/arcasHLA")


#' Read genotype output from OptiType
#'
#' @param dir Path to OptiType output directory
#'
#' @return data frame
#' @export
#'
#' @examples
# No example
read_OptiType <- function(dir){
  arg_col <- checkmate::makeAssertCollection()
  checkmate::assertDirectoryExists(dir)

  files <- list.files(dir, pattern = "_result.tsv")
  if (dir.exists(dir) & length(files)==0) {arg_col$push("No OptiType results.tsv files detected")}
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

  tibble::tibble(file = files) %>%
    dplyr::mutate(data = purrr::map(file, function(f){
      sample_path <- sprintf("%s/%s", dir, f)
      # suppressMessages(readr::read_tsv(sample_path, )) %>%
      suppressMessages(data.table::fread(sample_path)) %>%
        dplyr::select(A1:C2) %>%
        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "allele_id", values_to = "allele")%>%
        dplyr::mutate(allele_id = stringr::str_sub(allele_id, -1),
               genotyper = "optitype")
    })) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(file = gsub("_result\\.tsv", "", file)) %>%
    dplyr::rename(sample = file) %>%
    tidyr::separate(allele, into = c("locus"), sep = "\\*", remove = F, extra = "drop") %>%
    dplyr::select(sample, locus, allele_id, allele, genotyper)
}
# read_OptiType("/labs/khatrilab/solomonb/covid/isb/optitype")


#' Read genotype output from PHLAT
#'
#' @param dir Path to PHLAT output directory
#'
#' @return data frame
#' @export
#'
#' @examples
#' # No example
read_PHLAT <- function(dir){
  arg_col <- checkmate::makeAssertCollection()
  checkmate::assertDirectoryExists(dir)

  files <- list.files(dir, pattern = "_HLA.sum")
  if (dir.exists(dir) & length(files)==0) {arg_col$push("No PHLAT HLA.sum files detected")}
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

  tibble::tibble(file = files) %>%
    dplyr::mutate(data = purrr::map(file, function(f){
      sample_path <- sprintf("%s/%s", dir, f)
      # suppressMessages(readr::read_tsv(sample_path, )) %>%
      suppressMessages(data.table::fread(sample_path)) %>%
        dplyr::select(dplyr::contains("Allele")) %>%
        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "allele_id", values_to = "allele")%>%
        dplyr::mutate(allele_id = stringr::str_sub(allele_id, -1),
                      genotyper = "phlat")
    })) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(file = gsub("_HLA\\.sum", "", file)) %>%
    dplyr::rename(sample = file) %>%
    tidyr::separate(allele, into = c("locus"), sep = "\\*", remove = F, extra = "drop") %>%
    dplyr::select(sample, locus, allele_id, allele, genotyper)
}
# read_PHLAT("/labs/khatrilab/solomonb/covid/isb/phlat")


#' Read HLAminer output for a single file
#' Separated single file read into a single function since parsing
#' HLAminer output is complicated and would be difficult to read inside
#' of the function to parse the entire directory.
#' Meant to be used internally.
#'
#' @param path Path to HLAminer output directory
#'
#' @return data frame
#'
#' @examples
#' # No example
read_HLAminer_single_sample <- function(path){
  if (file.exists(path)){
    # tidyr::tibble(lines = readr::read_lines(path)) %>%
    data.table::fread(path, sep="?", col.names = c("lines")) %>%
      dplyr::filter(grepl("^HLA|\t", lines)) %>%
      tidyr::separate(lines, into = c("hla", "prediction", "allele"), sep = "\t",
               fill = "right", extra = "drop") %>%
      dplyr::mutate_all(function(x) ifelse(x=="",NA,x)) %>%
      tidyr::fill(hla, prediction, .direction = "down") %>%
      tidyr::drop_na() %>%
      tidyr::separate(allele, into = c("allele", "score", "expected", "confidence"),
               sep = ",", fill = "right", extra = ) %>%
      dplyr::group_by(hla, prediction) %>%
      dplyr::mutate_at(dplyr::vars(c(score, expected, confidence)), as.numeric) %>%
      # Clean up names
      dplyr::mutate(allele_id = purrr::map_chr(prediction, ~stringr::str_extract(., "[1-9]"))) %>%
      dplyr::mutate(allele = gsub("[A-Z]$", "", allele)) %>%
      # Select highest score
      dplyr::group_by(hla, allele_id) %>%
      dplyr::top_n(1, score) %>%
      dplyr::ungroup() %>%
      dplyr::select(allele_id, allele) %>%
      dplyr::mutate(genotyper = "hlaminer")
  } else {
    NULL
  }
}
# read_HLAminer_single_sample("/labs/khatrilab/solomonb/covid/isb/hla_miner/HLAminer_HPRA_INCOV001-CV.csv")


#' Read genotype output from PHLAT
#'
#' @param dir File path to a single HLAminer output file
#'
#' @return data frame
#' @export
#'
#' @examples
#' # No examples
read_HLAminer <- function(dir){
  arg_col <- checkmate::makeAssertCollection()
  checkmate::assertDirectoryExists(dir)

  files <- list.files(dir, pattern = "HLAminer_HPRA_.*\\.csv")
  if (dir.exists(dir) & length(files)==0) {arg_col$push("No HLAminer_HPRA csv files detected")}
  if (arg_col$isEmpty()==F) {purrr::map(arg_col$getMessages(),print);checkmate::reportAssertions(arg_col)}

  tibble::tibble(file = files) %>%
    dplyr::mutate(data = purrr::map(file, function(f){
      sample_path <- sprintf("%s/%s", dir, f)
      read_HLAminer_single_sample(sample_path)
    })) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(file = gsub("(HLAminer_HPRA_|\\.csv)", "", file)) %>%
    dplyr::rename(sample = file) %>%
    tidyr::separate(allele, into = c("locus"), sep = "\\*", remove = F, extra = "drop") %>%
    dplyr::select(sample, locus, allele_id, allele, genotyper)
}
# read_HLAminer("/labs/khatrilab/solomonb/covid/isb/hla_miner")

# "/labs/khatrilab/solomonb/covid/isb"

