#' Generate amino acid pairs
#'
#' This function generates combinations of amino acids in pairs. By default, it
#' generates all pair combinations of the 20 standard amino acids.
#'
#' @param x A character vector of amino acids (three-letter codes).
#' @param y Another character vector of amino acids (three-letter codes).
#' @param keep_self Whether to keep pairs involving the same amino acid.
#' @param keep_duplicates Whether to keep duplicated pairs.
#' @param keep_reverses Whether to keep pairs that are reversed versions of
#'   others. E.g. if `keep_reverses` is `TRUE` the pairs `"Ser"`-`"Arg"` and
#'   `"Arg"`-`"Ser"` will be kept in the returned tibble; however, if
#'   `keep_reverses` is `FALSE`, only the first pair is preserved in the output.
#'
#' @return A [tibble][tibble::tibble-package] of amino acid pairs.
#'
#' @examples
#' # Generate all pairs of the 20 standard amino acids
#' amino_acid_pairs()
#'
#' # Remove the self-to-self pairs
#' amino_acid_pairs(keep_self = FALSE)
#'
#' # Generate specific combinations of Ser against Ala and Trp.
#' amino_acid_pairs(x = 'Ser', y = c('Ala', 'Trp'))
#' @md
#' @importFrom rlang .data
#' @export
amino_acid_pairs <-
  function(x = amino_acids(),
           y = amino_acids(),
           keep_self = TRUE,
           keep_duplicates = TRUE,
           keep_reverses = TRUE) {

    if(!all_amino_acids(x))
      stop('`x` must be a vector of three-letter code amino acids')

    if (!all_amino_acids(y))
      stop('`y` must be a vector of three-letter code amino acids'
      )

  tbl <- tidyr::expand_grid(x = x, y = y)
  tbl <- `if`(keep_self, tbl, dplyr::filter(tbl, x != y))
  tbl <- `if`(keep_duplicates, tbl, dplyr::distinct(tbl))

  tbl <-
    if (keep_reverses) {
      tbl # do nothing
    } else {
      tbl %>%
        dplyr::rowwise() %>%
        dplyr::mutate(key = paste(sort(c(x, y)), collapse = '-')) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.data$key, .keep_all = TRUE) %>%
        dplyr::select(-'key')
    }

  return(tbl)
}
