#' Convert three-letter amino acid codes to one-letter codes
#'
#' Converts three-letter amino acid abbreviations to one-letter codes, e.g., Leu
#' to `r as_one_letter("Leu")`. The accepted codes in the input include the 20
#' standard amino acids and also Asx (Asparagine or Aspartic acid), converted
#' to B, and Glx (Glutamine or Glutamic acid) converted to Z.
#'
#' @param x A character vector of three-letter amino acid codes, e.g. `"Ser"`,
#'   `"Arg"`, `"Leu"`, or `"Asx"`.
#'
#' @return A character vector of one-letter amino acid codes, e.g. `"S"`, `"R"`,
#'   `"L"`, or `"B"`.
#'
#' @md
#' @examples
#' # Convert Ser to S, Arg to R and Pro to P.
#' as_one_letter(c('Ser', 'Arg', 'Pro'))
#'
#' # The function `as_one_letter()` is case insensitive on the input but will
#' # always return the one-letter codes in uppercase.
#' as_one_letter(c('ser', 'ArG', 'PRO'))
#'
#' # Convert the codes of the 20 standard amino acids. Note that the function
#' # `amino_acids()` returns the three-letter codes of the 20 standard amino
#' # acids.
#' as_one_letter(amino_acids())
#'
#' # Convert also special case codes Asx (Asparagine or Aspartic acid) and Glx
#' # (Glutamine or Glutamic acid)
#' as_one_letter(c('Asx', 'Glx'))
#'
#' # Invalid codes in the input are converted to NA.
#' # "Ser" is correctly mapped to "S" but "Serine" is not as it is not a
#' # three-letter amino acid code (the same applies to "Glucose").
#' as_one_letter(c('Ser', 'Serine', 'Glucose'))
#'
#' @export
as_one_letter <- function(x) {

  three_to_one_letter_codes <- stats::setNames(one_letter_codes, nm = three_letter_codes)
  unname(three_to_one_letter_codes[stringr::str_to_title(x)])
}
