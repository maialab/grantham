#' Convert one-letter amino acid codes to three-letter codes
#'
#' Converts amino acid one-letter abbreviations to three-letter codes, e.g., L
#' to `r as_three_letter("L")`. The accepted codes in the input include the 20
#' standard amino acids and also B (Asparagine or Aspartic acid), converted
#' to Asx, and Z (Glutamine or Glutamic acid) converted to Glx.
#'
#' @param x A character vector of one-letter amino acid codes, e.g. `"S"`,
#'   `"R"`, `"L"`, or `"P"`.
#'
#' @return A character vector of three-letter amino acid codes, e.g. `"Ser"`,
#'   `"Arg"`, `"Leu"`, or `"Pro"`.
#'
#' @md
#' @examples
#' # Convert S to Ser, R to Arg and P to Pro.
#' as_three_letter(c('S', 'R', 'P'))
#'
#' # The function `as_three_letter()` is case insensitive on the input but will
#' # always return the three-letter codes with the first letter in uppercase.
#' as_three_letter(c('S', 's', 'p', 'P'))
#'
#' # Convert also special case codes B (Asparagine or Aspartic acid) and Z
#' # (Glutamine or Glutamic acid)
#' as_three_letter(c('B', 'Z'))
#'
#' # Invalid codes in the input are converted to NA.
#' # "S" is correctly mapped to "Ser" but "Ser" and "Serine" are not
#' # one-letter amino acid codes and are therefore converted to NA.
#' as_three_letter(c('S', 's', 'Ser', 'Serine'))
#'
#' @export
as_three_letter <- function(x) {

  one_to_three_letter_codes <- stats::setNames(three_letter_codes, nm = one_letter_codes)
  unname(one_to_three_letter_codes[toupper(x)])
}
