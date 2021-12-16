#' @keywords internal
is_amino_acid <- function(x) {
  x %in% amino_acids()
}

#' @keywords internal
all_amino_acids <- function(x) {
  all(is_amino_acid(x))
}
