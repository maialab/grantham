#' The 20 standard amino acids
#'
#' The 20 amino acids that are encoded directly by the codons of the universal
#' genetic code.
#'
#' @return Three-letter codes of the standard amino acids.
#'
#' @examples
#' amino_acids()
#'
#' @export
amino_acids <- function() {
  c("Ser", "Arg", "Leu", "Pro", "Thr", "Ala", "Val", "Gly", "Ile",
    "Phe", "Tyr", "Cys", "His", "Gln", "Asn", "Lys", "Asp", "Glu",
    "Met", "Trp")
}
