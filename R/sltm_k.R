#' Linear positions of the entries of a strictly lower triangular matrix
#'
#' Returns the linear indices of the non-zero entries of a strictly lower
#' triangular matrix.
#'
#' @param n Dimension of a `n` by `n` square matrix.
#'
#' @return An integer vector of linear positions in column-major order.
#' @md
#'
#' @examples
#' sltm_k(3)
#'
#' @noRd
#' @keywords internal
sltm_k <- function(n) {
  if(!(n > 1)) stop('`n` must be greater than 1')

  utils::combn(seq_len(n), 2, function(ij) {ij2k(i = ij[2], j = ij[1], n)})
}
