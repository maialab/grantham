#' Grantham distance
#'
#' @description
#' This function calculates Grantham's distance \eqn{d_{i,j}} between two
#' amino acids (\eqn{i} and \eqn{j}) based on their chemical properties:
#'
#' \deqn{d_{i,j} = \rho ((\alpha (c_i-c_j)^2 + \beta (p_i-p_j)^2 + \gamma (v_i-v_j)^2)^\frac{1}{2}}
#'
#' This calculation is based on three amino acid side chain properties that were
#' found to be the three strongest correlators with the relative substitution
#' frequency (RSF) (references cited in Grantham (1974)), namely:
#'
#' - composition \eqn{c}, meaning the atomic weight ratio of hetero (noncarbon)
#' elements in end groups or rings to carbons in the side chain.
#' - polarity \eqn{p};
#' - molecular volume \eqn{v}.
#'
#' Each property difference is weighted by dividing by the mean distance found
#' with it alone in the formula. The constants \eqn{\alpha}, \eqn{\beta} and
#' \eqn{\gamma} are squares of the inverses of mean distances of each property,
#' respectively.
#'
#' The distances reported by Grantham (1972) are further scaled by a factor
#' ---here coined \eqn{\rho}--- such that the mean of all distances is 100.
#' Although this factor is not explicitly included in Grantham's distance
#' formula, it is actually used for calculating the amino acid pair distances
#' reported in Table 2 of Grantham's paper. So, for all intents and purposes,
#' this factor should be regarded as part of the formula used to calculate
#' Grantham distance, and therefore we include it explicitly in the equation
#' above.
#'
#' If you want to calculate Grantham's distance right off from the identity of
#' the amino acids, instead of using their chemical properties, then use
#' [grantham_distance()].
#'
#' @param c_i composition value for the _ith_ amino acid.
#' @param c_j composition value for the _jth_ amino acid.
#' @param p_i polarity value for the _ith_ amino acid.
#' @param p_j polarity value for the _jth_ amino acid.
#' @param v_i molecular volume value for the _ith_ amino acid.
#' @param v_j molecular volume value for the _jth_ amino acid.
#' @param alpha The constant \eqn{\alpha} in the equation of Grantham's
#'   paper, in page 863.
#' @param beta The constant \eqn{\beta} in the equation of Grantham's
#'   paper, in page 863.
#' @param gamma The constant \eqn{\gamma} in the equation of Grantham's
#'   paper, in page 863.
#' @param rho Grantham's distances reported in Table 2, Science (1974).
#'   185(4154): 862--4 by R. Grantham, are scaled by a factor (here named
#'   \eqn{\rho}) such that the mean value of all distances are 100. The `rho`
#'   parameter allows this factor \eqn{\rho} to be changed. By default
#'   \eqn{\rho=50.723}, the same value used by Grantham. This value is
#'   originally mentioned in the caption of Table 2 of the aforementioned paper.
#'
#' @return A double vector of Grantham's distances.
#'
#' @seealso Check [amino_acids_properties] for a table of the three property
#'   values that can be used with this formula. This data set is from Table 1,
#'   Science (1974). 185(4154): 862--4 by R. Grantham.
#'
#' @md
#' @export
grantham_equation <-
  function(c_i,
           c_j,
           p_i,
           p_j,
           v_i,
           v_j,
           alpha = 1.833,
           beta = 0.1018,
           gamma = 0.000399,
           rho = 50.723) {

    d_ij <- rho *
      (alpha * (c_i - c_j) ^ 2 +
         beta * (p_i - p_j) ^ 2 +
         gamma * (v_i - v_j) ^ 2) ^ 0.5

    return(d_ij)
  }

#' Grantham distance
#'
#' @description
#' This function calculates the Grantham distance for pairs of amino acids.
#' Amino acid identities should be provided as three-letter codes in `x` and
#' `y`. Amino acids identified in `x` and `y` are matched element-wise, i.e. the
#' first element of `x` is paired with the first element of `y`, and so on.
#'
#' The Grantham distance attempts to provide a proxy for the evolutionary
#' distance between two amino acids based on three key chemical properties:
#' composition, polarity and molecular volume. In turn, evolutionary distance is
#' used as a proxy for the impact of missense substitutions. The higher the
#' distance, the more deleterious the substitution is.
#'
#' The distance calculation is provided by two methods. The so-called _original_
#' method, meaning that the amino acid distances used are the ones provided by
#' Grantham in his original publication in Table 2. This is the default method.
#' In addition, you may choose the _exact_ method, which uses the chemical
#' properties provided in Grantham's Table 1 to compute the amino acid
#' differences anew. The distances calculated with the _exact_ method are not
#' rounded to the nearest integer and will differ by ~1 unit for some amino acid
#' pairs from the _original_ method.
#'
#' If you want to calculate Grantham's distance by providing the values of the
#' amino acid properties explicitly, then use [grantham_equation()] instead.
#'
#' @param x A character vector of amino acid three-letter codes.
#' @param y A character vector of amino acid three-letter codes.
#' @param method Either `"original"` (default) or `"exact"`, see description for
#'   more details.
#' @param alpha The constant \eqn{\alpha} in the equation of Grantham's
#'   paper, in page 863.
#' @param beta The constant \eqn{\beta} in the equation of Grantham's
#'   paper, in page 863.
#' @param gamma The constant \eqn{\gamma} in the equation of Grantham's
#'   paper, in page 863.
#' @param rho Grantham's distances reported in Table 2, Science (1974).
#'   185(4154): 862--4 by R. Grantham, are scaled by a factor (here named
#'   \eqn{\rho}) such that the mean value of all distances are 100. The `rho`
#'   parameter allows this factor \eqn{\rho} to be changed. By default
#'   \eqn{\rho=50.723}, the same value used by Grantham. This value is
#'   originally mentioned in the caption of Table 2 of the aforementioned paper.
#'
#' @return A [tibble][tibble::tibble-package] of Grantham's distances for each
#'   amino acid pair.
#'
#' @md
#'
#' @source \doi{10.1126/science.185.4154.862}.
#'
#' @examples
#' # Grantham's distance between Serine (Ser) and Glutamate (Glu)
#' grantham_distance('Ser', 'Glu')
#'
#' # Grantham's distance between Serine (Ser) and Glutamate (Glu)
#' # with the "exact" method
#' grantham_distance('Ser', 'Glu', method = 'exact')
#'
#' # `grantham_distance()` is vectorised
#' # amino acids are paired element-wise between `x` and `y`
#' grantham_distance(x = c('Pro', 'Gly'), y = c('Glu', 'Arg'))
#'
#' # Use `amino_acid_pairs()` to generate pairs (by default generates all pairs)
#' aa_pairs <- amino_acid_pairs()
#' grantham_distance(x = aa_pairs$x, y = aa_pairs$y)
#'
#' @export
grantham_distance <-
  function(x,
           y,
           method = c('original', 'exact'),
           alpha = 1.833,
           beta = 0.1018,
           gamma = 0.000399,
           rho = 50.723) {

  if(!all_amino_acids(x))
    stop('`x` should contain only amino acid three-letter codes.')

  if(!all_amino_acids(y))
    stop('`y` should contain only amino acid three-letter codes.')

  # `rec`: recycled vectors `x` and `y`:
  rec <- vctrs::vec_recycle_common(x = x, y = y)

  # Check that `method` is either 'original' or 'exact'.
  method <- match.arg(method)

  if(identical(method, 'original'))
    return(grantham_distance_original(x = rec$x,
                                      y = rec$y))
  else
    return(
      grantham_distance_exact(
        x = rec$x,
        y = rec$y,
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        rho = rho
      )
    )
}

#' Grantham's distance (original)
#'
#' This function calculates the Grantham's distance for pairs of amino acids. It
#' uses the pre-calculated distances for each amino acid pair as published in
#' Table 2 of Science (1974). 185(4154): 862--4 by R. Grantham.
#'
#' @param x A character vector of amino acid three-letter codes.
#' @param y A character vector of amino acid three-letter codes.
#'
#' @return A [tibble][tibble::tibble-package] of Grantham's distances for each
#'   amino acid pair.
#'
#' @md
#' @source \doi{10.1126/science.185.4154.862}.
#' @keywords internal
#' @export
grantham_distance_original <- function(x, y) {

  amino_acid_pairs <- matrix(c(aa_idx(x), aa_idx(y)), ncol = 2)
  tbl <- tibble::tibble(x = x, y = y, d = grantham_distances_matrix[amino_acid_pairs])

  return(tbl)
}

#' Grantham's distance (exact)
#'
#' @md
#'
#' @description
#' This function calculates the Grantham's distance for pairs of amino acids. It
#' uses the values for the amino acid properties as published in Table 1 of
#' Science (1974). 185(4154): 862--4 by R. Grantham.
#'
#' @details
#' Contrary to Grantham's distances presented in Table 2 of Grantham's paper, the
#' distances returned by this function are calculated anew starting from the
#' amino acid properties (composition, polarity and molecular volume). No
#' rounding to nearest integer is performed.
#'
#' @param x A character vector of amino acid three-letter codes, e.g. `"Ala"`
#'   (Alanine).
#' @param y A character vector of amino acid three-letter codes.
#' @param alpha The constant \eqn{\alpha} in the equation of Grantham's
#'   paper, in page 863.
#' @param beta The constant \eqn{\beta} in the equation of Grantham's
#'   paper, in page 863.
#' @param gamma The constant \eqn{\gamma} in the equation of Grantham's
#'   paper, in page 863.
#' @param rho Grantham's distances reported in Table 2, Science (1974).
#'   185(4154): 862--4 by R. Grantham, are scaled by a factor (here named
#'   \eqn{\rho}) such that the mean value of all distances are 100. The `rho`
#'   parameter allows this factor \eqn{\rho} to be changed. By default
#'   \eqn{\rho=50.723}, the same value used by Grantham. This value is
#'   originally mentioned in the caption of Table 2 of the aforementioned paper.
#'
#' @return A [tibble][tibble::tibble-package] of Grantham's distances for each
#'   amino acid pair.
#' @source \doi{10.1126/science.185.4154.862}.
#'
#' @seealso [grantham_equation()]
#'
#' @examples
#' grantham_distance_exact(c('Ser', 'Ser'), c('Pro', 'Trp'))
#'
#' @keywords internal
#' @export
grantham_distance_exact <- function(x,
                                    y,
                                    alpha = 1.833,
                                    beta = 0.1018,
                                    gamma = 0.000399,
                                    rho = 50.723) {

  # Filter the properties table for the queried amino acids
  x_tbl <- amino_acids_properties[aa_idx(x), ]
  y_tbl <- amino_acids_properties[aa_idx(y), ]

  # Grantham's distance computed from the amino acids' properties as provided in
  # Table 1 of Grantham (1974).
  d <- grantham_equation(c_i = x_tbl$c,
                    c_j = y_tbl$c,
                    p_i = x_tbl$p,
                    p_j = y_tbl$p,
                    v_i = x_tbl$v,
                    v_j = y_tbl$v,
                    alpha = alpha,
                    beta = beta,
                    gamma = gamma,
                    rho = rho
                    )

  tbl <- tibble::tibble(x = x, y = y, d = d)

  return(tbl)
}
