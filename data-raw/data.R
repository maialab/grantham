library(readr)
library(here)
library(grantham)

# Grantham distances' matrix
grantham_distances_matrix <-
  readr::read_csv(
    file = here::here('data-raw', 'grantham_distance_matrix.csv'),
    col_types = 'ciiiiiiiiiiiiiiiiiii',
    col_select = -1
  ) %>%
  as.matrix() %>%
  `rownames<-`(., colnames(.))

# Sort the rows and columns by the order present in `amino_acids()`. This
# ordering should already be as in the return value of `amino_acids()`, but just
# in case...
grantham_distances_matrix <-
  grantham_distances_matrix[amino_acids(), amino_acids()]

# The values for the amino acid properties in "amino_acids_properties.csv" were
# directly obtained from Table 1 of Grantham (1974).
amino_acids_properties <-
  readr::read_csv(
    file = here::here('data-raw', 'amino_acids_properties.csv'),
    col_types = 'cdd'
  ) %>% # Next line is just ensure that the order comes out the same as in `amino_acids()`.
  dplyr::left_join(tibble::tibble(amino_acid = amino_acids()), ., by = 'amino_acid')

# The 20 amino acids.
n_amino_acids <- length(amino_acids())

mean_chemical_distance <-
  with(amino_acids_properties,
       c(
         'c' = mean(outer(c, c, function(x, y) abs(x - y))[grantham:::sltm_k(n_amino_acids)]),
         'p' = mean(outer(p, p, function(x, y) abs(x - y))[grantham:::sltm_k(n_amino_acids)]),
         'v' = mean(outer(v, v, function(x, y) abs(x - y))[grantham:::sltm_k(n_amino_acids)])
       )
  ) %>%
  signif(digits = 4) %>%
  round(digits = 3)

# The mean weighting factors (as they are referred to in the caption of Table 1
# of R. Grantham (1974) as used as indicated in that caption. If we were to
# calculate them here from the `mean_chemical_distance` one would find that the
# alpha value (1.833) is slightly off by a small percentage 0.11% (calculated
# value is 1.831.)
# As the difference is relatively minor, we stick with the values reported in
# the original paper to avoid confusion.
mean_weighting_factors <- c('alpha' = 1.833, 'beta' = 0.1018, 'gamma' = 0.000399)

one_letter_codes <- c(
  'A', # Alanine
  'R', # Arginine
  'N', # Asparagine
  'D', # Aspartic acid
  'B', # Asparagine or aspartic acid
  'C', # Cysteine
  'E', # Glutamic acid
  'Q', # Glutamine
  'Z', # Glutamine or glutamic acid
  'G', # Glycine
  'H', # Histidine
  'I', # Isoleucine
  'L', # Leucine
  'K', # Lysine
  'M', # Methionine
  'F', # Phenylalanine
  'P', # Proline
  'S', # Serine
  'T', # Threonine
  'W', # Tryptophan
  'Y', # Tyrosine
  'V'  # Valine
)

three_letter_codes <- c(
  'Ala', # Alanine
  'Arg', # Arginine
  'Asn', # Asparagine
  'Asp', # Aspartic acid
  'Asx', # Asparagine or aspartic acid
  'Cys', # Cysteine
  'Glu', # Glutamic acid
  'Gln', # Glutamine
  'Glx', # Glutamine or glutamic acid
  'Gly', # Glycine
  'His', # Histidine
  'Ile', # Isoleucine
  'Leu', # Leucine
  'Lys', # Lysine
  'Met', # Methionine
  'Phe', # Phenylalanine
  'Pro', # Proline
  'Ser', # Serine
  'Thr', # Threonine
  'Trp', # Tryptophan
  'Tyr', # Tyrosine
  'Val'  # Valine
)

# These variables end up in R/sysdata.rda
usethis::use_data(
  amino_acids_properties,
  grantham_distances_matrix,
  mean_chemical_distance,
  mean_weighting_factors,
  one_letter_codes,
  three_letter_codes,
  internal = TRUE,
  overwrite = TRUE
)

# These end up in data/*.rda
usethis::use_data(amino_acids_properties, overwrite = TRUE)
usethis::use_data(grantham_distances_matrix, overwrite = TRUE)

