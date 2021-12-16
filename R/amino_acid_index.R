amino_acid_index <- function(amino_acid) {
  match(amino_acid, amino_acids())
}

# aa_idx: abbreviated form of `amino_acid_index`.
aa_idx <- amino_acid_index
