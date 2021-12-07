Rcpp::cppFunction(code = "
    IntegerVector recursive_count (arma::mat Y) {

    // Initialize count
    int n = 0;

    // Determine number of columns of input
    size_t n_cols_y = Y.n_cols;

    // Count number of distinct values in first column
    arma::colvec first_col = Y.col(0);
    arma::colvec unique_first_col = unique(first_col);
    n += unique_first_col.n_elem;

    // If there are multiple columns,
    // recursively call the function
    // for each distinct value in first column

    if (n_cols_y > 1) {
      for (arma::uword i=0; i < unique_first_col.n_elem; ++i ) {
         // Make a copy of Y with all columns except first,
         // and all the rows corresponding to the i-th unique value from first column
         arma::uvec i_row_indices = arma::find(first_col == unique_first_col(i));
         arma::mat Y_i = Y.rows(i_row_indices);
         Y_i = Y_i.tail_cols(n_cols_y - 1);

         // Call the function on the submatrix
         // and add the result
         IntegerVector n_from_submatrix(1);
         n_from_submatrix = recursive_count(Y_i);
         n += n_from_submatrix(0);
      }
    }
    IntegerVector return_value = wrap(n);
    return return_value;
    }
", depends = "RcppArmadillo")

# Create an example matrix
example_matrix <- matrix(
  c(1, 1, 1,
    1, 1, 2,
    1, 2, 1,
    2, 1, 1,
    2, 2, 1),
  byrow = TRUE,
  ncol = 3)

# Use the Rcpp function
# expected output is
recursive_count(example_matrix)
