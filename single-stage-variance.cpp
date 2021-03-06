//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <unordered_set>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//------------------------------------------------------------------------------
//
//
// The function `group_sums()` calculates column totals by group
// The functions `submat()` and `rowvec()` are helpers used inside `group_sums()`
//
//------------------------------------------------------------------------------

//[[Rcpp::export]]
arma::mat submat(NumericMatrix X, NumericVector T, int TestVal) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false);
  arma::mat y = Xmat.rows(find(tIdx == TestVal));
  return y;
}

// [[Rcpp::export]]
arma::mat remove_col(arma::mat x, int which){
  x.shed_col(which);
  return(x);
}

// [[Rcpp::export]]
arma::rowvec col_sums(arma::mat x){
  arma::mat X = arma::mat(x.begin(), x.n_rows, x.n_cols, false);
  return arma::sum(X, 0);
}

//[[Rcpp::export]]
NumericMatrix rcpp_group_sums(NumericMatrix X, NumericVector by_var) {
  NumericVector by_var_levels = sort_unique(by_var);
  int n = by_var_levels.size();
  int m = X.ncol();
  arma::mat out(n, m);
  for (int i(0); i < n; i++) {
    int level = by_var_levels(i);
    arma::mat sub = submat(X, by_var, level);
    arma::rowvec colsums = col_sums(sub);
    out.row(i) = colsums;
  }

  NumericMatrix named_out = wrap(out);
  rownames(named_out) = by_var_levels;
  colnames(named_out) = colnames(X);
  return named_out;
}

//[[Rcpp::export]]
arma::mat arma_group_sums(NumericMatrix X, NumericVector by_var) {
  NumericVector by_var_levels = sort_unique(by_var);
  int n = by_var_levels.size();
  int m = X.ncol();
  arma::mat out(n, m);
  for (int i(0); i < n; i++) {
    int level = by_var_levels(i);
    arma::mat sub = submat(X, by_var, level);
    arma::rowvec colsums = col_sums(sub);
    out.row(i) = colsums;
  }

  return out;
}
//------------------------------------------------------------------------------
//
//
// Below are helper functions used for obtaining single stage variance
//
//
//------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix covRcpp(Rcpp::NumericMatrix & X,
                      const int norm_type = 0) {

  const int n = X.nrow();
  const int m = X.ncol();
  const int df = n - 1 + norm_type;

  Rcpp::NumericMatrix cov(m, m);

  // Degenerate case
  if (n == 0) {
    std::fill(cov.begin(), cov.end(), Rcpp::NumericVector::get_na());
    return cov;
  }

  // Centering the matrix!
  for (int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }

  // Computing the covariance matrix
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      cov(i,j) = Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
  }

  return cov;
}

// [[Rcpp::export]]
arma::mat arma_cov(arma::mat Y,
                   arma::uvec ids) {
  arma::mat cov_mat(Y.n_cols, Y.n_cols, arma::fill::zeros);

  arma::uvec unique_ids = unique(ids);

  arma::uword n = unique_ids.n_elem;

  arma::uword df;
  if (n > 1) {
    df = n - 1;
  } else {
    df = 1;
  }

  arma::rowvec mean_Yi = sum(Y, 0) / n;
  for (arma::uword i=0; i < n; ++i) {
    arma::uvec s = arma::find(ids == unique_ids[i]);
    arma::rowvec Yi = sum(Y.rows(s), 0);
    Yi = Yi - mean_Yi;
    cov_mat += (arma::trans(Yi)*Yi);
  }
  cov_mat = cov_mat / df;
  return cov_mat;
}

// [[Rcpp::export]]
NumericMatrix cpp_subset_rows(NumericMatrix m, IntegerVector rows){

  int rl = rows.length();
  int cl = m.ncol();
  NumericMatrix out(rl, cl);

  for (int j=0; j<cl; j++){
    NumericMatrix::Column org_c = m(_, j);
    NumericMatrix::Column new_c = out(_, j);
    for (int i=0; i<rl; i++){
      new_c[i] = org_c[rows[i]-1];
    }
  }
  return(out);
}

// [[Rcpp::export]]
IntegerVector get_matching_indices(CharacterVector x, CharacterVector target) {

  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);

  for(int i = 0; i < nx; i++) {
    if (x[i] == target[0]) y.push_back(i+1);
  }

  return wrap(y);
}


//------------------------------------------------------------------------------
//
//
// The following function estimates single stage variance
// of a total assuming stratified, single-stage SRSWOR.
//
//
//------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix single_stage_variance(Rcpp::NumericMatrix Y,
                                          Rcpp::NumericVector samp_unit_ids,
                                          Rcpp::CharacterVector strata,
                                          Rcpp::IntegerVector strata_pop_sizes,
                                          Rcpp::CharacterVector singleton_method) {

  //char singleton_method = singleton_method;

  // Determine dimensions of result
  const int n_col_y = Y.ncol();
  Rcpp::NumericMatrix result(n_col_y, n_col_y);
  colnames(result) = colnames(Y);
  rownames(result) = colnames(Y);

  // Get distinct strata ids and their length, H
  Rcpp::CharacterVector distinct_strata_ids = Rcpp::unique(strata);

  const int H = distinct_strata_ids.length();

  // Get number of distinct sampling units
  Rcpp::NumericVector distinct_samp_unit_ids= Rcpp::unique(samp_unit_ids);
  const int n_samp_units = distinct_samp_unit_ids.length();

  // Recenter Y around grand mean of sampling units
  // if singleton_method is set accordingly
  if (singleton_method[0] == "adjust") {
    for (int j = 0; j < n_col_y; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - (Rcpp::sum(Y(Rcpp::_, j))/n_samp_units);
    }
  }

  // Initialize vectors giving number of PSUs per stratum
  // and vector giving sampling fraction
  Rcpp::IntegerVector n_h (H, 1);
  Rcpp::IntegerVector N_h (H, 1);
  Rcpp::NumericVector f_h (H, 0.0);

  //Initialize count of singleton strata
  int n_singleton_strata = 0;

  // Get information from each stratum:
  // - Number of sampling units, in sample and population
  // - Sampling fraction, if applicable
  // - Contribution to sampling variance
  for (int h = 0; h < H; ++h) {
    // Determine which rows of data correspond to the current stratum
    Rcpp::CharacterVector current_stratum = as<CharacterVector>((distinct_strata_ids[h]));
    Rcpp::IntegerVector h_indices = get_matching_indices(strata, current_stratum);

    // Get count of sampling units in stratum
    Rcpp::NumericVector h_samp_unit_ids = samp_unit_ids[h_indices - 1];
    Rcpp::NumericVector h_distinct_samp_unit_ids = Rcpp::unique(h_samp_unit_ids);
    n_h[h] = h_distinct_samp_unit_ids.length();

    // Get finite population correction
    N_h[h] = as<IntegerVector>(strata_pop_sizes[h_indices - 1])[0];
    f_h[h] = static_cast<double>(n_h[h]) /  static_cast<double>(N_h[h]);

    // Get weighted sampling unit totals in stratum
    Rcpp::NumericMatrix h_observations = cpp_subset_rows(Y, h_indices);
    arma::mat samp_unit_totals = arma_group_sums(h_observations, h_samp_unit_ids);

    // Estimate variance of weighted sampling unit totals
    Rcpp::NumericMatrix cov_mat(n_col_y,n_col_y);
    if (n_h[h] > 1) {
      cov_mat = wrap(arma::cov(samp_unit_totals));
    } else {
      n_singleton_strata += 1;
      if (singleton_method[0] == "adjust") {
        cov_mat = wrap(samp_unit_totals.t() * samp_unit_totals);
      }
    }

    // Add variance contribution
    result += (1 - f_h[h]) * n_h[h] * cov_mat;
  }

  if ((n_singleton_strata > 0) & (singleton_method[0] == "average")) {
    result += n_singleton_strata * result;
  }

  return result;
}

// [[Rcpp::export]]
arma::mat arma_onestage(arma::mat Y,
                        arma::colvec samp_unit_ids,
                        arma::colvec strata,
                        arma::colvec strata_pop_sizes) {

  // Determine dimensions of result
  size_t n_col_y = Y.n_cols;
  arma::mat result(n_col_y, n_col_y, arma::fill::zeros);

  // Get distinct strata ids and their length, H
  arma::colvec distinct_strata_ids = unique(strata);

  arma::uword H = distinct_strata_ids.n_elem;

  // Initialize vectors giving number of PSUs per stratum
  // and vector giving sampling fraction
  arma::uvec n_h (H, arma::fill::ones);
  arma::uvec N_h (H, arma::fill::ones);
  arma::dvec f_h (H, arma::fill::zeros);

  // Get information from each stratum:
  // - Number of sampling units, in sample and population
  // - Sampling fraction, if applicable
  // - Contribution to sampling variance
  for (arma::uword h = 0; h < H; ++h) {

    // Determine which rows of data correspond to the current stratum
    arma::uvec h_indices = arma::find(strata==distinct_strata_ids[h]);

    // Get count of sampling units in stratum
    arma::colvec h_distinct_samp_unit_ids = unique(samp_unit_ids.elem(h_indices));
    n_h[h] = h_distinct_samp_unit_ids.n_elem;
    N_h[h] = min(strata_pop_sizes.elem(h_indices));
    f_h[h] = static_cast<double>(n_h[h]) /  static_cast<double>(N_h[h]);

    arma::mat Y_h = Y.rows(h_indices);
    arma::rowvec mean_Yhi = arma::sum(Y_h, 0) / n_h[h];

    arma::uword df;
    if (n_h[h] > 1) {
      df = n_h[h] - 1;
    } else {
      df = 1;
    }
    arma::mat cov_mat(n_col_y, n_col_y, arma::fill::zeros);
    for (arma::uword i=0; i < h_distinct_samp_unit_ids.n_elem; ++i ) {
      arma::uvec unit_indices = arma::find(samp_unit_ids.elem(h_indices) == h_distinct_samp_unit_ids[i]);
      arma::rowvec Yhi = sum(Y_h.rows(unit_indices), 0);
      Yhi = Yhi - mean_Yhi;

      cov_mat += (arma::trans(Yhi)*Yhi);
    }
    cov_mat = cov_mat / df;

    result += ((1 - f_h[h]) * n_h[h]) * cov_mat;
  }

  return result;
}

// [[Rcpp::export]]
arma::mat arma_multistage(arma::mat Y,
                          arma::mat samp_unit_ids,
                          arma::mat strata,
                          arma::mat strata_pop_sizes) {

  size_t n_stages = samp_unit_ids.n_cols;

  arma::mat later_stage_ids;
  arma::mat later_stage_strata;
  arma::mat later_stage_strata_pop_sizes;

  if (n_stages > 1) {
    later_stage_ids = samp_unit_ids.tail_cols(n_stages - 1);
    later_stage_strata = strata.tail_cols(n_stages - 1);
    later_stage_strata_pop_sizes = strata_pop_sizes.tail_cols(n_stages-1);
  }

  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);

  // Calculate first-stage variance
  arma::mat V = arma_onestage(Y = Y,
                              samp_unit_ids = first_stage_ids,
                              strata = first_stage_strata,
                              strata_pop_sizes = first_stage_strata_pop_sizes);

  // For each first-stage unit, get variance contribution from next stage
  if (n_stages > 1) {

    // Get distinct first-stage strata ids and their length, H
    arma::colvec distinct_strata_ids = unique(first_stage_strata);
    arma::uword H = distinct_strata_ids.n_elem;

    for (arma::uword h = 0; h < H; ++h) {

      // Determine which rows of data correspond to the current first-stage stratum
      arma::uvec h_indices = arma::find(first_stage_strata==distinct_strata_ids(h));

      // Get submatrices of inputs corresponding to the current first-stage stratum
      arma::mat Y_h = Y.rows(h_indices);

      arma::mat h_samp_unit_ids = later_stage_ids.rows(h_indices);

      arma::mat h_strata = later_stage_strata.rows(h_indices);

      arma::mat h_strata_pop_sizes = later_stage_strata_pop_sizes.rows(h_indices);

      // Get count of first-stage sampling units in first-stage stratum
      // and finite population correction
      arma::colvec h_first_stage_units = first_stage_ids.elem(h_indices);
      arma::colvec h_unique_psus = unique(h_first_stage_units);
      arma::uword n_h = h_unique_psus.n_elem;
      arma::uword N_h = min(first_stage_strata_pop_sizes.elem(h_indices));
      double f_h = static_cast<double>(n_h) /  static_cast<double>(N_h);

      for (arma::uword i=0; i < n_h; ++i ) {
        arma::uvec unit_indices = arma::find(first_stage_ids.elem(h_indices) == h_unique_psus(i));
        arma::mat Y_hi = Y_h.rows(unit_indices);
        arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
        //hi_samp_unit_ids = remove_col(h_samp_unit_ids, 0);
        arma::mat hi_strata = h_strata.rows(unit_indices);
        arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);

        arma::mat V_hi = f_h * arma_multistage(Y_hi,
                                               hi_samp_unit_ids,
                                               hi_strata,
                                               hi_strata_pop_sizes);
        V += V_hi;

      }
    }
  }
  return V;
}

// [[Rcpp::export]]
arma::mat arma_multistage2(arma::mat Y,
                          arma::mat samp_unit_ids,
                          arma::mat strata,
                          arma::mat strata_pop_sizes) {

  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);

  // Calculate first-stage variance
  arma::mat V = arma_onestage(Y = Y,
                              samp_unit_ids = first_stage_ids,
                              strata = first_stage_strata,
                              strata_pop_sizes = first_stage_strata_pop_sizes);

  return V;
}

// [[Rcpp::export]]
Rcpp::List cLists(List x, List y) {
  int nsize = x.size();
  int msize = y.size();
  Rcpp::List out(nsize + msize);

  for(int i = 0; i < nsize; i++) {
    out[i] = x[i];
  }
  for(int i = 0; i < msize; i++) {
    out[nsize+i] = y[i];
  }

  return(out);
}

// [[Rcpp::export]]
Rcpp::List append_list_to_list(List x, List y) {
  int nsize = x.size();
  Rcpp::List out(nsize + 1);
  for(int i = 0; i < nsize; i++) {
    out[i] = x[i];
  }
  out[nsize] = y;
  return out;
}


// [[Rcpp::export]]
Rcpp::List arma_multistage_debug (arma::mat Y,
                                  arma::mat samp_unit_ids,
                                  arma::mat strata,
                                  arma::mat strata_pop_sizes) {

  size_t n_stages = samp_unit_ids.n_cols;

  arma::mat later_stage_ids = samp_unit_ids.tail_cols(n_stages - 1);
  arma::mat later_stage_strata = strata.tail_cols(n_stages - 1);
  arma::mat later_stage_strata_pop_sizes = strata_pop_sizes.tail_cols(n_stages-1);

  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);

  // Calculate first-stage variance
  arma::mat V = arma_onestage(Y = Y,
                              samp_unit_ids = first_stage_ids,
                              strata = first_stage_strata,
                              strata_pop_sizes = first_stage_strata_pop_sizes);

  Rcpp::List result = Rcpp::List::create(n_stages);
  result[0] = wrap(samp_unit_ids);

  //arma::mat result;
  // For each first-stage unit, get variance contribution from next stage
  if (n_stages > 1) {

    // Get distinct first-stage strata ids and their length, H
    arma::colvec distinct_strata_ids = unique(first_stage_strata);
    arma::uword H = distinct_strata_ids.n_elem;

    for (arma::uword h = 0; h < H; ++h) {

      // Determine which rows of data correspond to the current first-stage stratum
      arma::uvec h_indices = arma::find(first_stage_strata==distinct_strata_ids(h));

      // Get submatrices of inputs corresponding to the current first-stage stratum
      arma::mat Y_h = Y.rows(h_indices);

      //result = Y_h;
      arma::mat h_samp_unit_ids = later_stage_ids.rows(h_indices);

      arma::mat h_strata = later_stage_strata.rows(h_indices);

      arma::mat h_strata_pop_sizes = later_stage_strata_pop_sizes.rows(h_indices);

      // Get count of first-stage sampling units in first-stage stratum
      // and finite population correction
      arma::colvec h_first_stage_units = first_stage_ids.elem(h_indices);
      arma::colvec h_unique_psus = unique(h_first_stage_units);
      arma::uword n_h = h_unique_psus.n_elem;
      arma::uword N_h = min(first_stage_strata_pop_sizes.elem(h_indices));
      double f_h = static_cast<double>(n_h) /  static_cast<double>(N_h);
      if (FALSE) {
        f_h = f_h * 2;
      }

      for (arma::uword i=0; i < n_h; ++i ) {
        arma::uvec unit_indices = arma::find(first_stage_ids.elem(h_indices) == h_unique_psus(i));
        arma::mat Y_hi = Y_h.rows(unit_indices);
        arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
        //hi_samp_unit_ids = remove_col(h_samp_unit_ids, 0);
        arma::mat hi_strata = h_strata.rows(unit_indices);
        arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);

        arma::mat V_hi = arma_multistage(Y_hi,
                                         hi_samp_unit_ids,
                                         hi_strata,
                                         hi_strata_pop_sizes);

        Rcpp::List result_hi = Rcpp::List::create(Rcpp::Named("Y_hi") = Y_hi,
                                                  Rcpp::Named("hi_samp_unit_ids") = hi_samp_unit_ids,
                                                  Rcpp::Named("hi_strata") = hi_strata,
                                                  Rcpp::Named("hi_strata_pop_sizes") = hi_strata_pop_sizes,
                                                  Rcpp::Named("V_hi") = V_hi);



        result = append_list_to_list(result, result_hi);

      }
      // // Get matrix of data for stratum
      //
      // for (arma::uword i=0; i < n_h; ++i ) {
      //   arma::uvec unit_indices = arma::find(h_first_stage_units == h_distinct_samp_unit_ids[i]);
      //
      //   // Get submatrices of inputs corresponding to the current first-stage sampling-unit
      //   arma::mat Y_hi = Y_h.rows(unit_indices);
      //   arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
      //   hi_samp_unit_ids = hi_samp_unit_ids.tail_cols(n_stages - 1);
      //   arma::mat hi_strata = h_strata.rows(unit_indices);
      //   hi_strata = hi_strata.tail_cols(n_stages - 1);
      //   arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);
      //   hi_strata_pop_sizes = hi_strata_pop_sizes.tail_cols(n_stages - 1);
      //
      //   // Get variance from later stages of sampling within current first-stage sampling unit
      //   if (TRUE) {
      //     V += f_h * arma_multistage(Y = Y_hi,
      //                                samp_unit_ids = hi_samp_unit_ids,
      //                                strata = hi_strata,
      //                                strata_pop_sizes = hi_strata_pop_sizes);
      //   }
      // }
    }
  }
  return result;
}


/*** R
if (FALSE) {
# Generate example data,
# starting from a small data frame
# and replicating its rows
df <- data.frame(
  Stratum = c(1, 1, 2, 2, 2, 2, 2),
  Stratum_Pop_Size = c(100, 100, 300, 300, 300, 300, 300),
  PSU = c(1, 2, 3, 4, 5, 5, 5),
  Design_Weight = c(10.5, 10.5, 20, 20, 40, 40, 40),
  Sex = c("M", "F", "F", "M", "M", "F", "M"),
  Age = c(30.5, 40.5, 35, 52, 60, 48, 33),
  Height = c(6.2, 5.0, 5.3, 5.7, 6.1, 5.2, 5.8)
)

large_df <- df
increment <- 1
for (i in 1:2000) {
  new_df <- df
  new_df[['Stratum']] <- new_df[['Stratum']] + (2*i)
  large_df <- rbind(large_df, new_df)
}

# Extract necessary objects to supply
# to R functions for variance estimation
Y = as.matrix(large_df[,c("Age", "Height")])
samp_unit_ids = large_df[['PSU']]
clusters = matrix(large_df$PSU, ncol = 1)
strata = large_df[['Stratum']]
strata_pop_sizes = large_df[['Stratum_Pop_Size']]
stratas = matrix(large_df[['Stratum']], ncol = 1)
fpcs = survey::svydesign(ids = ~ PSU, strata = ~ Stratum,
                         nest = TRUE, fpc = ~ Stratum_Pop_Size,
                         data = large_df)$fpc

single_stage_variance(Y = Y,
                      samp_unit_ids = samp_unit_ids,
                      strata = strata,
                      strata_pop_sizes = strata_pop_sizes,
                      singleton_method = 'remove')

arma_onestage(
  Y = Y,
  samp_unit_ids = samp_unit_ids,
  strata = strata,
  strata_pop_sizes = strata_pop_sizes
)

arma_multistage(
  Y = Y,
  samp_unit_ids = as.matrix(samp_unit_ids),
  strata = as.matrix(strata),
  strata_pop_sizes = as.matrix(strata_pop_sizes)
)

if (TRUE) {

  # Check that results match
  testthat::expect_equal(
    object = single_stage_variance(
      Y = Y,
      samp_unit_ids = samp_unit_ids,
      strata = strata,
      strata_pop_sizes = strata_pop_sizes,
      singleton_method = 'fail') |> unname(),
    expected = survey:::onestage(x = Y, strata = large_df[['Stratum']],
                                 clusters = large_df[['PSU']],
                                 nPSU = fpcs$sampsize,
                                 fpc = fpcs$popsize,
                                 stage = 1) |> unname()
  )

  # Check that results match
  testthat::expect_equal(
    object = arma_onestage(
      Y = Y,
      samp_unit_ids = samp_unit_ids,
      strata = strata,
      strata_pop_sizes = strata_pop_sizes) |> unname(),
    expected = survey:::onestage(x = Y, strata = large_df[['Stratum']],
                                 clusters = large_df[['PSU']],
                                 nPSU = fpcs$sampsize,
                                 fpc = fpcs$popsize,
                                 stage = 1) |> unname()
  )

  # Compare speed
  microbenchmark::microbenchmark(
    'Rcpp' = single_stage_variance(Y = Y,
                                   samp_unit_ids = samp_unit_ids,
                                   strata = strata,
                                   strata_pop_sizes = strata_pop_sizes,
                                   singleton_method = 'fail'),
    'arma' = arma_onestage(Y = Y,
                           samp_unit_ids = samp_unit_ids,
                           strata = strata,
                           strata_pop_sizes = strata_pop_sizes),
    'survey:::onestage()' = survey:::onestage(x = Y, strata = large_df[['Stratum']],
                                              clusters = large_df[['PSU']],
                                              nPSU = fpcs$sampsize,
                                              fpc = fpcs$popsize,
                                              stage = 1),
    'survey::svyrecvar()' = survey::svyrecvar(x = Y,
                                              clusters = clusters,
                                              stratas = stratas,
                                              fpcs = fpcs)
  )
}
}
*/
