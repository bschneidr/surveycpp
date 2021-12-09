//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat arma_onestage(arma::mat Y,
                        arma::colvec samp_unit_ids,
                        arma::colvec strata,
                        arma::colvec strata_pop_sizes,
                        Rcpp::CharacterVector singleton_method) {

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
                          arma::mat strata_pop_sizes,
                          Rcpp::CharacterVector singleton_method) {

  size_t n_stages = samp_unit_ids.n_cols;

  // If there are later stages of sampling,
  // obtain the necessary columns from inputs,
  // which will be used recursively

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
                              strata_pop_sizes = first_stage_strata_pop_sizes,
                              singleton_method = singleton_method);

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
        // Create subsets of inputs specific to current first-stage sampling unit
        arma::uvec unit_indices = arma::find(first_stage_ids.elem(h_indices) == h_unique_psus(i));
        arma::mat Y_hi = Y_h.rows(unit_indices);
        arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
        arma::mat hi_strata = h_strata.rows(unit_indices);
        arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);

        // Estimate later-stage variance contribution
        arma::mat V_hi = f_h * arma_multistage(Y_hi,
                                               hi_samp_unit_ids,
                                               hi_strata,
                                               hi_strata_pop_sizes,
                                               singleton_method);
        V += V_hi;

      }
    }
  }
  return V;
}

/*** R
# Create an example survey design ----

  set.seed(1999)

  library(survey)
  options("survey.ultimate.cluster" = FALSE)
  data(api, package = 'survey')

##_ Create a fake population to sample from
  population <- MASS::mvrnorm(n = 100000,
                              mu = colMeans(apipop[,c('api00', 'api99')]),
                              Sigma = cov(apipop[,c('api00', 'api99')])) |>
    apply(MARGIN = 2, FUN = round) |> `colnames<-`(c("api00", "api99")) |>
    as.data.frame()

  population <- cbind(population,
                      'stratum' = sample(x = c(1:10),
                                         size = nrow(population),
                                         replace = TRUE),
                      'psu_id' = sample(x = c(1:500),
                                        size = nrow(population),
                                        replace = TRUE),
                      'ssu_id' = sample(x = c(1:50),
                                        size = nrow(population),
                                        replace = TRUE))
  rownames(population) <- 1:nrow(population)
  population <- population[,c('stratum', 'psu_id', 'ssu_id',
                              'api00', 'api99')]

##_ Add columns giving population sizes needed for FPCs
  population <- aggregate(x = population$psu_id,
                          by = population[,'stratum', drop = FALSE],
                          FUN = function(x) length(unique(x))) |>
    setNames(c("stratum", 'N_psus')) |>
    merge(x = population,
          by = c("stratum"))

  population <- aggregate(x = population$ssu_id,
                          by = population[,c('stratum', 'psu_id'), drop = FALSE],
                          FUN = function(x) length(unique(x))) |>
    setNames(c("stratum", 'psu_id' , 'N_ssus')) |>
    merge(x = population,
          by = c("stratum", "psu_id"))

##_ Draw stratified multistage sample
  population$is_sampled <- FALSE

  for (stratum_h in unique(population$stratum)) {
    stratum_psus <- population |>
      subset(stratum == stratum_h) |>
      getElement("psu_id") |>
      unique()

    sampled_psus <- sample(stratum_psus, size = 100, replace = FALSE)

    for (h_psu in sampled_psus) {

      ssus_in_psu_of_stratum_h <- population |>
        subset(stratum == stratum_h & psu_id == h_psu) |>
        getElement("ssu_id") |>
        unique()

      sampled_ssus <- sample(ssus_in_psu_of_stratum_h, size = 5, replace = FALSE)

      sample_indices <- which(population$stratum == stratum_h &
                              population$psu_id == h_psu &
                              population$ssu_id %in% sampled_ssus)

      population[['is_sampled']][sample_indices] <- TRUE

    }
  }

  sample_data <- population[population$is_sampled,]

##_Declare survey design ----
  multistage_design <- svydesign(strata = ~ stratum,
                                 id = ~ psu_id + ssu_id,
                                 fpc = ~ N_psus + N_ssus,
                                 nest = TRUE,
                                 data = sample_data)

# Extract inputs ----
  Y = as.matrix(multistage_design$variables[,c('api00', 'api99')])
  Y_wtd <- Y |> apply(MARGIN = 2,
                      function(x) x/multistage_design$prob)
  ##_ Strata, clusters, and pop sizes
  ##_ need to be represented as numeric matrices
  strata = lapply(multistage_design$strata,
                  MARGIN = 2, FUN = as.numeric) |>
    do.call(what = cbind)
  clusters <- multistage_design$cluster %>%
    lapply(as.numeric) %>% Reduce(f = cbind)

  strata_pop_sizes = as.matrix(multistage_design$fpc$popsize)


  fpcs <- multistage_design$fpc

# Check that results match ----

  testthat::expect_equal(
    object = arma_multistage(Y = Y_wtd,
                             samp_unit_ids = clusters,
                             strata = strata,
                             strata_pop_sizes = strata_pop_sizes,
                             singleton_method = 'average'),
    expected =   survey:::multistage(x = Y_wtd,
                                     clusters = clusters,
                                     stratas = strata,
                                     nPSUs = fpcs$sampsize,
                                     fpcs = fpcs$popsize,
                                     cal = NULL) |> unname()
  )

# Benchmark ----

  microbenchmark::microbenchmark(
    'arma_multistage' = arma_multistage(Y = Y_wtd,
                                        samp_unit_ids = clusters,
                                        strata = strata,
                                        strata_pop_sizes = strata_pop_sizes,
                                        singleton_method = 'average'),
    'survey:::multistage' =   survey:::multistage(x = Y_wtd,
                                                  clusters = clusters,
                                                  stratas = strata,
                                                  nPSUs = fpcs$sampsize,
                                                  fpcs = fpcs$popsize,
                                                  cal = NULL)
  )
*/
