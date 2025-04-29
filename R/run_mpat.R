#' Run MixFisher Pleiotropy Test (Davies Approximation)
#'
#' This function applies the `MixFisher` method from the `MPAT` package using the Davies approximation to assess pleiotropy significance.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `MixFisher` method is then applied using the Davies approximation to calculate p-values.
#'
#' @return A `data.frame` with `MixFisher`-based p-values using the Davies method for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_mixfisher_davies(pleio)
#' head(result)
#'
#' @export
run_mixfisher_davies = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixFisher, Sigma = zscore_sigma,  method = "davies")[1,]

  pleio_p = data.frame(mixfisher_davies_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixFisher Pleiotropy Test (Liu Approximation)
#'
#' Applies the `MixFisher` method using the Liu approximation to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' Similar to `run_mixfisher_davies()`, but uses the Liu approximation instead of Davies.
#'
#' @return A `data.frame` with `MixFisher`-based p-values using the Liu method for each variant.
#'
#' @examples
#' result <- run_mixfisher_liu(pleio)
#' head(result)
#'
#' @export
run_mixfisher_liu = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixFisher, Sigma = zscore_sigma,  method = "liu")[1,]

  pleio_p = data.frame(mixfisher_liu_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixFisher Pleiotropy Test (Liu-modified Approximation)
#'
#' Applies the `MixFisher` method using the Liu-modified approximation to test for pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' Uses a modified Liu approximation for estimating p-values in the `MixFisher` method.
#'
#' @return A `data.frame` with `MixFisher`-based p-values using the Liu-modified method for each variant.
#'
#' @examples
#' result <- run_mixfisher_liumod(pleio)
#' head(result)
#'
#' @export
run_mixfisher_liumod = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixFisher, Sigma = zscore_sigma,  method = "liumod")[1,]

  pleio_p = data.frame(mixfisher_liumod_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}


#' Run MixSD Pleiotropy Test (Davies Approximation)
#'
#' Implements the `MixSD` test for pleiotropy using the Davies approximation.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function applies the `MixSD` method to compute p-values.
#'
#' @return A `data.frame` with `MixSD`-based p-values using the Davies method for each variant.
#'
#' @examples
#' result <- run_mixsd_davies(pleio)
#' head(result)
#'
#' @export
run_mixsd_davies = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixSD, Sigma = zscore_sigma,  method = "davies")

  pleio_p = data.frame(mixsd_davies_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixTippett Pleiotropy Test (Liu Approximation)
#'
#' Applies the `MixTippett` method using the Liu approximation to test for pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' This function applies the `MixTippett` method using the Liu approximation.
#'
#' @return A `data.frame` with `MixTippett`-based p-values using the Liu method for each variant.
#'
#' @examples
#' result <- run_mixtippett_liu(pleio)
#' head(result)
#'
#' @export
run_mixsd_liu = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixSD, Sigma = zscore_sigma,  method = "liu")

  pleio_p = data.frame(mixsd_liu_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixSD Pleiotropy Test (Liu-modified Approximation)
#'
#' Applies the `MixSD` method using the Liu-modified approximation to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' This function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `MixSD` method is then applied using the Liu-modified approximation to calculate p-values.
#'
#' @return A `data.frame` with `MixSD`-based p-values using the Liu-modified method for each variant.
#'
#' @examples
#' result <- run_mixsd_liumod(pleio)
#' head(result)
#'
#' @export
run_mixsd_liumod = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixSD, Sigma = zscore_sigma,  method = "liumod")

  pleio_p = data.frame(mixsd_liumod_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixTippett Pleiotropy Test (Davies Approximation)
#'
#' Applies the `MixTippett` method from the `MPAT` package using the Davies approximation to evaluate pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `MixTippett` method is then applied using the Davies approximation to compute p-values.
#'
#' @return A `data.frame` with `MixTippett`-based p-values using the Davies method for each variant.
#'
#' @examples
#' result <- run_mixtippett_davies(pleio)
#' head(result)
#'
#' @export
run_mixtippett_davies = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixTippett, Sigma = zscore_sigma,  method = "davies")[1,]

  pleio_p = data.frame(mixtippett_davies_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixTippett Pleiotropy Test (Liu Approximation)
#'
#' Applies the `MixTippett` method using the Liu approximation to test for pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function applies the `MixTippett` method using the Liu approximation.
#'
#' @return A `data.frame` with `MixTippett`-based p-values using the Liu method for each variant.
#'
#' @examples
#' result <- run_mixtippett_liu(pleio)
#' head(result)
#'
#' @export
run_mixtippett_liu = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixTippett, Sigma = zscore_sigma,  method = "liu")[1,]

  pleio_p = data.frame(mixtippett_liu_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixTippett Pleiotropy Test (Liu-modified Approximation)
#'
#' Applies the `MixTippett` method from the `MPAT` package using the Liu-modified approximation to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `MixTippett` method is a statistical approach for combining p-values across multiple correlated traits.
#' This function applies the Liu-modified approximation to compute pleiotropy significance while accounting for trait correlations.
#'
#' @return A `data.frame` with `MixTippett`-based p-values using the Liu-modified method for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_mixtippett_liumod(pleio)
#' head(result)
#'
#' @export
run_mixtippett_liumod = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixTippett, Sigma = zscore_sigma,  method = "liumod")[1,]

  pleio_p = data.frame(mixtippett_liumod_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixVar Pleiotropy Test (Davies Approximation)
#'
#' Implements the `MixVar` test for pleiotropy using the Davies approximation.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function applies the `MixVar` method to compute p-values using the Davies approximation.
#'
#' @return A `data.frame` with `MixVar`-based p-values using the Davies method for each variant.
#'
#' @examples
#' result <- run_mixvar_davies(pleio)
#' head(result)
#'
#' @export
run_mixvar_davies = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixVar, Sigma = zscore_sigma,  method = "davies")

  pleio_p = data.frame(mixvar_davies_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixVar Pleiotropy Test (Liu Approximation)
#'
#' Applies the `MixVar` method from the `MPAT` package using the Liu approximation to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `MixVar` method aggregates association signals across multiple correlated traits by modeling their variance.
#' This function applies the Liu approximation to estimate the statistical significance of pleiotropy.
#'
#' @return A `data.frame` with `MixVar`-based p-values using the Liu method for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_mixvar_liu(pleio)
#' head(result)
#'
#' @export
run_mixvar_liu = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixVar, Sigma = zscore_sigma,  method = "liu")

  pleio_p = data.frame(mixvar_liu_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixVar Pleiotropy Test (Liu-modified Approximation)
#'
#' Applies the `MixVar` method from the `MPAT` package using the Liu-modified approximation to test for pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `MixVar` method combines variance-based association signals across multiple correlated traits.
#' This function uses a modified Liu approximation to compute pleiotropy significance.
#'
#' @return A `data.frame` with `MixVar`-based p-values using the Liu-modified method for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_mixvar_liumod(pleio)
#' head(result)
#'
#' @export
run_mixvar_liumod = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixVar, Sigma = zscore_sigma,  method = "liumod")

  pleio_p = data.frame(mixvar_liumod_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run DSUM Pleiotropy Test
#'
#' This function applies the `DSUM` method from the `MPAT` package to evaluate pleiotropy significance.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `DSUM` method is then applied to compute p-values.
#'
#' @return A `data.frame` with `DSUM`-based p-values for each variant.
#'
#' @examples
#' result <- run_dsum(pleio)
#' head(result)
#'
#' @export
run_dsum = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::DSUM, Sigma = zscore_sigma)

  pleio_p = data.frame(dsum_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MinP Pleiotropy Test
#'
#' Applies the `MinP` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' This function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `MinP` method is then applied to compute the minimum p-value across phenotypes.
#'
#' @return A `data.frame` with `MinP`-based p-values for each variant.
#'
#' @examples
#' result <- run_minp(pleio)
#' head(result)
#'
#' @export
run_minp = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::MinP, Sigma = zscore_sigma)

  pleio_p = data.frame(minp_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run MixAda Pleiotropy Test
#'
#' Applies the `mixAda` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' Uses adaptive combination of association statistics across multiple traits.
#'
#' @return A `data.frame` with `mixAda`-based p-values for each variant.
#'
#' @examples
#' result <- run_mixada(pleio)
#' head(result)
#'
#' @export
run_mixada = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::mixAda, Sigma = zscore_sigma)

  pleio_p = data.frame(mixada_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run PC-Fisher Pleiotropy Test
#'
#' Applies the `PCFisher` method to evaluate pleiotropy significance.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' This function extracts z-score matrices and applies the `PCFisher` method.
#'
#' @return A `data.frame` with `PCFisher`-based p-values for each variant.
#'
#' @examples
#' result <- run_pcfisher(pleio)
#' head(result)
#'
#' @export
run_pcfisher = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::PCFisher, Sigma = zscore_sigma)

  pleio_p = data.frame(pcfisher_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run PCLC Pleiotropy Test
#'
#' Applies the `PCLC` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `PCLC` method (Principal Component-based Likelihood Combination) combines p-values across multiple correlated phenotypes using principal component analysis.
#' It accounts for the correlation structure among phenotypes and provides an aggregated p-value for each variant.
#'
#' @return A `data.frame` with `PCLC`-based p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_pclc(pleio)
#' head(result)
#'
#' @export
run_pclc = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::PCLC, Sigma = zscore_sigma)

  pleio_p = data.frame(pclc_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run PCMinP Pleiotropy Test
#'
#' Applies the `PCMinP` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `PCMinP` method (Principal Component Minimum P-value) uses principal component analysis to combine association signals across multiple correlated traits.
#' It identifies the minimum p-value across principal components and adjusts for correlation.
#'
#' @return A `data.frame` with `PCMinP`-based p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_pcminp(pleio)
#' head(result)
#'
#' @export
run_pcminp = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::PCMinP, Sigma = zscore_sigma)

  pleio_p = data.frame(pcminp_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run SUM Pleiotropy Test
#'
#' Computes the sum of association signals across multiple phenotypes.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts effect sizes and standard errors, computes z-scores, and applies the `SUM` method.
#'
#' @return A `data.frame` with `SUM`-based p-values for each variant.
#'
#' @examples
#' result <- run_sum(pleio)
#' head(result)
#'
#' @export
run_sum = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::SUM, Sigma = zscore_sigma)

  pleio_p = data.frame(sum_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run VC Pleiotropy Test
#'
#' Applies the `VC` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `VC` (Variance Component) method models genetic effects as a variance component across multiple correlated traits.
#' It estimates the shared genetic contribution across phenotypes and tests for significant pleiotropy.
#'
#' @return A `data.frame` with `VC`-based p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_vc(pleio)
#' head(result)
#'
#' @export
run_vc = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::VC, Sigma = zscore_sigma)

  pleio_p = data.frame(vc_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run Wald Pleiotropy Test
#'
#' Applies the `Wald` test for pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts z-scores and applies the `Wald` method.
#'
#' @return A `data.frame` with `Wald`-based p-values for each variant.
#'
#' @examples
#' result <- run_wald(pleio)
#' head(result)
#'
#' @export
run_wald = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::Wald, Sigma = zscore_sigma)

  pleio_p = data.frame(wald_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run WI Pleiotropy Test
#'
#' Applies the `WI` method from the `MPAT` package to assess pleiotropy.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The `WI` (Weighted Integration) method combines association signals across multiple correlated traits using a weighted approach.
#' It accounts for correlation and assigns weights to different phenotypic associations.
#'
#' @return A `data.frame` with `WI`-based p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_wi(pleio)
#' head(result)
#'
#' @export
run_wi = function(pleio){
  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::WI, Sigma = zscore_sigma)

  pleio_p = data.frame(wi_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' Run PCO Pleiotropy Test
#'
#' Applies the `PCO` method to test pleiotropy across correlated phenotypes.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' Computes correlation structure and applies `PCO` for pleiotropy testing.
#'
#' @return A `data.frame` with `PCO`-based p-values for each variant.
#'
#' @examples
#' result <- run_pco(pleio)
#' head(result)
#'
#' @export
run_pco = function(pleio){

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)
  sigma_estimate = SigmaOEstimate(zscore_sigma, simNum = 1000)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::PCO, Sigma = zscore_sigma, sigma_estimate)

  pleio_p = data.frame(pleio_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

#' @export
run_pcaq = function(pleio){

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)
  sigma_estimate = SigmaXEstimate(zscore_sigma, simNum = 1000)

  pleio_p = apply(zscore_matrix, MARGIN = 1, MPAT::PCAQ, Sigma = zscore_sigma, sigma_estimate)

  pleio_p = data.frame(pleio_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)
}

