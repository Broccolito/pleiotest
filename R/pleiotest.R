pleiotest = function(pleio_object,
                     run_amatz_test = TRUE,
                     run_cmats_test = TRUE,
                     run_emats_test = TRUE,
                     run_dsum_test = TRUE,
                     run_metacca_test = TRUE,
                     run_minp_test = TRUE,
                     run_mixada_test = TRUE,
                     run_mixfisher_davies_test = TRUE,
                     run_mixfisher_liu_test = TRUE,
                     run_mixfisher_liumod_test = TRUE,
                     run_mixsd_davies_test = TRUE,
                     run_mixsd_liu_test = TRUE,
                     run_mixsd_liumod_test = TRUE,
                     run_mixtippett_davies_test = TRUE,
                     run_mixtippett_liu_test = TRUE,
                     run_mixtippett_liumod_test = TRUE,
                     run_mixvar_davies_test = TRUE,
                     run_mixvar_liu_test = TRUE,
                     run_mixvar_liumod_test = TRUE,
                     run_pcaq_test = FALSE,
                     run_pcfisher_test = TRUE,
                     run_pclc_test = TRUE,
                     run_pcminp_test = TRUE,
                     run_pco_test = FALSE,
                     run_sum_test = TRUE,
                     run_tates_test = TRUE,
                     run_vc_test = TRUE,
                     run_wald_test = TRUE,
                     run_wi_test = TRUE){

  pleioanalyze_result = tibble::tibble(snp_name = rownames(pleio_object@summary_stats_matrix[[1]]))
  execution_time = numeric()

  if(run_amatz_test){
    message("Running amatz test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_amatz(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["amatz"] = time_taken
  }

  if(run_cmats_test){
    message("Running cmats test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_cmats(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["cmats"] = time_taken
  }

  if(run_emats_test){
    message("Running emats test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_emats(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["emats"] = time_taken
  }

  if(run_dsum_test){
    message("Running dsum test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_dsum(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["dsum"] = time_taken
  }

  if(run_metacca_test){
    message("Running metacca test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_metacca(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["metacca"] = time_taken
  }

  if(run_minp_test){
    message("Running minp test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_minp(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["minp"] = time_taken
  }

  if(run_mixada_test){
    message("Running mixada test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixada(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixada"] = time_taken
  }

  if(run_mixfisher_liu_test){
    message("Running mixfisher_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixfisher_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixfisher_liu"] = time_taken
  }

  if(run_mixfisher_liumod_test){
    message("Running mixfisher_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixfisher_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixfisher_liumod"] = time_taken
  }

  if(run_mixsd_davies_test){
    message("Running mixsd_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_davies"] = time_taken
  }

  if(run_mixsd_liu_test){
    message("Running mixsd_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_liu"] = time_taken
  }

  if(run_mixsd_liumod_test){
    message("Running mixsd_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_liumod"] = time_taken
  }

  if(run_mixtippett_davies_test){
    message("Running mixtippett_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_davies"] = time_taken
  }

  if(run_mixtippett_liu_test){
    message("Running mixtippett_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_liu"] = time_taken
  }

  if(run_mixtippett_liumod_test){
    message("Running mixtippett_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_liumod"] = time_taken
  }

  if(run_mixvar_davies_test){
    message("Running mixvar_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_davies"] = time_taken
  }

  if(run_mixvar_liu_test){
    message("Running mixvar_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_liu"] = time_taken
  }

  if(run_mixvar_liumod_test){
    message("Running mixvar_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_liumod"] = time_taken
  }

  rownames(pleioanalyze_result) = NULL

  pleiot = new("pleiot",
               pleio_object = pleio_object,
               pleiotest_stats = pleioanalyze_result,
               execution_time = execution_time
  )

  return(pleiot)
}
