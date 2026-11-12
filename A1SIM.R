# Appendix Setting A1: prediction-weighted trace criteria T1-T4 (Figure S3TRACE).
# Use the shared data generator and the four trace calculations in ossrmp_preS3.R.
# The CPU column records elapsed fitting time in seconds for all four criteria.

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(ggplot2)
})

a1_project_dir <- getOption("A1TRACE.project_dir", "I:/softmax/OSSRMP")

# Helper functions

a1_trace_criteria <- function() paste0("T", 1:4)

make_a1_seed <- function(seed, offset = 0) {
  if (is.null(seed)) return(NULL)
  value <- as.double(seed) + offset
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value > .Machine$integer.max || value != floor(value)) stop("Invalid seed.")
  as.integer(value)
}

a1_distribution_args <- function(Dist) {
  switch(as.character(Dist),
    "1" = list(DistLabel = "D1-Gaussian-AR1", x_dist = "gaussian",
               cov_type = "ar1", rho = 0.5),
    "2" = list(DistLabel = "D2-t4-EXCH", x_dist = "t",
               cov_type = "exchangeable", rho = 0.5, df_t = 4),
    "3" = list(DistLabel = "D3-SkewQuad-AR1", x_dist = "skew_quadratic",
               cov_type = "ar1", rho = 0.5),
    "4" = list(DistLabel = "D4-GMix-EXCH", x_dist = "gaussian_mixture",
               cov_type = "exchangeable", rho = 0.5, mix_prob = 0.3),
    stop("Dist must be one of 1, 2, 3, 4."))
}

a1_method_files <- function(project_dir) {
  files <- file.path(project_dir, c("Syn.R", "ossrmp_preS3.R"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing source files: ", paste(missing, collapse = ", "))
  files
}

load_a1_methods <- function(files) {
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  methods <- new.env(parent = globalenv())
  suppressPackageStartupMessages(for (file in files) sys.source(file, methods))
  methods
}

a1_trace_mean <- function(x) if (length(x)) mean(x) else NA_real_
a1_trace_sd <- function(x) if (length(x) > 1L) sd(x) else NA_real_

a1_trace_rows <- function(i, Dist, cfg, values, cpu, errors) {
  criteria <- a1_trace_criteria()
  value <- unname(values[criteria])
  success <- is.finite(value)
  message <- unname(errors[criteria])
  message[!success & (is.na(message) | message == "")] <- "Non-finite trace returned."
  message[success] <- NA_character_
  data.frame(rep = i, Dist = Dist, DistLabel = a1_distribution_args(Dist)$DistLabel,
             K = cfg$K, N = cfg$N, n0 = cfg$n0, n_dag = cfg$n_dag, rho_N = cfg$rho_N,
             Criterion = criteria, Trace = ifelse(success, value, NA_real_),
             CPU = cpu, Success = success, ErrorMsg = message, stringsAsFactors = FALSE)
}

# One replication: generate one dataset per example and fit the four criteria.

fit_a1_trace <- function(i, Dist, cfg, methods) {
  criteria <- a1_trace_criteria()
  cpu <- NA_real_
  out <- tryCatch({
    args <- a1_distribution_args(Dist)
    args$DistLabel <- NULL
    # Use distinct seeds for data generation and fitting.
    seed_data <- make_a1_seed(cfg$seed_base, 100000 * i + 1000 * Dist)
    seed_fit <- make_a1_seed(seed_data, 1)
    # Rfast with seed = NULL does not follow set.seed(); use the seeded R path.
    dat <- do.call(methods$DataGenSM,
                   c(list(n = cfg$N, n_class = cfg$K + 1L, standardize_X = FALSE,
                          seed = seed_data, use_rfast = is.null(seed_data)), args))
    start <- proc.time()[["elapsed"]]
    values <- tryCatch(methods$ossrmp_preS3(
      X = dat$X, y = dat$y, n0 = cfg$n0, n_dag = cfg$n_dag, rho = cfg$rho_N,
      seed = seed_fit, tol = cfg$tol, maxit = cfg$maxit, verbose = cfg$verbose_fit,
      add_pilot_variance = cfg$add_pilot_variance,
      add_full_variance = cfg$add_full_variance, stable_bread = cfg$stable_bread),
      finally = cpu <- proc.time()[["elapsed"]] - start)
    if (!is.numeric(values) || length(values) != length(criteria) ||
        !setequal(names(values), criteria)) stop("Expected four named trace values, T1-T4.")
    errors <- attr(values, "error_messages")
    if (is.null(errors)) errors <- setNames(rep(NA_character_, length(criteria)), criteria)
    list(values = values, errors = errors)
  }, error = function(e) {
    list(values = setNames(rep(NA_real_, length(criteria)), criteria),
         errors = setNames(rep(conditionMessage(e), length(criteria)), criteria))
  })
  a1_trace_rows(i, Dist, cfg, out$values, cpu, out$errors)
}

run_a1_replication <- function(i, cfg, methods = a1_trace_worker_methods) {
  raw <- dplyr::bind_rows(lapply(cfg$Dist, function(Dist) fit_a1_trace(i, Dist, cfg, methods)))
  failures <- raw[!raw$Success, c("Dist", "Criterion", "ErrorMsg"), drop = FALSE]
  tag <- sprintf("K%d_rep%04d", cfg$K, i)
  saveRDS(raw, file.path(cfg$inner_out_dir, paste0("A1TRACE_", tag, ".rds")))
  if (nrow(failures))
    write.table(failures, file.path(cfg$inner_out_dir, paste0("errors_", tag, ".tsv")),
                sep = "\t", row.names = FALSE, quote = TRUE)
  raw
}

# Run independent replications; PSOCK workers load the same method files.

run_a1_repetitions <- function(cfg, m, cores, backend, files, methods) {
  workers <- min(as.integer(cores), as.integer(m))
  if (workers == 1L) {
    results <- lapply(seq_len(m), run_a1_replication, cfg = cfg, methods = methods)
  } else if (backend == "psock") {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(files, init) {
      assign("a1_trace_worker_methods", init(files), envir = .GlobalEnv)
      TRUE
    }, files, load_a1_methods)
    exports <- c("run_a1_replication", "fit_a1_trace", "a1_trace_rows",
                 "make_a1_seed", "a1_distribution_args", "a1_trace_criteria")
    parallel::clusterExport(cl, exports, envir = environment(run_a1_repetitions))
    results <- parallel::parLapply(cl, seq_len(m), run_a1_replication, cfg = cfg)
  } else {
    results <- parallel::mclapply(seq_len(m), run_a1_replication, cfg = cfg, methods = methods,
                                 mc.cores = workers, mc.set.seed = is.null(cfg$seed_base))
  }
  bad <- which(!vapply(results, is.data.frame, logical(1)))
  if (length(bad)) stop("Worker failed in replication(s): ", paste(bad, collapse = ", "))
  bind_rows(results)
}

summarize_a1_results <- function(raw) {
  raw %>%
    group_by(Dist, DistLabel, K, N, n0, n_dag, rho_N, Criterion) %>%
    summarise(B = sum(Success), Fail = sum(!Success),
              Mean = a1_trace_mean(Trace[Success]), SD = a1_trace_sd(Trace[Success]),
              SE = if (B > 0L) SD / sqrt(B) else NA_real_,
              CPU = a1_trace_mean(CPU[Success]), .groups = "drop") %>%
    arrange(K, Dist, Criterion)
}

# Run all K settings and save replicate, cumulative, and final results.

run_setting_a1 <- function(
    N = 1e6, K = c(3L, 5L, 7L), Dist = 1:4, m = 200L, n0 = 3000L, n_dag = 3000L,
    rho_N = 0.999, tol = 1e-8, maxit = 200L, verbose_fit = FALSE, cores = 4L,
    backend = c("auto", "psock", "multicore"), project_dir = a1_project_dir,
    out_dir = file.path(project_dir, "A1TRACE_runs"), seed_base = 12345L,
    add_pilot_variance = TRUE, add_full_variance = FALSE, stable_bread = FALSE) {
  is_count <- function(x) length(x) == 1L && is.finite(x) && x >= 1 &&
    x == floor(x) && x <= .Machine$integer.max
  stopifnot(is_count(N), is_count(n0), is_count(n_dag), is_count(m), is_count(cores),
            is_count(maxit), length(K) > 0L, all(vapply(K, is_count, logical(1))),
            !anyDuplicated(K), n0 > max(K), n0 < N, n_dag <= N - n0,
            length(Dist) > 0L, all(Dist %in% 1:4), !anyDuplicated(Dist),
            length(rho_N) == 1L, is.finite(rho_N), rho_N > 0, rho_N < 1,
            length(tol) == 1L, is.finite(tol), tol > 0)
  make_a1_seed(seed_base)
  make_a1_seed(seed_base, 100000 * m + 1000 * max(Dist) + 2 * length(K) + 1)
  backend <- match.arg(backend)
  if (backend == "auto") backend <- if (.Platform$OS.type == "windows") "psock" else "multicore"
  if (backend == "multicore" && .Platform$OS.type == "windows") stop("Use backend = 'psock' on Windows.")
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  files <- a1_method_files(project_dir)
  methods <- load_a1_methods(files)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  log <- function(...) {
    message <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    cat(message, "\n")
    cat(message, "\n", file = file.path(out_dir, "run.log"), append = TRUE)
  }
  cfg <- list(N = N, Dist = Dist, n0 = n0, n_dag = n_dag, rho_N = rho_N,
              tol = tol, maxit = maxit, verbose_fit = verbose_fit, seed_base = seed_base,
              add_pilot_variance = add_pilot_variance, add_full_variance = add_full_variance,
              stable_bread = stable_bread)
  metadata <- c(cfg, list(K = K, m = m, cores = cores, backend = backend,
                          source_md5 = tools::md5sum(files), session = sessionInfo()))
  saveRDS(metadata, file.path(out_dir, "A1TRACE_metadata.rds"))
  start <- Sys.time()
  results <- vector("list", length(K))
  for (j in seq_along(K)) {
    tag <- sprintf("K%d", K[j])
    task_cfg <- c(cfg, list(K = K[j], inner_out_dir = file.path(out_dir, "replicates", tag)))
    task_cfg$seed_base <- make_a1_seed(seed_base, 2 * j)
    dir.create(task_cfg$inner_out_dir, recursive = TRUE, showWarnings = FALSE)
    log("START ", tag, " | ", m, " replications")
    results[[j]] <- tryCatch(run_a1_repetitions(task_cfg, m, cores, backend, files, methods),
                            error = function(e) {
      log("ERROR ", tag, ": ", conditionMessage(e))
      structure(list(error = conditionMessage(e)), class = "sim_error")
    })
    saveRDS(results[[j]], file.path(out_dir, paste0("A1TRACE_", tag, ".rds")))
    ok <- vapply(results, is.data.frame, logical(1))
    if (any(ok)) saveRDS(bind_rows(results[ok]), file.path(out_dir, "A1TRACE_cumulative_raw.rds"))
    if (ok[j]) log("DONE ", tag, " | ", sum(results[[j]]$Success), "/", nrow(results[[j]]),
                   " finite traces")
  }
  ok <- vapply(results, is.data.frame, logical(1))
  if (!any(ok)) stop("All tasks failed; see run.log and task records.")
  raw <- bind_rows(results[ok])
  if (!any(raw$Success)) stop("No finite trace values; see error records in the replicates directory.")
  final <- list(raw = raw, summary = summarize_a1_results(raw), failed_K = K[!ok],
                elapsed_hours = as.numeric(difftime(Sys.time(), start, units = "hours")))
  saveRDS(final, file.path(out_dir, "A1TRACE_final.rds"))
  write.csv(final$summary, file.path(out_dir, "A1TRACE_summary.csv"), row.names = FALSE)
  log(sprintf("ALL DONE: %d/%d K tasks completed in %.2f hours.", sum(ok), length(K), final$elapsed_hours))
  final
}

# Boxplots for Figure S3TRACE.

plot_setting_a1 <- function(res_a1) {
  data <- res_a1$raw[res_a1$raw$Success, , drop = FALSE]
  if (!nrow(data)) stop("No finite traces to plot.")
  criteria <- a1_trace_criteria()
  data$Example <- factor(data$Dist, levels = 1:4, labels = paste("Example", 1:4))
  k_values <- sort(unique(data$K))
  data$KLabel <- factor(data$K, levels = k_values, labels = paste("K =", k_values))
  data$Criterion <- factor(data$Criterion, levels = criteria)
  ggplot2::ggplot(data, ggplot2::aes(Example, Trace, fill = Criterion)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge2(preserve = "single"),
                          width = 0.7, linewidth = 0.4, outlier.size = 0.7) +
    ggplot2::facet_wrap(~KLabel, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(T1 = "#1B9E77", T2 = "#D95F02",
                                        T3 = "#7570B3", T4 = "#E7298A"),
                               breaks = criteria, labels = expression(hat(T)[1], hat(T)[2],
                                                                      hat(T)[3], hat(T)[4])) +
    ggplot2::labs(x = NULL, y = "Trace value", fill = "Criterion") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))
}

# Simulation settings and calls
# Set options(A1TRACE.autorun = FALSE) before sourcing to load only the functions.

if (isTRUE(getOption("A1TRACE.autorun", TRUE))) {
  N <- 1e6
  K <- c(3, 5, 7)
  Dist <- 1:4
  m <- 200
  n0 <- 3000
  n_dag <- 3000
  rho_N <- 0.999
  cores <- 4L
  backend <- "auto"
  seed_base <- 12345
  project_dir <- a1_project_dir
  out_dir <- file.path(project_dir, "A1TRACE_runs")

  res_a1_trace <- run_setting_a1(
    N = N, K = K, Dist = Dist, m = m, n0 = n0, n_dag = n_dag, rho_N = rho_N,
    cores = cores, backend = backend, project_dir = project_dir, out_dir = out_dir,
    seed_base = seed_base)
  p_a1_trace <- plot_setting_a1(res_a1_trace)
  if (interactive()) print(p_a1_trace)
  ggplot2::ggsave(file.path(out_dir, "S3TRACE.pdf"), p_a1_trace, width = 13, height = 6)
}
