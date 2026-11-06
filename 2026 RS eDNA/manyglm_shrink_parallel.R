## manyglm_shrink_parallel.R
## Table 1 multivariate test (anova.manyglm) with cor.type = "shrink", bootstraps split
## across cores. anova.manyglm is single-threaded and "shrink" takes ~26 h per 999
## bootstraps on one M2 core (865 taxa), so the bootstraps are generated up front from a
## fixed seed (bootID = resampled sample ids), split into blocks, run in parallel and
## combined. mvabund's p = (1 + n bootstrap stats >= observed) / (1 + n bootstraps), so each
## block's exceedance count is p * (1 + n) - 1; counts are summed over blocks. The result
## is identical to one run with the same bootID (checked 2026-10-25), and reproducible.
##
## Input : manyglm_inputs.rds (export_mva_inputs chunk in Red_sea_analysis_plots_euk.Rmd)
## Output: manyglm_shrink_results.rds (read by table1_manyglm and sensitivity); saved after
##         each dataset, so a finished dataset is kept if the run is stopped. Datasets
##         already in the file are skipped (delete the file to re-run everything).
##
## Run from the terminal, not RStudio (forked workers are unreliable inside RStudio):
##   cd "<project>/2026 RS eDNA"
##   caffeinate -ims Rscript manyglm_shrink_parallel.R [n_cores] [n_boot] [datasets]
## defaults: 4 cores, 999 bootstraps, all datasets in the inputs file. Keys: pct01 = both
## primers, >= 0.1% presence (Table 1 since 2026-10-27); main = both primers, all
## detections (historic name); v9 / v7 = single primer, all detections; v9_pct01 /
## v7_pct01 = single primer, >= 0.1%. Progress -> manyglm_shrink_log.txt
## Time on an M2 with 4 cores: main 8.8 h, pct01 5.7 h, v9 1.3 h, v7 1.6 h (2026-10-26/27).

suppressMessages({
  library(mvabund)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
n_boot  <- if (length(args) >= 2) as.integer(args[2]) else 999
sets    <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else NULL   # NULL = all
seed    <- 20261025

## files sit next to this script
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))))
in_file  <- file.path(script_dir, "manyglm_inputs.rds")
out_file <- file.path(script_dir, "manyglm_shrink_results.rds")
log_file <- file.path(script_dir, "manyglm_shrink_log.txt")
if (n_boot < 99) out_file <- file.path(script_dir, "manyglm_shrink_results_TEST.rds")

logmsg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}

inp <- readRDS(in_file)
if (is.null(sets)) sets <- names(inp$comm)
env <- inp$env_mva
env$habitat <- relevel(factor(as.character(env$habitat)), ref = "Coral")
stopifnot(all(sets %in% names(inp$comm)),
          all(vapply(inp$comm, function(m) identical(rownames(m), as.character(env$sample_id)),
                     logical(1))))
rhs <- "habitat + sTemp + snitrate + spp_range + sland_distance"
n_samp <- nrow(env)

## bootstrap sample ids (pit.trap resamples whole samples with replacement)
set.seed(seed)
boot_id <- matrix(sample.int(n_samp, n_samp * n_boot, replace = TRUE),
                  nrow = n_boot)
storage.mode(boot_id) <- "integer"
## 4 blocks per core, handed out as cores free up (mc.preschedule = FALSE), so faster
## cores do more blocks (efficiency cores, thermal slow-down on a MacBook Air)
n_blocks <- min(n_boot, 4 * n_cores)
blocks <- split(seq_len(n_boot), cut(seq_len(n_boot), n_blocks, labels = FALSE))

res <- if (file.exists(out_file)) readRDS(out_file) else
  list(tables = list(), n_boot = n_boot, seed = seed, input_date = inp$date)
stopifnot(res$n_boot == n_boot, res$seed == seed)
todo <- setdiff(sets, names(res$tables))
logmsg("start: ", n_samp, " samples; ", n_boot, " bootstraps on ", n_cores, " cores; ",
       "to run: ", paste(todo, collapse = ", "),
       if (length(todo) < length(sets)) paste0(" (already done: ",
         paste(intersect(sets, names(res$tables)), collapse = ", "), ")"))

for (m in todo) {
  t0 <- Sys.time()
  Y  <- mvabund(as.matrix(inp$comm[[m]]))
  logmsg("dataset ", m, ": ", ncol(Y), " taxa; ", rhs)
  fm <- as.formula(paste("Y ~", rhs))

  ## shrink parameters are estimated by random cross-validation: same seed in every
  ## worker so all blocks use the same parameters
  set.seed(seed)
  mod <- manyglm(fm, data = env, family = binomial("cloglog"), cor.type = "shrink")

  run_block <- function(b) {
    set.seed(seed)
    a <- anova(mod, resamp = "pit.trap", test = "score", cor.type = "shrink",
               p.uni = "none", bootID = boot_id[blocks[[b]], , drop = FALSE],
               show.time = "none")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " dataset ", m, " block ", b, " of ",
        length(blocks), " done\n", file = log_file, append = TRUE, sep = "")
    a
  }
  out <- mclapply(seq_along(blocks), run_block, mc.cores = n_cores,
                  mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("block(s) failed: ", paste(which(bad), collapse = ", "), "\n",
                     paste(unlist(out[bad]), collapse = "\n"))

  ## combine: observed statistics are identical across blocks; sum exceedance counts
  ## (keep.boot's bootStat holds only one comparison, so counts come from the p-values)
  tb    <- as.data.frame(out[[1]]$table)
  n_b   <- vapply(out, function(a) a$n.bootsdone, numeric(1))
  count <- Reduce(`+`, lapply(seq_along(out), function(b)
    round(out[[b]]$table$`Pr(>score)` * (n_b[b] + 1) - 1)))
  tb$p  <- (1 + count) / (1 + sum(n_b))
  tb$`Pr(>score)` <- NULL
  res$n_boot_done[[m]] <- sum(n_b)
  res$tables[[m]] <- tb
  res$minutes[[m]] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(res, out_file)
  logmsg("dataset ", m, " done in ", round(res$minutes[[m]]), " min")
  print(tb)
}
logmsg("all done -> ", basename(out_file))
