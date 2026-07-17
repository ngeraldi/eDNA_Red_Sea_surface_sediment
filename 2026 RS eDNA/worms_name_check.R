# worms_name_check.R
# Run from RStudio console (requires internet access).
# Checks whether every species name in the filtered dataset already matches
# the WoRMS accepted name — i.e. verifies the reference library update worked.
#
# Output: prints a mismatch table to the console and writes
#         worms_name_check_results.csv to the '2026 RS eDNA' folder.
# -------------------------------------------------------------------------

library(tidyverse)
library(here)
library(worrms)
library(data.table)

dir      <- paste0(here::here("2026 RS eDNA"), "/")
dat_name <- "RS_surface_filtered_data_all.csv"

# ---------- 1. Build species list (mirrors tidy_sp + taxa2) ----------

dat <- data.table::fread(paste0(dir, dat_name), sep = ",")

tl <- c("Superkingdom", "Kingdom", "Phylum", "Class", "Order",
        "Family", "Genus", "Species", "taxon_id")

taxa <- dat %>%
  filter(primer == "18smini" | primer == "euka02") %>%
  select(lineage) %>%
  filter(!duplicated(lineage)) %>%
  tidyr::separate(lineage, into = tl, sep = ";", remove = FALSE)

taxa2 <- taxa %>%
  filter(Species != "NA" & !is.na(Species))

cat(sprintf("Species to check: %d\n\n", nrow(taxa2)))

# ---------- 2. Query WoRMS in batches of 50 ----------

xx  <- taxa2$Species
z   <- 50
stt <- seq(1, length(xx), by = z)

worm_info_list <- vector("list", length(stt))
for (idx in seq_along(stt)) {
  i   <- stt[idx]
  j   <- if (i == last(stt)) length(xx) else i + (z - 1)
  g   <- xx[i:j]
  gg  <- worrms::wm_records_names(name = g, fuzzy = FALSE, marine_only = FALSE)
  gg  <- Filter(Negate(is.null), gg)   # drop NULL entries (names not found)
  worm_info_list[[idx]] <- if (length(gg) > 0) bind_rows(gg) else NULL
  cat(sprintf("  checked %d / %d\n", min(j, length(xx)), length(xx)))
}
worm_info <- bind_rows(worm_info_list)   # safe: ignores NULLs, aligns columns

worm_info <- worm_info %>%
  select(AphiaID, scientificname, valid_AphiaID, valid_name,
         kingdom, phylum, class, order, family, genus, isMarine, status) %>%
  rename(AphiaID_old    = AphiaID,
         species_input  = scientificname,
         AphiaID        = valid_AphiaID,
         species_worrms = valid_name)

# ---------- 3. Compare input name vs WoRMS accepted name ----------

check <- taxa2 %>%
  select(Species) %>%
  rename(species_input = Species) %>%
  left_join(worm_info, by = "species_input") %>%
  mutate(
    name_matches = trimws(tolower(species_input)) == trimws(tolower(coalesce(species_worrms, ""))),
    flag = case_when(
      is.na(species_worrms)  ~ "NOT FOUND IN WORMS",
      !name_matches          ~ "NAME MISMATCH",
      TRUE                   ~ "OK"
    )
  )

# ---------- 4. Report ----------

n_ok       <- sum(check$flag == "OK")
n_mismatch <- sum(check$flag == "NAME MISMATCH")
n_missing  <- sum(check$flag == "NOT FOUND IN WORMS")

cat(sprintf("\n=== WoRMS name check results ===\n"))
cat(sprintf("  Total checked    : %d\n", nrow(check)))
cat(sprintf("  Names match (OK) : %d\n", n_ok))
cat(sprintf("  Name mismatches  : %d\n", n_mismatch))
cat(sprintf("  Not found        : %d\n\n", n_missing))

if (n_mismatch > 0) {
  cat("--- NAME MISMATCHES (input vs WoRMS accepted) ---\n")
  check %>%
    filter(flag == "NAME MISMATCH") %>%
    select(species_input, species_worrms, status) %>%
    print(n = Inf)
  cat("\n")
}

if (n_missing > 0) {
  cat("--- NOT FOUND IN WORMS ---\n")
  check %>%
    filter(flag == "NOT FOUND IN WORMS") %>%
    select(species_input) %>%
    print(n = Inf)
  cat("\n")
}

# ---------- 5. Save full results ----------

data.table::fwrite(check, paste0(dir, "worms_name_check_results.csv"), sep = ",")
cat(sprintf("Full results written to: %sworms_name_check_results.csv\n", dir))
