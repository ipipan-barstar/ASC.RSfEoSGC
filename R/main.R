
RESULTS_DIR <- "../Results/"
HASHTAGS_10 <- c("#1", "#90dayfiance", "#anjisalvacion", "#bbnaija", "#lolinginlove", "#nowplaying", "#puredoctrinesofchrist", "#tejasswiprakash", "#tejran", "#ukraine")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (grepl("-", args[1])) {
    SELECTION_NAME <- args[1]
    selection <- as.numeric(strsplit(SELECTION_NAME, "-")[[1]])
    TAG_COUNT <- length(selection)
  } else {
    TAG_COUNT <- as.integer(args[1])
    if (is.na(TAG_COUNT) | TAG_COUNT < 3 | TAG_COUNT > 10) {
      stop(paste("The provided integer is invalid:", TAG_COUNT))
    }
  }
} else {
  TAG_COUNT <- 3
}

ntries <- choose(10, TAG_COUNT)
if (!exists("SELECTION_NAME")) {
  repeat {
      if (ntries < 1) {
        stop(paste("Results for all subsets of length ", TAG_COUNT, " are calculated"))
      } else {
        ntries <- ntries - 1
      }
    selection <- sort(sample(0:9, TAG_COUNT, replace = FALSE))
    SELECTION_NAME <- paste(selection, collapse = "-")
    FIG_DIR <- paste0(RESULTS_DIR, SELECTION_NAME, ".figs/")
    if (!file.exists(FIG_DIR)) {
      break
    } else {
      cat(FIG_DIR, " exists, trying again \n")
    }
  }
} else {
  FIG_DIR <- paste0(RESULTS_DIR, SELECTION_NAME, ".figs/")
}

TARGET_DIR_PREFIX <- paste0(RESULTS_DIR, SELECTION_NAME, ".")
SELECTED_HASHTAGS <- HASHTAGS_10[selection + 1]

cat("Selected ids[", sep = "", TAG_COUNT, "]: ", paste(selection, collapse = " "), '\n')
cat("Selected hashtags:", SELECTED_HASHTAGS, '\n')
cat("Saving results to: ", TARGET_DIR_PREFIX, sep = "", "*\n")

# Function to pause and wait for the user to press Enter
execute <- function(proc_fname) {
  fancyline <- "  *********************  "
  invisible(readline(prompt = cat(fancyline, "Press [Enter] to execute ", proc_fname, fancyline)))
  source(proc_fname)
  paste(fancyline, "DONE  ", proc_fname, fancyline) }

execute("TWT_read_exp.R");
execute("TWT_doc_matrix.R");
execute("S_describeBis.R");

S_BAK <- S
nnn_cats_BAK <- nnn_cats
RES0_BAK <- RES0

for (lb in c(0, 0.1, 0.2)) {
  if (lb > 0) {
    cat("REDUCTION OF S by ", lb);

    print(sum((aDVtop100 - aDVlow100) < lb))
    S = S_BAK[(aDVtop100 - aDVlow100) >= lb, (aDVtop100 - aDVlow100) >= lb]
    nnn_cats = nnn_cats_BAK[(aDVtop100 - aDVlow100) >= lb]; RES0 = RES0_BAK[(aDVtop100 - aDVlow100) >= lb,]
  } else {
    cat("S WITHOUT REDUCTION\n");
  }
  print(table(nnn_cats));
  execute("find_Rcut.R");
  execute("find_NRcut.R");
  execute("NNN_spectral_clusterings.R");
}
