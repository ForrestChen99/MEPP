list.files('250630_深度测序_原位转移/4.SNP_Analysis/1.single_sample')
getwd()

library(readxl)
library(dplyr)
library(stringr)
library(purrr)

# ── directory that holds the *.xlsx files ────────────────────────────────
in_target_window <- function(type_string) {
  if (is.na(type_string) || type_string == "") return(FALSE)
  tokens <- str_split(type_string, ";", simplify = TRUE)
  pos    <- as.numeric(str_extract(tokens, "^[0-9]+"))
  any(pos %in% 104:113, na.rm = TRUE) ###Trp53 target window
}

in_target_window_pten <- function(type_string) {
  if (is.na(type_string) || type_string == "") return(FALSE)
  tokens <- str_split(type_string, ";", simplify = TRUE)
  pos    <- as.numeric(str_extract(tokens, "^[0-9]+"))
  any(pos %in% 25:50, na.rm = TRUE) ###Pten target window
}

one_file_pct <- function(fn) {
  df <- read_xlsx(file.path(path, fn))
  n_total <- nrow(df)
  
  n_window <- df %>%
    rowwise() %>%
    filter(in_target_window(Type)) %>%
    ungroup() %>%
    nrow()
  
  tibble(
    sample      = tools::file_path_sans_ext(fn),
    total_rows  = n_total,
    window_rows = n_window,
    percent     = round(100 * n_window / n_total, 2)
  )
}


one_file_pct_pten <- function(fn) {
  df <- read_xlsx(file.path(path, fn))
  n_total <- nrow(df)
  
  n_window <- df %>%
    rowwise() %>%
    filter(in_target_window_pten(Type)) %>%
    ungroup() %>%
    nrow()
  
  tibble(
    sample      = tools::file_path_sans_ext(fn),
    total_rows  = n_total,
    window_rows = n_window,
    percent     = round(100 * n_window / n_total, 2)
  )
}

path <- "250812整理/trp53"

result_tbl <- list.files(path, pattern = "\\.xlsx$", full.names = FALSE) %>% 
  map_dfr(one_file_pct)

print(result_tbl)

###Pten 25-50 abnormal
path <- "250812整理/pten"

result_tbl_pten <- list.files(path, pattern = "\\.xlsx$", full.names = FALSE) %>% 
  map_dfr(one_file_pct_pten)

print(result_tbl_pten)

res = rbind(result_tbl,result_tbl_pten)
writexl::write_xlsx(res,'res_250812.xlsx')



````````Frameshift mutation````````````
# ---- packages ----
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tools)

# ---- helper: net frameshift from a "Type" string ----
# expects tokens like "106D;107I;..."  (pos + I/D; other ops ignored)
net_shift <- function(type_string){
  if (is.na(type_string) || type_string == "" || type_string == "-") return(0L)
  toks <- str_split(type_string, ";", simplify = TRUE)
  len  <- suppressWarnings(as.integer(str_extract(toks, "^[0-9]+")))
  op   <- str_extract(toks, "[A-Za-z]$")
  sum(ifelse(op == "I", len,
             ifelse(op == "D", -len, 0L)), na.rm = TRUE)
}

# ---- helper: summarise ONE data frame into frameshift metrics ----
frameshift_from_df <- function(df){
  # guess columns present in your sheet
  type_col  <- intersect(c("Type","CIGAR","Mutation","Edit"), names(df))[1]
  count_col <- intersect(c("AltDepth","Count","Reads","Readcount","n"), names(df))[1]
  if (is.na(type_col) || is.na(count_col))
    stop("Couldn't find 'Type' and/or read-count columns in this file.")
  
  tmp <- df %>%
    transmute(
      type       = .data[[type_col]],
      n_reads    = as.numeric(.data[[count_col]]),
      net_indel  = vapply(.data[[type_col]], net_shift, integer(1)),
      has_indel  = str_detect(.data[[type_col]], "I|D"),
      frameshift = has_indel & (abs(net_indel) %% 3 != 0)
    )
  
  frameshift_reads <- sum(tmp$n_reads[tmp$frameshift], na.rm = TRUE)
  indel_reads      <- sum(tmp$n_reads[tmp$has_indel],  na.rm = TRUE)
  total_reads      <- sum(tmp$n_reads,                 na.rm = TRUE)
  
  tibble(
    total_reads  = total_reads,
    indel_reads  = indel_reads,
    frameshift_reads = frameshift_reads,
    pct_frameshift_of_indels = ifelse(indel_reads > 0,
                                      100 * frameshift_reads / indel_reads, NA_real_),
    pct_frameshift_of_total  = ifelse(total_reads > 0,
                                      100 * frameshift_reads / total_reads, NA_real_)
  )
}

# ---- helper: process ONE .xlsx file (first sheet) ----
one_file <- function(path_xlsx){
  df <- read_xlsx(path_xlsx, sheet = 1)
  out <- frameshift_from_df(df)
  out$file  <- basename(path_xlsx)
  out$gene  <- ifelse(str_detect(tolower(out$file), "pten"), "pten", "trp53")
  out %>% relocate(file, gene)
}

# ================= RUN on a folder =================
# set to your folder (e.g., "trp53" or "pten")

folder <- "250812整理/trp53"

files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)

summary_tbl <- map_dfr(files, ~ tryCatch(one_file(.x),
                                         error = function(e) {
                                           message("Skip ", .x, " :: ", e$message)
                                           tibble(file = basename(.x), gene = NA,
                                                  total_reads = NA, indel_reads = NA,
                                                  frameshift_reads = NA,
                                                  pct_frameshift_of_indels = NA,
                                                  pct_frameshift_of_total = NA)
                                         }))

print(summary_tbl %>% arrange(file))

##pten
folder <- "250812整理/pten"

files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)

summary_tbl_pten <- map_dfr(files, ~ tryCatch(one_file(.x),
                                         error = function(e) {
                                           message("Skip ", .x, " :: ", e$message)
                                           tibble(file = basename(.x), gene = NA,
                                                  total_reads = NA, indel_reads = NA,
                                                  frameshift_reads = NA,
                                                  pct_frameshift_of_indels = NA,
                                                  pct_frameshift_of_total = NA)
                                         }))

print(summary_tbl_pten %>% arrange(file))

res = rbind(summary_tbl,summary_tbl_pten)



writexl::write_xlsx(res, "frameshift_summary.xlsx")
