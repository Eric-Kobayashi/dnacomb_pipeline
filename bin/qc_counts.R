#!/usr/bin/env Rscript
# Produce simple QC results for read count tables
library(argparser)
library(magrittr)
library(tidyverse)
library(ggpubr)
library(patchwork)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

category_colours <- c(
  "match" = "green", "exact_match" = "green", "nearest_match" = "darkgreen",
  "mismatch" = "orange", "nonmatch" = "red",
  "recombination" = "blue", "exact_recombination" = "blue", "nearest_recombination" = "darkblue",
  "multimatch" = "purple"
)

filter_colours = c(
  unfiltered = "#4daf4a", empty_read = "#e41a1c", short_read = "#ff7f00",
  bad_alignment = "#984ea3", low_mean_quality = "#377eb8"
)

### Parse args
p <- arg_parser("Generate read count QC plots", name = "qc_counts.R", hide.opts = TRUE) %>%
  add_argument("--roots", "Root paths", nargs = Inf) %>%
  add_argument("--library", "Library TSV") %>%
  add_argument("--out", "Output directory", default = ".")

args <- parse_args(p)

lib <- read_tsv(args$library)
if (!"_id" %in% names(lib)) {
  lib <- mutate(lib, `_id` = 1:n())
}

root_names <- basename(args$roots)

# Don't actually use counts for now - save mem by not importing
# counts <- str_c(args$roots, ".counts.tsv") %>%
#   set_names(root_names) %>%
#   map(read_tsv) %>%
#   bind_rows(.id = "_sample") %>%
#   mutate(combination_status = factor(combination_status, levels = rev(names(category_colours))))

library_counts <- str_c(args$roots, ".library_counts.tsv") %>%
  set_names(root_names) %>%
  map(read_tsv) %>%
  bind_rows(.id = "_sample") %>%
  mutate(combination_status = factor(combination_status, levels = rev(names(category_colours))))

summary <- str_c(args$roots, ".summary.tsv") %>%
  set_names(root_names) %>%
  map(read_tsv) %>%
  bind_rows(.id = "_sample")

samples <- unique(library_counts$`_sample`)
n_samples <- length(samples)
regions <- setdiff(names(library_counts), c("_sample", "group", "combination_status", "combinations_in_library", "combination_id", "count"))
lib_regions <- regions[regions %in% names(lib)]

# Filtering Summary
p_filtering_abs <- filter(summary, group == "filtered" & metric != "total" | group == "unfiltered" & metric == "total") %>%
  mutate(metric = if_else(group == "filtered", metric, "unfiltered"),
         metric = factor(metric, levels = rev(names(filter_colours)))) %>%
  ggplot(aes(x = `_sample`, y = count, fill = metric)) +
  geom_col(width = 0.5) +
  scale_fill_manual(name = "", values = filter_colours) +
  labs(x = "", y = "Reads") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

p_filtering_prop <- filter(summary, group == "filtered" & metric != "total" | group == "unfiltered" & metric == "total") %>%
  mutate(metric = if_else(group == "filtered", metric, "unfiltered"),
         metric = factor(metric, levels = rev(names(filter_colours)))) %>%
  ggplot(aes(x = `_sample`, y = overall_proportion, fill = metric)) +
  geom_col(width = 0.5) +
  scale_fill_manual(name = "", values = filter_colours) +
  labs(x = "", y = "Proportion of reads") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

p_filtering <- p_filtering_abs + p_filtering_prop + guide_area() +
  plot_layout(heights = c(1, 0.2), widths = c(0.5, 0.5), design = "12\n33", guides = "collect")
h_filtering <- min(0.5 * n_samples, 5)

p_matches_abs <- filter(summary, group == "unfiltered", metric != "total") %>%
  mutate(metric = factor(metric, levels = rev(names(category_colours)))) %>%
  ggplot(aes(x = `_sample`, y = count, fill = metric)) +
  geom_col(width = 0.5) +
  scale_fill_manual(name = "", values = category_colours) +
  labs(x = "", y = "Reads") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

p_matches_prop <- filter(summary, group == "unfiltered", metric != "total") %>%
  mutate(metric = factor(metric, levels = rev(names(category_colours)))) %>%
  ggplot(aes(x = `_sample`, y = group_proportion, fill = metric)) +
  geom_col(width = 0.5) +
  scale_fill_manual(name = "", values = category_colours) +
  labs(x = "", y = "Proportion of reads") +
  coord_flip() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

p_matches <- p_matches_abs + p_matches_prop + guide_area() +
  plot_layout(heights = c(1, 0.2), widths = c(0.5, 0.5), design = "12\n33", guides = "collect")
h_matches <- min(0.5 * n_samples, 5)

# Library completeness
completeness <- select(library_counts, `_sample`, group, combination_id, count) %>%
  count(`_sample`, combination_id, wt = count, name = "count") %>%
  left_join(lib, ., by = join_by(`_id` == combination_id), relationship = "many-to-many") %>%
  mutate(`_sample` = factor(`_sample`, levels = samples)) %>%
  complete(`_sample`, nesting(`_id`, !!!rlang::syms(lib_regions))) %>%
  drop_na(`_sample`) %>%
  group_by(`_sample`) %>%
  summarise(p = sum(!is.na(count)) / n())

p_completeness <- ggplot(completeness, aes(y = `_sample`, x = p)) +
  geom_col(width = 0.5) +
  lims(x = c(0, 1)) +
  labs(x = "Completeness", y = "") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        axis.ticks.y = element_blank())
h_completeness <- min(0.5 * n_samples, 5)

# Library representation
calc_library_representation <- function(tbl, ...) {
  if (nrow(tbl) == 0) {
    return(mutate(tbl, rank = numeric(), prop_reads = numeric(),
                  prop_rank = numeric(), cumprop_reads = numeric()))
  }

  group_by(tbl, ...) %>%
    arrange(count) %>%
    mutate(rank = 1:n(),
           prop_reads = count/sum(count),
           prop_rank = rank/max(rank),
           cumprop_reads = cumsum(prop_reads),
           auc = sum(diff(prop_rank) * (head(cumprop_reads,-1) + tail(cumprop_reads,-1)))/2) %>%
    ungroup()
}

rep_all <- count(library_counts, `_sample`, !!!rlang::syms(regions), combination_status, combinations_in_library,
                 combination_id, wt = count, name = "count") %>%
  calc_library_representation(`_sample`)

rep_lib <- filter(library_counts, combinations_in_library == 1) %>%
  count(`_sample`, !!!rlang::syms(regions), combination_status, combinations_in_library,
        combination_id, wt = count, name = "count") %>%
  calc_library_representation(`_sample`)

if (n_samples > 8) {
  top_lib <- distinct(rep_all, `_sample`, auc) %>%
    arrange(auc) %>%
    slice(c(1, 2, 3, 4, n() - 3, n() - 2, n() - 1, n())) %>%
    pull(`_sample`)

  p_rep_all <- ggplot(mapping = aes(x = prop_rank, y = cumprop_reads)) +
    geom_line(data = filter(rep_all, !`_sample` %in% top_lib), mapping = aes(group = `_sample`)) +
    geom_line(data = filter(rep_all, `_sample` %in% top_lib), mapping = aes(colour = `_sample`)) +
    labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "All Combinations") +
    scale_colour_brewer(palette = "Dark2", name = "") +
    guides(colour = guide_legend(direction = "horizontal", nrow = 3))

  p_rep_lib <- ggplot(mapping = aes(x = prop_rank, y = cumprop_reads)) +
    geom_line(data = filter(rep_lib, !`_sample` %in% top_lib), mapping = aes(group = `_sample`)) +
    geom_line(data = filter(rep_lib, `_sample` %in% top_lib), mapping = aes(colour = `_sample`)) +
    labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "In Library") +
    scale_colour_brewer(palette = "Dark2", name = "") +
    guides(colour = guide_legend(direction = "horizontal", nrow = 3))
} else {
  leg_rows <- case_when(
    n_samples > 6 ~ 3,
    n_samples > 3 ~ 2,
    TRUE ~ 1
  )

  p_rep_all <- ggplot(rep_all, aes(x = prop_rank, y = cumprop_reads, colour = `_sample`)) +
    geom_line() +
    labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "All Combinations") +
    scale_colour_brewer(palette = "Dark2", name = "") +
    guides(colour = guide_legend(direction = "horizontal", nrow = leg_rows))

  p_rep_lib <- ggplot(rep_lib, aes(x = prop_rank, y = cumprop_reads, colour = `_sample`)) +
    geom_line() +
    labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "In Library") +
    scale_colour_brewer(palette = "Dark2", name = "") +
    guides(colour = guide_legend(direction = "horizontal", nrow = leg_rows))
}

p_rep <- p_rep_all + p_rep_lib + guide_area() +
  plot_layout(heights = c(1, 0.2), widths = c(0.5, 0.5), design = "12\n33", guides = "collect", axis_titles = "collect")
h_rep <- 15

# Assemble overall plot
heights <- c(h_filtering, h_matches, h_completeness, h_rep)
p <- wrap_plots(p_filtering, p_matches, p_completeness, p_rep,
                heights = heights, widths = c(1)) +
  plot_annotation(title = "Read Counts QC")

ggsave(str_c(args$out, "/count_qc.pdf"), units = "cm", height = sum(heights), width = 25)
ggsave(str_c(args$out, "/count_qc.png"), units = "cm", height = sum(heights), width = 25)
