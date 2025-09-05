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
  "low_mean_quality" = "brown", "bad_alignment" = "grey", "multimatch" = "purple"
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

counts <- str_c(args$roots, ".counts.tsv") %>%
  set_names(root_names) %>%
  map(read_tsv) %>%
  bind_rows(.id = "_sample")

library_counts <- str_c(args$roots, ".library_counts.tsv") %>%
  set_names(root_names) %>%
  map(read_tsv) %>%
  bind_rows(.id = "_sample")

summary <- str_c(args$roots, ".summary.tsv") %>%
  set_names(root_names) %>%
  map(read_tsv) %>%
  bind_rows(.id = "_sample")

samples <- unique(counts$`_sample`)
regions <- str_remove(names(counts)[str_detect(names(counts), "_nearest")], "_nearest")

# Filtering Summary
p_filtering_abs <- filter(summary, group == "filtered" & metric != "total" | group == "unfiltered" & metric == "total") %>%
  mutate(metric = if_else(group == "filtered", metric, "unfiltered")) %>%
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
  mutate(metric = if_else(group == "filtered", metric, "unfiltered")) %>%
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

p_matches_abs <- filter(summary, group == "unfiltered", metric != "total") %>%
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

# Library completeness
completeness <- select(library_counts, `_sample`, group, combination_id, count) %>%
  count(`_sample`, combination_id, wt = count, name = "count") %>%
  left_join(lib, .,
          by = join_by(`_id` == combination_id)) %>%
  mutate(`_sample` = factor(`_sample`, levels = samples)) %>%
  complete(`_sample`, nesting(`_id`, !!!rlang::syms(regions))) %>%
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
           cumprop_reads = cumsum(prop_reads)) %>%
    ungroup()
}

p_rep_all <- count(library_counts, `_sample`, !!!rlang::syms(regions), combination_status, combinations_in_library,
      combination_id, wt = count, name = "count") %>%
  calc_library_representation(`_sample`) %>%
  ggplot(aes(x = prop_rank, y = cumprop_reads, colour = `_sample`)) +
  geom_line() +
  labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "All Combinations") +
  scale_colour_brewer(palette = "Dark2", name = "") +
  guides(colour = guide_legend(direction = "horizontal"))

p_rep_lib <- filter(library_counts, combinations_in_library == 1) %>%
  count(`_sample`, !!!rlang::syms(regions), combination_status, combinations_in_library,
                   combination_id, wt = count, name = "count") %>%
  calc_library_representation(`_sample`) %>%
  ggplot(aes(x = prop_rank, y = cumprop_reads, colour = `_sample`)) +
  geom_line() +
  labs(x = "Proportional Rank", y = "Cumulative Proportion of Reads", title = "In Library") +
  scale_colour_brewer(palette = "Dark2", name = "") +
  guides(colour = guide_legend(direction = "horizontal"))

p_rep <- p_rep_all + p_rep_lib + guide_area() +
  plot_layout(heights = c(1, 0.2), widths = c(0.5, 0.5), design = "12\n33", guides = "collect")

# Assemble overall plot
p <- wrap_plots(p_filtering, p_matches, p_completeness, p_rep,
                heights = c(1, 1, 1, 1), widths = c(1)) +
  plot_annotation(title = "Read Counts QC")

ggsave(str_c(args$out, "/count_qc.pdf"), units = "cm", height = 40, width = 20)
