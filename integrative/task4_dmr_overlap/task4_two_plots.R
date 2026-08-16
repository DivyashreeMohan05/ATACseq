library(ggplot2)

# ============================================================
# Integrative Task 4
# Plot 1: DMR overlap with heterochromatin
# Plot 2: DMR overlap with differential ATAC peaks
# ============================================================

total_dmrs <- 3254

# ============================================================
# PLOT 1 — DMR overlap with ChromHMM Het states
# ============================================================

het_data <- data.frame(
  Het_category = c(
    "Kidney Het only",
    "Liver Het only",
    "Het in both tissues"
  ),
  count = c(
    87,   # 128 kidney Het minus 41 shared
    201,  # 242 liver Het minus 41 shared
    41
  )
)

het_data$percentage <- 100 * het_data$count / total_dmrs

het_data$label <- paste0(
  het_data$count,
  "\n(",
  sprintf("%.2f", het_data$percentage),
  "%)"
)

het_data$Het_category <- factor(
  het_data$Het_category,
  levels = c(
    "Kidney Het only",
    "Liver Het only",
    "Het in both tissues"
  )
)

plot_het <- ggplot(
  het_data,
  aes(
    x = Het_category,
    y = count,
    fill = Het_category
  )
) +
  geom_col(
    width = 0.65,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 4.5,
    lineheight = 1
  ) +
  scale_y_continuous(
    limits = c(0, 235),
    breaks = seq(0, 200, 50),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "DMRs overlapping heterochromatic states",
    subtitle = paste0(
      "329 of 3,254 significant DMRs (10.11%) ",
      "overlapped Het in at least one tissue"
    ),
    x = NULL,
    y = "Number of DMRs",
    caption = paste0(
      "Het includes the ChromHMM states Het_P and Het_S. ",
      "Categories are mutually exclusive."
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 12.5,
      lineheight = 1.05,
      margin = margin(b = 12)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      size = 11.5,
      lineheight = 1.05
    ),
    axis.text.y = element_text(size = 11),
    plot.caption = element_text(
      size = 9.5,
      hjust = 0,
      margin = margin(t = 10)
    ),
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 15,
      l = 15
    )
  )

# ============================================================
# PLOT 2 — DMR overlap with tissue-specific ATAC peaks
# ============================================================

atac_data <- data.frame(
  ATAC_peak_type = c(
    "Kidney-specific",
    "Kidney-specific",
    "Liver-specific",
    "Liver-specific"
  ),
  DMR_direction = c(
    "Kidney hypermethylated",
    "Liver hypermethylated",
    "Kidney hypermethylated",
    "Liver hypermethylated"
  ),
  count = c(
    0,
    76,
    459,
    3
  )
)

atac_data$ATAC_peak_type <- factor(
  atac_data$ATAC_peak_type,
  levels = c(
    "Kidney-specific",
    "Liver-specific"
  )
)

atac_data$DMR_direction <- factor(
  atac_data$DMR_direction,
  levels = c(
    "Kidney hypermethylated",
    "Liver hypermethylated"
  )
)

# Percentages within each ATAC peak group
atac_data$percentage <- ave(
  atac_data$count,
  atac_data$ATAC_peak_type,
  FUN = function(x) 100 * x / sum(x)
)

# Labels for the main segments
atac_data$label <- ifelse(
  atac_data$count > 0,
  paste0(
    atac_data$count,
    "\n(",
    sprintf("%.1f", atac_data$percentage),
    "%)"
  ),
  ""
)

# Replace 100.0% with 100%
atac_data$label[
  atac_data$percentage == 100
] <- paste0(
  atac_data$count[
    atac_data$percentage == 100
  ],
  "\n(100%)"
)

# Data for the two large labels
large_atac_labels <- atac_data[
  atac_data$count >= 10,
]

# Manually define their vertical positions
large_atac_labels$label_y <- c(
  38,      # 76-DMR segment
  232.5    # centre of the 459-DMR segment above the 3-DMR segment
)

# Small exceptional group
small_atac_label <- atac_data[
  atac_data$ATAC_peak_type == "Liver-specific" &
    atac_data$DMR_direction == "Liver hypermethylated",
]

plot_atac <- ggplot(
  atac_data,
  aes(
    x = ATAC_peak_type,
    y = count,
    fill = DMR_direction
  )
) +
  geom_col(
    width = 0.65,
    color = "white",
    linewidth = 0.4
  ) +

  # Labels for 76 and 459
  geom_text(
    data = large_atac_labels,
    aes(
      y = label_y,
      label = label
    ),
    size = 4.5,
    lineheight = 1
  ) +

  # Label for the small group of 3 DMRs
  geom_segment(
    data = small_atac_label,
    aes(
      x = ATAC_peak_type,
      xend = ATAC_peak_type,
      y = count,
      yend = 28
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = small_atac_label,
    aes(
      x = ATAC_peak_type,
      y = 38,
      label = "3\n(0.6%)"
    ),
    inherit.aes = FALSE,
    size = 4,
    lineheight = 1
  ) +

  scale_fill_manual(
    values = c(
      "Kidney hypermethylated" = "#F8766D",
      "Liver hypermethylated" = "#00BFC4"
    )
  ) +

  scale_x_discrete(
    labels = c(
      "Kidney-specific" =
        "Kidney-specific\nATAC peaks",
      "Liver-specific" =
        "Liver-specific\nATAC peaks"
    )
  ) +

  scale_y_continuous(
    limits = c(0, 520),
    breaks = seq(0, 500, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +

  labs(
    title = "DMRs overlapping differential ATAC peaks",
    subtitle = paste0(
      "538 of 3,254 significant DMRs (16.53%) ",
      "overlapped tissue-specific accessible peaks"
    ),
    x = NULL,
    y = "Number of DMRs",
    fill = "DMR direction",
    caption = paste0(
      "535 of 538 overlaps (99.44%) showed an inverse ",
      "methylation–accessibility pattern"
    )
  ) +

  coord_cartesian(clip = "off") +

  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 12.5,
      lineheight = 1.05,
      margin = margin(b = 12)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      size = 11.5,
      lineheight = 1.05
    ),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 10.5),
    plot.caption = element_text(
      size = 9.5,
      hjust = 0,
      margin = margin(t = 10)
    ),
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 15,
      l = 15
    )
  )

# ============================================================
# Display plots
# ============================================================

print(plot_het)
print(plot_atac)

# ============================================================
# Save individual plots
# ============================================================

ggsave(
  filename = "task4_plot1_DMR_Het_overlap.png",
  plot = plot_het,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "task4_plot1_DMR_Het_overlap.pdf",
  plot = plot_het,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "task4_plot2_DMR_ATAC_overlap.png",
  plot = plot_atac,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "task4_plot2_DMR_ATAC_overlap.pdf",
  plot = plot_atac,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

# ============================================================
# Optional: combine both plots into one presentation figure
# ============================================================

if (requireNamespace("patchwork", quietly = TRUE)) {

  combined_plot <- (
    plot_het + plot_atac
  ) +
    patchwork::plot_layout(
      widths = c(1, 1.15)
    ) +
    patchwork::plot_annotation(
      title = paste0(
        "Kidney–liver DMR overlap with chromatin states ",
        "and accessibility"
      ),
      theme = theme(
        plot.title = element_text(
          size = 20,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 12)
        )
      )
    )

  print(combined_plot)

  ggsave(
    filename = "task4_two_panel_figure.png",
    plot = combined_plot,
    width = 15,
    height = 7,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = "task4_two_panel_figure.pdf",
    plot = combined_plot,
    width = 15,
    height = 7,
    units = "in",
    bg = "white"
  )

  cat("\nCombined two-panel figure created.\n")

} else {

  cat(
    "\nThe two individual plots were created successfully.\n",
    "The package 'patchwork' was not available, so the ",
    "combined figure was not created.\n",
    sep = ""
  )
}

cat("\nFiles created:\n")
cat("task4_plot1_DMR_Het_overlap.png\n")
cat("task4_plot1_DMR_Het_overlap.pdf\n")
cat("task4_plot2_DMR_ATAC_overlap.png\n")
cat("task4_plot2_DMR_ATAC_overlap.pdf\n")
