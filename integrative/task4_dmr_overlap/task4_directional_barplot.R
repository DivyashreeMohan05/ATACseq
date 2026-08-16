library(ggplot2)

plot_data <- data.frame(
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
  count = c(0, 76, 459, 3)
)

# Percentage within each ATAC peak type
plot_data$percentage <- ave(
  plot_data$count,
  plot_data$ATAC_peak_type,
  FUN = function(x) 100 * x / sum(x)
)

plot_data$label <- ifelse(
  plot_data$count > 0,
  paste0(
    plot_data$count,
    "\n(",
    round(plot_data$percentage, 1),
    "%)"
  ),
  ""
)

p <- ggplot(
  plot_data,
  aes(
    x = ATAC_peak_type,
    y = percentage,
    fill = DMR_direction
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "DNA methylation is inversely related to chromatin accessibility",
    subtitle = "Kidney–liver DMRs overlapping tissue-specific ATAC-seq peaks",
    x = NULL,
    y = "Percentage of overlapping DMRs",
    fill = "DMR direction",
    caption = "535 of 538 overlaps (99.44%) show the expected inverse pattern"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  "task4_DMR_ATAC_direction_percentage.pdf",
  p,
  width = 8,
  height = 6
)

ggsave(
  "task4_DMR_ATAC_direction_percentage.png",
  p,
  width = 8,
  height = 6,
  dpi = 300
)

print(p)

