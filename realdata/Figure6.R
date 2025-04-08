library(ggplot2)
library(dplyr)
library(patchwork)


setwd("~/gNODE/realdata/data_cdiff/trajectory/comparison")
base_path1 <- "~/gNODE/realdata/data_cdiff/trajectory/"

Group <- c("Control","Case")
base_folders <- c("NODE", "pNODE", "gNODE")

path <- paste0(base_path1, "comparison")
setwd(path)
all_rrmse <- readRDS("all_rrmse.rds")
all_rrmse$Group <- factor(all_rrmse$Group, levels = c("Control", "Case"))
all_rrmse$Method <- factor(all_rrmse$Method, levels = c("NODE", "pNODE","gNODE"))

# Calculate summary statistics (mean ± SD)
summary_data <- all_rrmse %>%
  group_by(Group, Method) %>%
  summarise(
    mean_rrmse = mean(rrmse, na.rm = TRUE),
    se_rrmse = sd(rrmse, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  )

# Create the plot
plot2 <- ggplot(summary_data, aes(x = Group, y = mean_rrmse, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.8),  
    width = 0.7, 
    color = NA,  
    linewidth = 0) +  
  geom_errorbar(
    aes(ymin = mean_rrmse - se_rrmse, ymax = mean_rrmse + se_rrmse),
    position = position_dodge(width = 0.8),  
    width = 0.2,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c("#F9C89B","#7CA3B8","#CE8A8D")) +
  xlab("Group") +
  ylab("Relative RMSE") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +  
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    text = element_text(size=28),
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 30, margin = margin(r = 15, unit = "pt")), 
    axis.text.x = element_text(size = 28, color = "black"),  
    axis.text.y = element_text(size = 28, color = "black"),
    panel.grid.major = element_line(colour=NA),
    panel.background = element_rect(fill="transparent", colour=NA),
    plot.background = element_rect(fill="transparent", colour=NA),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0.5, 'cm')
  ) +
  guides(fill = guide_legend(byrow = TRUE))

# Save the plot
ggsave(
  filename = "Figure6.png",
  plot = plot2,
  device = "png",
  dpi = 600,
  units = "in",
  width = 12,
  height = 8
)

